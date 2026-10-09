suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(shinycssloaders)
  library(cachem)
  library(DT)
})

# ============================================================
# Developing Mouse Visual Cortex Gene Expression Viewer
# ============================================================

app_title <- "Developing Mouse Visual Cortex Gene Expression Viewer"

s3_bucket_name <- "sea-ad-single-cell-profiling"
s3_region <- "us-west-2"
s3_prefix <- paste0(
  "Multiregion_2026/previous_objects/parquet/",
  "TEST_dev_mouse"
)

bucket <- arrow::s3_bucket(
  bucket = s3_bucket_name,
  region = s3_region,
  anonymous = TRUE
)

metadata_source <- bucket$path(
  file.path(s3_prefix, "metadata", "obs_metadata.parquet")
)
gene_map_source <- bucket$path(
  file.path(s3_prefix, "metadata", "gene_map.parquet")
)
cpm_source <- bucket$path(
  file.path(s3_prefix, "metadata", "CPM_scaling.csv")
)
statistics_source <- bucket$path(
  file.path(s3_prefix, "statistics", "mouse_gene_summary_statistics.csv")
)
value_sets_source <- bucket$path(
  file.path(s3_prefix, "metadata", "value_sets.csv")
)

read_startup_data <- function() {
  tryCatch(
    list(
      metadata = arrow::read_parquet(
        metadata_source,
        as_data_frame = TRUE
      ),
      gene_map = arrow::read_parquet(
        gene_map_source,
        as_data_frame = TRUE
      ),
      cpm = arrow::read_csv_arrow(
        cpm_source,
        as_data_frame = TRUE
      ),
      statistics = arrow::read_csv_arrow(
        statistics_source,
        as_data_frame = TRUE
      ),
      value_sets = arrow::read_csv_arrow(
        value_sets_source,
        as_data_frame = TRUE
      )
    ),
    error = function(e) {
      stop(
        paste(
          "Public AWS data connection failed.",
          "The application could not read the mouse visual cortex files.",
          "Please verify that the S3 prefix is public and available.",
          "Original error:",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )
}

startup <- read_startup_data()
metadata <- as.data.frame(startup$metadata, stringsAsFactors = FALSE)
gene_map <- as.data.frame(startup$gene_map, stringsAsFactors = FALSE)
cpm_scaling <- as.data.frame(startup$cpm, stringsAsFactors = FALSE)
gene_statistics <- as.data.frame(
  startup$statistics,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
value_sets <- as.data.frame(
  startup$value_sets,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
rm(startup)

table_column_definitions_path <- "table_column_definitions.csv"
if (!file.exists(table_column_definitions_path)) {
  stop(
    "Missing local table definition file: ",
    table_column_definitions_path
  )
}
table_column_definitions <- read.csv(
  table_column_definitions_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_definition_columns <- c(
  "column_names",
  "column_definitions"
)
missing_definition_columns <- setdiff(
  required_definition_columns,
  names(table_column_definitions)
)
if (length(missing_definition_columns) > 0) {
  stop(
    "table_column_definitions.csv is missing: ",
    paste(missing_definition_columns, collapse = ", ")
  )
}
table_column_definitions <- table_column_definitions |>
  transmute(
    column_names = as.character(column_names),
    column_definitions = as.character(column_definitions)
  ) |>
  filter(
    !is.na(column_names),
    nzchar(column_names),
    !is.na(column_definitions),
    nzchar(column_definitions)
  ) |>
  distinct(column_names, .keep_all = TRUE)

# ============================================================
# Precomputed ABC Atlas gene URLs
# ============================================================
abc_atlas_urls_path <- "abc_atlas_gene_urls.csv.gz"

if (!file.exists(abc_atlas_urls_path)) {
  stop(
    "Missing local ABC Atlas URL file: ",
    abc_atlas_urls_path,
    call. = FALSE
  )
}

abc_atlas_gene_urls <- read.csv(
  gzfile(abc_atlas_urls_path, open = "rt"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(abc_atlas_gene_urls) != 2L) {
  stop(
    abc_atlas_urls_path,
    " must contain exactly two columns: gene symbol and URL.",
    call. = FALSE
  )
}

names(abc_atlas_gene_urls) <- c(
  "gene_symbol",
  "abc_atlas_url"
)

abc_atlas_gene_urls$gene_symbol <- trimws(
  as.character(abc_atlas_gene_urls$gene_symbol)
)
abc_atlas_gene_urls$abc_atlas_url <- trimws(
  as.character(abc_atlas_gene_urls$abc_atlas_url)
)

valid_abc_atlas_rows <-
  !is.na(abc_atlas_gene_urls$gene_symbol) &
  nzchar(abc_atlas_gene_urls$gene_symbol) &
  !is.na(abc_atlas_gene_urls$abc_atlas_url) &
  nzchar(abc_atlas_gene_urls$abc_atlas_url) &
  startsWith(
    abc_atlas_gene_urls$abc_atlas_url,
    "https://knowledge.brain-map.org/abcatlas#"
  )

abc_atlas_gene_urls <- abc_atlas_gene_urls[
  valid_abc_atlas_rows,
  ,
  drop = FALSE
]
abc_atlas_gene_urls <- abc_atlas_gene_urls[
  !duplicated(abc_atlas_gene_urls$gene_symbol),
  ,
  drop = FALSE
]

if (nrow(abc_atlas_gene_urls) == 0L) {
  stop(
    abc_atlas_urls_path,
    " contains no valid ABC Atlas gene URLs.",
    call. = FALSE
  )
}

abc_atlas_url_lookup <- setNames(
  abc_atlas_gene_urls$abc_atlas_url,
  abc_atlas_gene_urls$gene_symbol
)

# ============================================================
# Optional searchable glossary
# ============================================================
tokens_path <- file.path("www", "tokens.csv")
tokens_available <- FALSE
tokens <- data.frame(
  token = character(),
  type = character(),
  name = character(),
  stringsAsFactors = FALSE
)
if (file.exists(tokens_path)) {
  tokens_candidate <- tryCatch(
    read.csv(
      tokens_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
  required_token_columns <- c("token", "type", "name")
  if (
    !is.null(tokens_candidate) &&
    all(required_token_columns %in% names(tokens_candidate))
  ) {
    tokens <- tokens_candidate |>
      transmute(
        token = trimws(as.character(token)),
        type = trimws(as.character(type)),
        name = trimws(as.character(name))
      ) |>
      filter(
        !is.na(token),
        nzchar(token),
        !is.na(type),
        nzchar(type),
        !is.na(name),
        nzchar(name)
      ) |>
      distinct(token, .keep_all = TRUE) |>
      arrange(tolower(token), token)
    tokens_available <- nrow(tokens) > 0
  }
}

# ============================================================
# Validate and align supporting data
# ============================================================

required_metadata <- c("sample_id", "library_label")
missing_metadata <- setdiff(required_metadata, names(metadata))
if (length(missing_metadata) > 0) {
  stop(
    "obs_metadata.parquet is missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}

required_gene_map <- c("gene_symbol", "gene_column")
missing_gene_map <- setdiff(required_gene_map, names(gene_map))
if (length(missing_gene_map) > 0) {
  stop(
    "gene_map.parquet is missing: ",
    paste(missing_gene_map, collapse = ", ")
  )
}

required_cpm <- c("library_label", "CPM_scaling_factor")
missing_cpm <- setdiff(required_cpm, names(cpm_scaling))
if (length(missing_cpm) > 0) {
  stop(
    "CPM_scaling.csv is missing: ",
    paste(missing_cpm, collapse = ", ")
  )
}

if (nrow(metadata) != nrow(cpm_scaling)) {
  stop("CPM_scaling.csv and obs_metadata.parquet have different row counts.")
}

if (!identical(
  as.character(metadata$library_label),
  as.character(cpm_scaling$library_label)
)) {
  stop("CPM_scaling.csv library order does not match obs_metadata.parquet.")
}

if (!"donor_age" %in% names(metadata)) {
  stop("obs_metadata.parquet must contain donor_age to derive stage.")
}

metadata <- metadata |>
  mutate(
    sample_id = as.character(sample_id),
    library_label = as.character(library_label),
    CPM_scaling_factor = as.numeric(cpm_scaling$CPM_scaling_factor),
    stage = case_when(
      grepl("^E[0-9]+(\\.[0-9]+)?$", trimws(as.character(donor_age))) ~ "embryonic",
      grepl("^P[0-9]+(\\.[0-9]+)?$", trimws(as.character(donor_age))) ~ "postnatal",
      TRUE ~ NA_character_
    )
  )

if (anyDuplicated(metadata$sample_id)) {
  stop("sample_id must be unique in obs_metadata.parquet.")
}

names(gene_statistics) <- make.unique(
  names(gene_statistics),
  sep = "_"
)

possible_statistics_gene_columns <- c(
  "gene_symbol", "Gene Symbol", "Gene", "gene", "symbol"
)
statistics_gene_column <- possible_statistics_gene_columns[
  possible_statistics_gene_columns %in% names(gene_statistics)
][1]

if (is.na(statistics_gene_column)) {
  stop(
    "mouse_gene_summary_statistics.csv must contain one of: ",
    paste(possible_statistics_gene_columns, collapse = ", ")
  )
}

gene_statistics[[statistics_gene_column]] <- as.character(
  gene_statistics[[statistics_gene_column]]
)

gene_map <- gene_map |>
  transmute(
    gene_symbol = as.character(gene_symbol),
    gene_column = as.character(gene_column)
  ) |>
  filter(
    !is.na(gene_symbol),
    nzchar(gene_symbol),
    !is.na(gene_column),
    nzchar(gene_column)
  ) |>
  distinct(gene_symbol, .keep_all = TRUE) |>
  arrange(gene_symbol)

available_genes <- gene_map$gene_symbol
gene_symbol_to_key <- setNames(
  gene_map$gene_column,
  gene_map$gene_symbol
)

# ============================================================
# Value sets: deterministic order and colors
# ============================================================

required_value_set_columns <- c(
  "label",
  "field",
  "color_hex_triplet",
  "order"
)
missing_value_set_columns <- setdiff(
  required_value_set_columns,
  names(value_sets)
)
if (length(missing_value_set_columns) > 0) {
  stop(
    "value_sets.csv is missing: ",
    paste(missing_value_set_columns, collapse = ", ")
  )
}

value_sets <- value_sets |>
  transmute(
    label = trimws(as.character(label)),
    field = trimws(as.character(field)),
    color = trimws(as.character(color_hex_triplet)),
    order = suppressWarnings(as.numeric(order))
  ) |>
  filter(
    !is.na(label),
    nzchar(label),
    !is.na(field),
    nzchar(field)
  ) |>
  arrange(field, order)

normalize_value_set_field <- function(field) {
  normalized <- tolower(trimws(as.character(field)))
  normalized <- gsub("[^a-z0-9]+", "_", normalized)
  normalized <- gsub("_label$", "", normalized)
  normalized <- gsub("_alias$", "", normalized)
  normalized <- gsub("^_|_$", "", normalized)
  
  aliases <- c(
    "brain_region" = "region_of_interest",
    "region" = "region_of_interest",
    "roi" = "region_of_interest",
    "age" = "donor_age",
    "age_label" = "donor_age",
    "developmental_age" = "donor_age",
    "donor_sex" = "sex"
  )
  
  if (normalized %in% names(aliases)) {
    unname(aliases[[normalized]])
  } else {
    normalized
  }
}

value_sets$field_key <- vapply(
  value_sets$field,
  normalize_value_set_field,
  character(1)
)

value_set_for <- function(field, observed_values) {
  observed <- unique(as.character(observed_values))
  observed <- observed[!is.na(observed) & nzchar(observed)]
  field_key <- normalize_value_set_field(field)
  
  specification <- value_sets |>
    filter(.data$field_key == field_key) |>
    arrange(order)
  
  if (nrow(specification) == 0) {
    return(NULL)
  }
  
  specification <- specification |>
    distinct(label, .keep_all = TRUE)
  
  defined_levels <- specification$label
  levels <- c(
    defined_levels[defined_levels %in% observed],
    observed[!observed %in% defined_levels]
  )
  
  defined_colors <- setNames(
    specification$color,
    specification$label
  )
  colors <- setNames(
    rep(NA_character_, length(levels)),
    levels
  )
  matched <- intersect(levels, names(defined_colors))
  colors[matched] <- defined_colors[matched]
  
  list(
    levels = levels,
    colors = colors
  )
}

# ============================================================
# Metadata fields and helpers
# ============================================================

first_existing <- function(candidates, data_names = names(metadata)) {
  hit <- candidates[candidates %in% data_names]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

age_field <- first_existing(c(
  "age", "age_label", "donor_age", "developmental_age",
  "age_in_days", "age_days"
))
region_field <- first_existing(c(
  "region_of_interest_label", "Brain Region", "brain_region",
  "region_label", "structure"
))
class_field <- first_existing(c("class", "class_label", "Class"))
subclass_field <- first_existing(c("subclass", "subclass_label", "Subclass"))
cluster_field <- first_existing(c("cluster_alias", "cluster", "cluster_label"))
subcluster_field <- first_existing(c(
  "subcluster_alias", "subcluster", "subcluster_label"
))

if (is.na(age_field)) {
  stop("No age field was found in obs_metadata.parquet.")
}
if (is.na(region_field)) {
  stop("No brain-region field was found in obs_metadata.parquet.")
}
if (is.na(cluster_field)) {
  stop("No cluster field was found in obs_metadata.parquet.")
}

preferred_fields <- unique(c(
  class_field,
  subclass_field,
  cluster_field,
  subcluster_field,
  region_field,
  age_field,
  "stage",
  "sex",
  "donor_sex",
  "library_label",
  "donor_label",
  "donor_id",
  "developmental_compartment",
  "sample_type"
))

preferred_fields <- preferred_fields[
  !is.na(preferred_fields) &
    preferred_fields %in% names(metadata)
]

categorical_fields <- preferred_fields[vapply(
  metadata[preferred_fields],
  function(x) {
    is.character(x) || is.factor(x) ||
      dplyr::n_distinct(x, na.rm = TRUE) <= 100
  },
  logical(1)
)]

omit_plot_fields <- c(
  "library_label"
)
omit_filter_fields <- c(
  "library_prep"
)

plot_fields <- setdiff(
  categorical_fields,
  omit_plot_fields
)
numeric_filter_fields <- names(metadata)[vapply(
  metadata,
  is.numeric,
  logical(1)
)]
numeric_filter_fields <- setdiff(
  numeric_filter_fields,
  c(
    "CPM_scaling_factor"
  )
)
filter_fields <- setdiff(
  unique(c(
    categorical_fields,
    numeric_filter_fields
  )),
  omit_filter_fields
)

filter_input_id <- function(field) {
  field_index <- match(field, filter_fields)
  if (is.na(field_index)) {
    stop("Unknown filter field: ", field)
  }
  paste0("stack_filter_", field_index)
}

apply_filter_specification <- function(data, filters) {
  if (length(filters) == 0) {
    return(data)
  }
  
  for (specification in filters) {
    field <- specification$field
    if (!field %in% names(data)) {
      next
    }
    
    if (identical(specification$type, "numeric")) {
      bounds <- as.numeric(specification$value)
      if (
        length(bounds) == 2 &&
        all(is.finite(bounds))
      ) {
        numeric_values <- suppressWarnings(
          as.numeric(data[[field]])
        )
        data <- data[
          is.finite(numeric_values) &
            numeric_values >= min(bounds) &
            numeric_values <= max(bounds),
          ,
          drop = FALSE
        ]
      }
    } else {
      included_values <- as.character(specification$value)
      if (length(included_values) > 0) {
        data <- data[
          as.character(data[[field]]) %in% included_values,
          ,
          drop = FALSE
        ]
      }
    }
  }
  
  data
}

cell_type_fields <- unique(na.omit(c(
  class_field,
  subclass_field,
  cluster_field,
  subcluster_field
)))

parse_age <- function(x) {
  text <- toupper(trimws(as.character(x)))
  out <- suppressWarnings(as.numeric(text))
  
  embryonic <- grepl("^E\\s*[0-9]", text)
  postnatal <- grepl("^P\\s*[0-9]", text)
  
  extract_number <- function(values) {
    suppressWarnings(as.numeric(
      sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", values)
    ))
  }
  
  out[embryonic] <- extract_number(text[embryonic]) - 100
  out[postnatal] <- extract_number(text[postnatal])
  out
}

capped_age_positions <- function(x, maximum_step = 3) {
  labels <- age_level_order(x)
  raw_positions <- parse_age(labels)
  
  valid <- is.finite(raw_positions)
  labels <- labels[valid]
  raw_positions <- raw_positions[valid]
  
  if (length(raw_positions) == 0) {
    return(setNames(numeric(), character()))
  }
  
  if (length(raw_positions) == 1) {
    return(setNames(0, labels))
  }
  
  raw_steps <- diff(raw_positions)
  capped_steps <- pmin(raw_steps, maximum_step)
  
  positions <- c(0, cumsum(capped_steps))
  setNames(positions, labels)
}

age_level_order <- function(x) {
  values <- unique(as.character(x))
  values <- values[!is.na(values) & nzchar(values)]
  numeric_values <- parse_age(values)
  values[order(numeric_values, values, na.last = TRUE)]
}

progression_positions <- function(data, field, maximum_step = 3) {
  values <- as.character(data[[field]])
  values <- values[!is.na(values) & nzchar(values)]
  labels <- field_levels(data, field)
  labels <- labels[labels %in% values]
  parsed <- parse_age(labels)
  recognized_age <- grepl(
    "^[EP]\\s*[0-9]",
    toupper(trimws(labels))
  )
  
  if (
    length(labels) > 0 &&
    all(recognized_age) &&
    all(is.finite(parsed))
  ) {
    raw_steps <- diff(parsed)
    positions <- if (length(raw_steps) == 0) {
      1
    } else {
      c(1, 1 + cumsum(pmin(raw_steps, maximum_step)))
    }
  } else {
    positions <- seq_along(labels)
  }
  
  names(positions) <- labels
  half_step <- if (length(positions) < 2) {
    0.5
  } else {
    max(0.5, min(diff(positions)) / 2)
  }
  
  list(
    labels = labels,
    breaks = unname(positions),
    lookup = positions,
    limits = c(
      min(positions) - half_step,
      max(positions) + half_step
    )
  )
}

field_levels <- function(data, field) {
  values <- unique(as.character(data[[field]]))
  values <- values[!is.na(values) & nzchar(values)]
  specification <- value_set_for(field, values)
  
  if (!is.null(specification)) {
    specification$levels
  } else if (identical(field, age_field)) {
    age_level_order(values)
  } else {
    sort(values)
  }
}

factor_field <- function(data, field) {
  data[[field]] <- factor(
    as.character(data[[field]]),
    levels = field_levels(data, field)
  )
  data
}

distinct_category_colors <- function(n) {
  if (n <= 0) {
    return(character())
  }
  grDevices::hcl.colors(
    n,
    palette = "Dark 3",
    alpha = 1,
    rev = FALSE
  )
}

field_colors <- function(data, field) {
  levels <- field_levels(data, field)
  specification <- value_set_for(field, data[[field]])
  
  if (!is.null(specification)) {
    colors <- specification$colors[levels]
    missing <- is.na(colors) | !grepl("^#[0-9A-Fa-f]{6}$", colors)
    
    if (any(missing)) {
      colors[missing] <- distinct_category_colors(sum(missing))
    }
    
    return(colors)
  }
  
  setNames(distinct_category_colors(length(levels)), levels)
}

manual_color_scale <- function(data, field, aesthetic = "color") {
  colors <- field_colors(data, field)
  if (identical(aesthetic, "fill")) {
    scale_fill_manual(
      values = colors,
      breaks = names(colors),
      drop = TRUE,
      na.value = "#808080"
    )
  } else {
    scale_color_manual(
      values = colors,
      breaks = names(colors),
      drop = TRUE,
      na.value = "#808080"
    )
  }
}

fields_within_limit <- function(data, fields, maximum) {
  fields[vapply(
    fields,
    function(field) {
      field %in% names(data) &&
        n_distinct(data[[field]], na.rm = TRUE) <= maximum
    },
    logical(1)
  )]
}

make_cache_key <- function(gene_symbol) {
  paste0(
    "gene",
    paste(as.integer(charToRaw(enc2utf8(gene_symbol))), collapse = "")
  )
}

gene_cache <- cachem::cache_mem(
  max_size = 1024 * 1024^2,
  max_age = Inf
)

default_gene <- if ("Reln" %in% available_genes) {
  "Reln"
} else {
  available_genes[[1]]
}
default_x <- if ("donor_age" %in% plot_fields) {
  "donor_age"
} else {
  age_field
}

default_second <- if (!is.na(subclass_field)) {
  subclass_field
} else {
  cluster_field
}
default_facet <- if (!is.na(subclass_field)) subclass_field else cluster_field
default_color <- region_field

# ============================================================
# Header and UI
# ============================================================

app_header <- div(
  class = "mouse-header",
  tags$a(
    class = "mouse-logo-link",
    href = "https://alleninstitute.org/",
    target = "_blank",
    tags$img(
      src = "allen_institute_logo.svg",
      alt = "alleninstitute.org",
      class = "mouse-logo"
    )
  ),
  tags$span(class = "mouse-title", app_title),
  div(
    class = "mouse-header-links",
    tags$a(
      href = "https://alleninstitute.github.io/abc_atlas_access/descriptions/Dev-Mouse-Vis-Cortex-dataset.html",
      target = "_blank",
      icon("brain"),
      tags$span("Dataset")
    ),
    tags$a(
      href = "https://brain-map.org",
      target = "_blank",
      icon("globe"),
      tags$span("Brain Map")
    )
  )
)

ui <- function(request) {
  page_sidebar(
    title = app_header,
    theme = bs_theme(
      version = 5,
      primary = "#214E68",
      success = "#4B9B58"
    ),
    
    tags$head(
      tags$title(app_title),
      uiOutput("plot_size_css"),
      tags$link(
        rel = "icon",
        type = "image/png",
        href = "allen_institute_logo.svg"
      ),
      tags$style(HTML(
        "
      html, body { width: 100%; overflow-x: hidden; }
      .navbar, header.navbar {
        background: #111827 !important;
        min-height: 64px !important;
        border: 0 !important;
      }
      .navbar-brand { width: 100% !important; padding: 0 !important; }
      .mouse-header {
        position: relative; display: flex; align-items: center;
        justify-content: space-between; width: 100%; height: 64px;
        padding: 0 16px;
      }
      .mouse-header-spacer {
        width: 210px;
      }
      
      .sidebar-logo-container {
        width: 100%;
        text-align: center;
        margin: 0 0 7px 0;
      }
      
      .sidebar-logo {
        display: block;
        width: 100%;
        max-width: 345px;
        max-height: 250px;
        object-fit: contain;
        margin: 0 auto;
      }
      .mouse-logo-link { width: 180px; display: flex; align-items: center; }
      .mouse-logo { height: 52px; max-width: 175px; object-fit: contain; }
      .mouse-title {
        position: absolute; left: 50%; transform: translateX(-50%);
        color: white; font-size: 21px; font-weight: 700; white-space: nowrap;
      }
      .mouse-header-links { display: flex; gap: 6px; width: 210px; justify-content: flex-end; }
      .mouse-header-links a {
        color: white !important; text-decoration: none; font-size: 13px;
        padding: 18px 8px;
      }
      .bslib-sidebar-layout > .sidebar,
      .bslib-sidebar-layout > aside {
        background: #214E68 !important; color: white !important;
        padding: 8px 10px !important;
      }
      .bslib-sidebar-layout > .sidebar label,
      .bslib-sidebar-layout > .sidebar h4,
      .bslib-sidebar-layout > .sidebar .form-check-label,
      .bslib-sidebar-layout > .sidebar .help-block {
        color: white !important;
      }
      .compact-controls h4, .gene-controls h4 {
        font-size: 1rem; margin: 4px 0;
      }
      .compact-controls hr { margin: 5px 0; border-color: rgba(255,255,255,.3); }
      .compact-controls .shiny-input-container { margin-bottom: 3px !important; }
      /* Standard input labels */
      .compact-controls .shiny-input-container > label {
        font-size: 0.75rem !important;
        line-height: 1.05 !important;
        margin-bottom: 1px !important;
      }
      /* Checkbox labels across Shiny and bslib versions */
      .compact-controls .form-check-label,
      .compact-controls .shiny-input-checkbox label,
      .compact-controls .checkbox label,
      .compact-controls input[type='checkbox'] + label {
        font-size: 0.75rem !important;
        line-height: 1.05 !important;
        margin-bottom: 1px !important;
        white-space: nowrap;
      }
      .compact-controls .form-check {
        min-height: 19px !important;
        margin: 1px 0 2px 0 !important;
      }
      .compact-controls .form-check-input {
        margin-top: 0.10rem !important;
      }
      .compact-controls .form-select,
      .compact-controls .form-control,
      .compact-controls .selectize-input {
        min-height: 28px !important; font-size: .78rem !important;
      }
      .repeat-plot-note {
        color: rgba(255, 255, 255, 0.85);
        font-size: 0.67rem;
        line-height: 1.1;
        text-align: center;
        margin: 3px 0 4px 0;
      }
      #filter_fields + .selectize-control .selectize-input,
      #filter_fields + .selectize-control .selectize-dropdown,
      #filter_fields-selectized + .selectize-dropdown,
      .stacked-filter-control .selectize-input,
      .stacked-filter-control .selectize-dropdown {
        background: white !important; color: #111827 !important;
      }
      .stacked-filter-panel {
        border: 1px solid rgba(255,255,255,.28);
        border-radius: 5px;
        padding: 5px 6px 2px;
        margin: 3px 0 5px;
        background: rgba(255,255,255,.06);
      }
      .stacked-filter-control {
        margin-bottom: 5px;
      }
      .stacked-filter-control:last-child {
        margin-bottom: 1px;
      }
      .stacked-filter-control .irs--shiny .irs-bar,
      .stacked-filter-control .irs--shiny .irs-single,
      .stacked-filter-control .irs--shiny .irs-from,
      .stacked-filter-control .irs--shiny .irs-to {
        background: #4B9B58;
        border-color: #4B9B58;
      }
      .gene-status { font-size: .76rem; min-height: 1.2rem; margin: 2px 0; }
      .btn-primary, .btn-outline-primary {
        background: #4B9B58 !important; border-color: white !important;
        color: white !important; font-weight: 600;
      }
      .filter-row { display: grid; grid-template-columns: 1fr 1fr; gap: 6px; }
      .gene-table-panel { padding: 8px; overflow-x: auto; }
      .statistics-value-filter-instructions {
        font-size: 1rem; line-height: 1.5; margin: 0 0 1rem;
        color: inherit;
      }
      .statistics-value-filter-panel {
        display: grid; grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: 8px; width: 100%; box-sizing: border-box;
      }
      .statistics-value-filter-panel .shiny-input-container {
        width: 100% !important; margin-bottom: 0 !important;
      }
      .statistics-value-filter-panel label {
        color: #003057; font-size: .78rem; font-weight: 600;
        margin-bottom: 2px;
      }
      @media (max-width: 1100px) {
        .statistics-value-filter-panel {
          grid-template-columns: repeat(2, minmax(0, 1fr));
        }
      }
      .gene-table-panel table.dataTable { font-size: .82rem; }
      .gene-table-panel table.dataTable thead th.has-definition {
        cursor: help;
        text-decoration: underline dotted;
        text-underline-offset: 2px;
      }
      label.dropdown-help {
        cursor: help;
        text-decoration: underline dotted;
        text-underline-offset: 2px;
      }
      .glossary-panel {
        margin: 0;
      }
      .glossary-panel .shiny-input-container {
        margin-bottom: 2px !important;
      }
      .glossary-definition {
        min-height: 18px;
        font-size: 0.75rem;
        line-height: 1.15;
        color: white;
        margin: 1px 0 0;
        overflow-wrap: anywhere;
      }
      .expression-scale-row {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 6px;
      }
      .expression-panel { width: 100%; overflow: hidden; }
      .expression-title {
        min-height: 42px; padding: 9px 13px; background: #111827;
        color: white; font-weight: 600;
      }
      .plot-wrapper {
        width: 100%; height: calc(95vh - 233px); min-height: 360px;
        padding: 6px; overflow: hidden; box-sizing: border-box;
      }
      .matrix-statistics-note {
        min-height: 30px; padding: 4px 12px 8px; color: #374151;
        font-size: 0.82rem; line-height: 1.25; text-align: center;
      }
      .plot-wrapper .shiny-plot-output,
      .plot-wrapper .shiny-spinner-output-container,
      .plot-wrapper .load-container {
        width: 100% !important; height: 100% !important;
        min-height: 0 !important;
      }
      @media (max-width: 1000px) {
        .mouse-title { font-size: 15px; max-width: 55vw; overflow: hidden; text-overflow: ellipsis; }
        .mouse-header-links span { display: none; }
        .mouse-header-links { width: 80px; }
        .mouse-logo-link { width: 110px; }
        .mouse-logo { max-width: 105px; }
        .mouse-header-spacer {
          width: 80px;
        }
      }
      "
      ))
    ),
    tags$script(HTML(paste(
      c(
        '(function() {',
        '  const dropdownHelp = {',
        '    gene: "Click the drop-down box to select or search for a gene symbol.",',
        '    filter_fields: "Select one or more metadata fields used to filter the loaded gene data.",',
        '    plot_type: "Click the drop-down box to select the type of plot to view.",',
        '    trajectory_measure: "Choose whether the plot displays gene expression or the fraction of cells expressing the gene.",',
        '    progression_variable: "Select the ordered developmental variable used on the trajectory axis.",',
        '    facet_variable: "Select the metadata field used to divide the trajectory into separate panels.",',
        '    color_variable: "Select the metadata field used to color trajectory observations and smoothers.",',
        '    smoother: "Select the smoothing method used to summarize the trajectory.",',
        '    facets_per_row: "Choose the maximum number of facet panels displayed in each row.",',
        '    facet_row_height: "Choose the vertical height assigned to each row of facet panels.",',
        '    violin_display: "Choose whether group distributions are displayed as violins or boxplots.",',
        '    violin_facets_per_page: "Choose how many facets are rendered on each trajectory or violin page.",',
        '    violin_page: "Choose which page of facets to display.",',
        '    comparison_gene: "Select or search for the second gene used in the correlation plot.",',
        '    correlation_color_variable: "Select the metadata field used to color points in the two-gene plot.",',
        '    x_variable: "Select the metadata field shown on the horizontal axis.",',
        '    second_dimension: "Select the metadata field shown on the vertical axis or used for violin facets."',
        '  };',
        '  function applyDropdownHelp() {',
        '    document.querySelectorAll("select[id]").forEach(function(select) {',
        '      const container = select.closest(".form-group, .shiny-input-container");',
        '      const label = container ? container.querySelector("label") : null;',
        '      if (!label) return;',
        '      const cleanLabel = label.textContent.trim();',
        '      const help = dropdownHelp[select.id] ||',
        '        ("Click the drop-down box to select " + cleanLabel.toLowerCase() + ".");',
        '      label.setAttribute("title", help);',
        '      label.classList.add("dropdown-help");',
        '    });',
        '  }',
        '  document.addEventListener("DOMContentLoaded", function() {',
        '    applyDropdownHelp();',
        '    new MutationObserver(applyDropdownHelp).observe(',
        '      document.body,',
        '      { childList: true, subtree: true }',
        '    );',
        '  });',
        '  document.addEventListener("shiny:connected", applyDropdownHelp);',
        '})();'
      ),
      collapse = "\n"
    ))),
    
    sidebar = sidebar(
      width = 320,
      div(
        class = "sidebar-logo-container",
        tags$a(
          href = "https://alleninstitute.github.io/abc_atlas_access/descriptions/Dev-Mouse-Vis-Cortex-dataset.html",
          target = "_blank",
          tags$img(
            src = "Developing_Mouse_Visual_Cortex_logo.png",
            alt = app_title,
            class = "sidebar-logo"
          )
        )
      ),
      helpText(
        "Choose a gene, retrieve its data, select plot options, and generate the plot."
      ),
      div(
        class = "gene-controls",
        selectizeInput(
          "gene",
          "Gene symbol",
          choices = NULL,
          selected = NULL,
          options = list(
            placeholder = paste0(
              "Type a gene symbol (e.g., ",
              default_gene,
              ")"
            ),
            maxOptions = 50,
            create = FALSE
          )
        ),
        actionButton(
          "get_gene_data",
          "Get gene data",
          class = "btn-outline-primary",
          width = "100%"
        ),
        div(class = "gene-status", textOutput("gene_status")),
        conditionalPanel(
          condition = "output.gene_loaded",
          uiOutput("plot_in_abc_atlas_button")
        )
      ),
      conditionalPanel(
        condition = "output.gene_loaded",
        div(
          class = "compact-controls",
          hr(),
          h4("Filter and scale"),
          selectizeInput(
            "filter_fields",
            "Filter metadata",
            choices = filter_fields,
            selected = character(),
            multiple = TRUE,
            options = list(
              placeholder = "Add one or more filters",
              maxOptions = 1000,
              closeAfterSelect = TRUE,
              plugins = list("remove_button")
            )
          ),
          div(
            class = "stacked-filter-panel",
            uiOutput("stacked_filter_controls")
          ),
          checkboxInput(
            "omit_zero_values",
            "Discard observations with zero counts",
            FALSE
          ),
          hr(),
          br(),
          h4("Plot"),
          selectInput(
            "plot_type",
            "Plot type",
            choices = c(
              "Trajectory across development" = "trajectory",
              "Heatmap of mean expression" = "heatmap",
              "Dot plot of mean expression" = "dot",
              "Violin plot with observations" = "violin",
              "Two-gene correlation" = "correlation"
            ),
            selected = "trajectory"
          ),
          conditionalPanel(
            condition = "input.plot_type != 'dot'",
            selectInput(
              "trajectory_measure",
              "Gene metric",
              choices = c(
                "Gene expression" = "expression",
                "Fraction of cells expressing" = "fraction"
              ),
              selected = "expression"
            )
          ),
          conditionalPanel(
            condition = "input.plot_type == 'trajectory'",
            selectInput(
              "progression_variable",
              "Developmental progression axis",
              choices = age_field,
              selected = age_field
            ),
            selectInput(
              "facet_variable",
              "Facet by (maximum 50 values)",
              choices = c(
                "Show all data together" = "none",
                cell_type_fields,
                region_field
              ),
              selected = default_facet
            ),
            selectInput(
              "color_variable",
              "Color by (maximum 50 values)",
              choices = c(
                "All data" = "none",
                region_field
              ),
              selected = default_color
            ),
            selectInput(
              "smoother",
              "Trend line",
              choices = c(
                "LOESS smoother" = "loess",
                "Linear fit" = "lm",
                "None" = "none"
              ),
              selected = "loess"
            )
          ),
          conditionalPanel(
            condition = "input.plot_type != 'trajectory' && input.plot_type != 'correlation'",
            selectInput(
              "x_variable",
              "Horizontal axis",
              choices = plot_fields,
              selected = default_x
            ),
            selectInput(
              "second_dimension",
              "Second dimension",
              choices = setdiff(plot_fields, default_x),
              selected = default_second
            )
          ),
          conditionalPanel(
            condition = "input.plot_type == 'violin' || (input.plot_type == 'trajectory' && input.facet_variable != 'none')",
            div(
              class = "filter-row",
              selectInput(
                "facets_per_row",
                "Facets per row",
                choices = c(
                  "1" = 1,
                  "2" = 2,
                  "3" = 3,
                  "4" = 4
                ),
                selected = 3
              ),
              selectInput(
                "facet_row_height",
                "Height per facet",
                choices = c(
                  "Compact" = 225,
                  "Standard" = 300,
                  "Tall" = 450,
                  "Extra tall" = 600
                ),
                selected = 300
              )
            )
          ),
          conditionalPanel(
            condition = "input.plot_type == 'violin'",
            selectInput(
              "violin_display",
              "Distribution display",
              choices = c(
                "Violin" = "violin",
                "Boxplot" = "boxplot"
              ),
              selected = "violin"
            ),
            checkboxInput(
              "show_violin_points",
              "Show individual observations",
              TRUE
            )
          ),
          conditionalPanel(
            condition = "input.plot_type == 'violin' || (input.plot_type == 'trajectory' && input.facet_variable != 'none')",
            div(
              class = "filter-row",
              selectInput(
                "violin_facets_per_page",
                "Facets per page",
                choices = c(
                  "8" = 8,
                  "12" = 12,
                  "20" = 20,
                  "All" = "all"
                ),
                selected = 12
              ),
              uiOutput("violin_page_control")
            )
          ),
          conditionalPanel(
            condition = "input.plot_type == 'correlation'",
            selectizeInput(
              "comparison_gene",
              "Comparison gene symbol",
              choices = NULL,
              selected = NULL,
              options = list(
                placeholder = "Type a comparison gene symbol",
                maxOptions = 50,
                create = FALSE
              )
            ),
            actionButton(
              "get_comparison_gene_data",
              "Gene 2nd gene data",
              class = "btn-outline-primary",
              width = "100%"
            ),
            div(
              class = "gene-status",
              textOutput("comparison_gene_status")
            ),
            selectInput(
              "correlation_color_variable",
              "Color by (maximum 50 values)",
              choices = c(
                "All data" = "none",
                region_field
              ),
              selected = default_color
            ),
            checkboxInput(
              "correlation_automatic_limits",
              "Use automatic axis limits",
              TRUE
            ),
            conditionalPanel(
              condition = "!input.correlation_automatic_limits",
              div(
                class = "expression-scale-row",
                numericInput(
                  "correlation_axis_minimum",
                  "Minimum",
                  value = 0,
                  min = 0
                ),
                numericInput(
                  "correlation_axis_maximum",
                  "Maximum",
                  value = 1,
                  min = 0
                )
              )
            ),
            checkboxInput(
              "show_orthogonal_fit",
              "Show orthogonal fit line",
              TRUE
            ),
            div(
              class = "expression-scale-row",
              numericInput(
                "correlation_size_minimum",
                "Min. dot size",
                value = 1.5,
                min = 0,
                step = 0.5
              ),
              numericInput(
                "correlation_size_maximum",
                "Max. dot size",
                value = 7,
                min = 0.5,
                step = 0.5
              )
            )
          ),
          conditionalPanel(
            condition = "input.plot_type != 'correlation'",
            conditionalPanel(
              condition = "input.trajectory_measure == 'expression'",
              checkboxInput(
                "log_scale",
                "Plot ln(CPM + 1)",
                TRUE
              )
            ),
            checkboxInput(
              "automatic_expression_limits",
              "Use automatic expression limits",
              TRUE
            ),
            conditionalPanel(
              condition = "!input.automatic_expression_limits",
              div(
                class = "expression-scale-row",
                numericInput(
                  "expression_minimum",
                  "Minimum",
                  value = 0,
                  min = 0
                ),
                numericInput(
                  "expression_maximum",
                  "Maximum",
                  value = 1,
                  min = 0
                )
              )
            ),
            conditionalPanel(
              condition = "input.plot_type == 'dot'",
              div(
                class = "expression-scale-row",
                numericInput(
                  "dot_size_minimum",
                  "Min. dot size",
                  value = 1,
                  min = 0,
                  step = 0.5
                ),
                numericInput(
                  "dot_size_maximum",
                  "Max. dot size",
                  value = 9,
                  min = 0.5,
                  step = 0.5
                )
              )
            ),
            conditionalPanel(
              condition = "input.plot_type == 'heatmap'",
              checkboxInput(
                "show_heatmap_counts",
                "Show number of observations",
                TRUE
              )
            ),
            conditionalPanel(
              condition = "input.plot_type == 'trajectory'",
              checkboxInput(
                "show_points",
                "Show individual observations",
                TRUE
              )
            )
          ),
          actionButton(
            "make_plot",
            "Generate plot",
            class = "btn-primary",
            width = "100%"
          ),
          actionButton(
            "share_view",
            "Share this view",
            icon = icon("share-nodes"),
            class = "btn-outline-light btn-sm mt-1",
            width = "100%"
          ),
          actionButton(
            "reset_defaults",
            "Reset plot options",
            class = "btn-outline-light btn-sm mt-1",
            width = "100%"
          )
        )
      ),
      if (tokens_available) {
        tagList(
          hr(style = "margin: 3px 0;"),
          div(
            class = "glossary-panel",
            selectizeInput(
              "glossary_token",
              "Glossary of terms",
              choices = setNames(tokens$token, tokens$token),
              selected = character(),
              multiple = FALSE,
              options = list(
                placeholder = "Type or select a term",
                maxOptions = nrow(tokens),
                create = FALSE
              )
            ),
            div(
              class = "glossary-definition",
              textOutput("glossary_definition", inline = TRUE)
            )
          ),
          hr(style = "margin: 3px 0;")
        )
      },
      
      h4(
        "Contribute",
        style = "font-size: 1.2rem; margin: 1px 0 0px; color: white;"
      ),
      div(
        style = paste(
          "display: flex; gap: 6px; width: 100%;",
          "align-items: stretch; margin-bottom: 1px;"
        ),
        tags$a(
          href = "mailto:jeremym@alleninstitute.org?body=''&subject='Dev mouse' app comments",
          icon("envelope", lib = "font-awesome"),
          tags$span("PROVIDE FEEDBACK"),
          style = paste(
            "display: inline-flex; flex: 1 1 0;",
            "align-items: center; justify-content: center; gap: 4px;",
            "min-height: 29px; padding: 4px 6px;",
            "background: #666666; border: 1px solid white;",
            "border-radius: 0.25rem; color: white !important;",
            "font-size: 0.68rem; font-weight: 600;",
            "line-height: 1; text-decoration: none; white-space: nowrap;"
          )
        ),
        tags$a(
          href = "https://github.com/AllenInstitute/mouse-developing-vis-app/",
          target = "_blank",
          rel = "noopener noreferrer",
          icon("code", lib = "font-awesome"),
          tags$span("ACCESS CODE"),
          style = paste(
            "display: inline-flex; flex: 1 1 0;",
            "align-items: center; justify-content: center; gap: 4px;",
            "min-height: 29px; padding: 4px 6px;",
            "background: #666666; border: 1px solid white;",
            "border-radius: 0.25rem; color: white !important;",
            "font-size: 0.68rem; font-weight: 600;",
            "line-height: 1; text-decoration: none; white-space: nowrap;"
          )
        )
      ),
      p(
        "We encourage contributions via the above links.",
        style = "font-size: 0.9rem; line-height: 1.15; margin: 1px 0 5px; color: white;"
      ),
      h4(
        "Acknowledgements",
        style = "font-size: 1.2rem; margin: 1px 0 0px; color: white;"
      ),
      p(
        "App developed by Jeremy Miller, coded with help from Microsoft Copilot (model GPT 5.6 Think). All data included in the app are from ",
        a(
          "Continuous cell-type diversification in mouse visual cortex development",
          href = "https://doi.org/10.1038/s41586-025-09644-1",
          target = "_blank",
          rel = "noopener noreferrer",
          style = "color: white; font-weight: bold;"
        ),
        " and are freely available on ",
        a(
          "ABC Atlas Access",
          href = "https://alleninstitute.github.io/abc_atlas_access/descriptions/Dev-Mouse-Vis-Cortex-dataset.html",
          target = "_blank",
          rel = "noopener noreferrer",
          style = "color: white; font-weight: bold;"
        ),
        ". Beta testers include: Irika Sinha.",
        style = "font-size: 0.9rem; line-height: 1.15; margin: 1px 0 5px; color: white;"
      )
    ),
    
    navset_card_tab(
      id = "main_tabs",
      nav_panel(
        "Select gene of interest",
        div(
          class = "gene-table-panel",
          DT::DTOutput("gene_statistics_table"),
          tags$hr(),
          p(
            class = "statistics-value-filter-instructions",
            paste(
              "Use this section to filter the table by specific category",
              "values rather than using a text-matching search."
            )
          ),
          div(
            class = "statistics-value-filter-panel",
            selectizeInput(
              "statistics_filter_gene_symbol",
              "gene_symbol",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = "Select exact gene symbols",
                maxOptions = 100,
                closeAfterSelect = FALSE,
                plugins = list("remove_button")
              )
            ),
            selectizeInput(
              "statistics_filter_max_ROI",
              "max_ROI",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = "Select max_ROI values",
                closeAfterSelect = FALSE,
                plugins = list("remove_button")
              )
            ),
            selectizeInput(
              "statistics_filter_max_subclass",
              "max_subclass",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = "Select max_subclass values",
                closeAfterSelect = FALSE,
                plugins = list("remove_button")
              )
            ),
            selectizeInput(
              "statistics_filter_gene_type",
              "gene_type",
              choices = NULL,
              selected = NULL,
              multiple = TRUE,
              options = list(
                placeholder = "Select gene_type values",
                closeAfterSelect = FALSE,
                plugins = list("remove_button")
              )
            )
          )
        )
      ),
      nav_panel(
        "Plot gene of interest",
        div(
          class = "expression-panel",
          div(
            class = "expression-title",
            textOutput("plot_title")
          ),
          div(
            class = "plot-wrapper",
            shinycssloaders::withSpinner(
              plotOutput(
                "expression_plot",
                width = "100%",
                height = "100%"
              ),
              type = 8,
              color = "#4B9B58",
              size = 1
            )
          ),
          div(class = "matrix-statistics-note", textOutput("plot_statistics_note"))
        )
      )
    )
  )
}

# ============================================================
# Server
# ============================================================

server <- function(input, output, session) {
  output$glossary_definition <- renderText({
    req(tokens_available)
    selected_token <- input$glossary_token
    if (
      is.null(selected_token) ||
      length(selected_token) != 1 ||
      is.na(selected_token) ||
      !nzchar(selected_token)
    ) {
      return("")
    }
    match_index <- match(selected_token, tokens$token)
    if (is.na(match_index)) return("")
    paste0(
      tokens$token[[match_index]],
      ": ",
      tokens$name[[match_index]],
      " (",
      tokens$type[[match_index]],
      ")"
    )
  })
  gene_data <- reactiveVal(NULL)
  loaded_gene <- reactiveVal(NULL)
  requested_gene <- reactiveVal(default_gene)
  gene_status_message <- reactiveVal("No gene retrieved yet.")
  abc_atlas_url <- reactiveVal(NULL)
  comparison_gene_data <- reactiveVal(NULL)
  loaded_comparison_gene <- reactiveVal(NULL)
  comparison_gene_status_message <- reactiveVal(
    "No comparison gene retrieved yet."
  )
  previous_dimension_plot_type <- reactiveVal(NULL)
  bookmark_plot_trigger <- reactiveVal(0L)
  pending_bookmark_state <- reactiveVal(NULL)
  restored_filter_values <- reactiveVal(NULL)
  restore_plot_pending <- reactiveVal(FALSE)
  bookmark_schema_version <- 1L
  
  bookmark_excluded_inputs <- unique(c(
    "gene",
    "comparison_gene",
    "filter_fields",
    "omit_zero_values",
    "plot_type",
    "trajectory_measure",
    "progression_variable",
    "facet_variable",
    "color_variable",
    "smoother",
    "x_variable",
    "second_dimension",
    "correlation_color_variable",
    "correlation_automatic_limits",
    "correlation_axis_minimum",
    "correlation_axis_maximum",
    "show_orthogonal_fit",
    "correlation_size_minimum",
    "correlation_size_maximum",
    "log_scale",
    "automatic_expression_limits",
    "expression_minimum",
    "expression_maximum",
    "dot_size_minimum",
    "dot_size_maximum",
    "show_heatmap_counts",
    "show_points",
    "violin_display",
    "violin_facets_per_page",
    "violin_page",
    "show_violin_points",
    "facets_per_row",
    "facet_row_height",
    "statistics_filter_gene_symbol",
    "statistics_filter_max_ROI",
    "statistics_filter_max_subclass",
    "statistics_filter_gene_type",
    "get_gene_data",
    "get_comparison_gene_data",
    "make_plot",
    "reset_defaults",
    "share_view",
    "bookmarked_url",
    "gene_statistics_table_rows_selected",
    "gene_statistics_table_search",
    "gene_statistics_table_rows_current",
    "gene_statistics_table_cell_clicked",
    vapply(filter_fields, filter_input_id, character(1))
  ))
  setBookmarkExclude(bookmark_excluded_inputs)
  
  observe({
    setBookmarkExclude(names(input))
  }, priority = 1000)
  
  onBookmark(function(state) {
    selected_filter_fields <- isolate(input$filter_fields)
    if (is.null(selected_filter_fields)) {
      selected_filter_fields <- character()
    }
    selected_filter_fields <- selected_filter_fields[
      selected_filter_fields %in% filter_fields
    ]
    
    stacked_filter_values <- setNames(
      lapply(selected_filter_fields, function(field) {
        isolate(input[[filter_input_id(field)]])
      }),
      selected_filter_fields
    )
    
    state$values$viewer_state <- list(
      schema_version = bookmark_schema_version,
      gene = isolate(loaded_gene()),
      comparison_gene = isolate(loaded_comparison_gene()),
      filter_fields = selected_filter_fields,
      stacked_filter_values = stacked_filter_values,
      omit_zero_values = isolate(input$omit_zero_values),
      plot_type = isolate(input$plot_type),
      trajectory_measure = isolate(input$trajectory_measure),
      progression_variable = isolate(input$progression_variable),
      facet_variable = isolate(input$facet_variable),
      color_variable = isolate(input$color_variable),
      smoother = isolate(input$smoother),
      x_variable = isolate(input$x_variable),
      second_dimension = isolate(input$second_dimension),
      correlation_color_variable = isolate(
        input$correlation_color_variable
      ),
      correlation_automatic_limits = isolate(
        input$correlation_automatic_limits
      ),
      correlation_axis_minimum = isolate(
        input$correlation_axis_minimum
      ),
      correlation_axis_maximum = isolate(
        input$correlation_axis_maximum
      ),
      show_orthogonal_fit = isolate(input$show_orthogonal_fit),
      correlation_size_minimum = isolate(input$correlation_size_minimum),
      correlation_size_maximum = isolate(input$correlation_size_maximum),
      log_scale = isolate(input$log_scale),
      automatic_expression_limits = isolate(
        input$automatic_expression_limits
      ),
      expression_minimum = isolate(input$expression_minimum),
      expression_maximum = isolate(input$expression_maximum),
      dot_size_minimum = isolate(input$dot_size_minimum),
      dot_size_maximum = isolate(input$dot_size_maximum),
      show_heatmap_counts = isolate(input$show_heatmap_counts),
      show_points = isolate(input$show_points),
      violin_display = isolate(input$violin_display),
      violin_facets_per_page = isolate(input$violin_facets_per_page),
      violin_page = isolate(input$violin_page),
      show_violin_points = isolate(input$show_violin_points),
      facets_per_row = isolate(input$facets_per_row),
      facet_row_height = isolate(input$facet_row_height)
    )
  })
  
  onRestore(function(state) {
    saved <- state$values$viewer_state
    pending_bookmark_state(saved)
    if (!is.null(saved)) {
      saved_fields <- as.character(saved$filter_fields)
      saved_fields <- saved_fields[saved_fields %in% filter_fields]
      restored_filter_values(list(
        fields = saved_fields,
        values = saved$stacked_filter_values
      ))
      restore_plot_pending(TRUE)
    }
  })
  
  onBookmarked(function(url) {
    showModal(modalDialog(
      title = "Share this view",
      p(
        "Copy this URL to share the current gene, filters, and plot settings."
      ),
      textInput(
        "bookmarked_url",
        NULL,
        value = url,
        width = "100%"
      ),
      tags$a(
        href = url,
        target = "_blank",
        rel = "noopener noreferrer",
        "Open bookmarked view"
      ),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  observeEvent(input$share_view, {
    if (is.null(loaded_gene())) {
      showNotification(
        "Retrieve a gene before sharing this view.",
        type = "warning"
      )
      return()
    }
    session$doBookmark()
  }, ignoreInit = TRUE)
  
  output$gene_loaded <- reactive({
    !is.null(gene_data()) &&
      !is.null(loaded_gene())
  })
  outputOptions(
    output,
    "gene_loaded",
    suspendWhenHidden = FALSE
  )
  
  statistics_value_filter_columns <- c(
    gene_symbol = statistics_gene_column,
    max_ROI = "max_ROI",
    max_subclass = "max_subclass",
    gene_type = "gene_type"
  )
  statistics_value_filter_columns <- statistics_value_filter_columns[
    statistics_value_filter_columns %in% names(gene_statistics)
  ]
  
  statistics_filter_choices <- function(column_name) {
    values <- unique(as.character(gene_statistics[[column_name]]))
    values <- values[!is.na(values) & nzchar(values)]
    sort(values, na.last = TRUE)
  }
  
  session$onFlushed(function() {
    for (display_name in names(statistics_value_filter_columns)) {
      column_name <- unname(
        statistics_value_filter_columns[[display_name]]
      )
      updateSelectizeInput(
        session,
        paste0("statistics_filter_", display_name),
        choices = statistics_filter_choices(column_name),
        selected = character(),
        server = TRUE
      )
    }
  }, once = TRUE)
  
  filtered_gene_statistics <- reactive({
    filtered_table <- gene_statistics
    
    for (display_name in names(statistics_value_filter_columns)) {
      column_name <- unname(
        statistics_value_filter_columns[[display_name]]
      )
      selected_values <- input[[
        paste0("statistics_filter_", display_name)
      ]]
      
      if (!is.null(selected_values) && length(selected_values) > 0) {
        filtered_table <- filtered_table[
          as.character(filtered_table[[column_name]]) %in%
            as.character(selected_values),
          ,
          drop = FALSE
        ]
      }
    }
    
    filtered_table
  })
  
  output$gene_statistics_table <- DT::renderDT({
    table_data <- filtered_gene_statistics()
    
    numeric_columns <- names(table_data)[vapply(
      table_data,
      is.numeric,
      logical(1)
    )]
    
    first_numeric_index <- if (length(numeric_columns) > 0) {
      match(numeric_columns[[1]], names(table_data)) - 1L
    } else {
      0L
    }
    
    gene_column_index <- match(
      statistics_gene_column,
      names(table_data)
    ) - 1L
    
    header_definitions <- setNames(
      table_column_definitions$column_definitions,
      table_column_definitions$column_names
    )
    header_definitions <- unname(
      header_definitions[names(table_data)]
    )
    header_definitions[is.na(header_definitions)] <- ""
    header_definitions_json <- jsonlite::toJSON(
      header_definitions,
      auto_unbox = TRUE
    )
    
    table <- DT::datatable(
      table_data,
      rownames = FALSE,
      filter = "top",
      selection = list(
        mode = "single",
        target = "row"
      ),
      options = list(
        pageLength = 10,
        lengthMenu = c(10, 20, 50, 100),
        scrollX = TRUE,
        autoWidth = TRUE,
        stateSave = FALSE,
        searchHighlight = TRUE,
        order = list(
          list(first_numeric_index, "desc")
        ),
        initComplete = DT::JS(
          "function(settings, json) {",
          paste0("  var definitions = ", header_definitions_json, ";"),
          "  var api = this.api();",
          "  api.columns().every(function(index) {",
          "    var definition = definitions[index] || '';",
          "    if (definition !== '') {",
          "      var header = $(this.header());",
          "      header.attr('title', definition);",
          "      header.attr('data-bs-toggle', 'tooltip');",
          "      header.attr('data-bs-placement', 'bottom');",
          "      header.addClass('has-definition');",
          "      if (window.bootstrap && bootstrap.Tooltip) {",
          "        bootstrap.Tooltip.getOrCreateInstance(header[0], {container: 'body'});",
          "      }",
          "    }",
          "  });",
          "}"
        ),
        columnDefs = list(
          list(
            targets = gene_column_index,
            render = DT::JS(
              "function(data, type, row, meta) {",
              "  if (type !== 'display' || data === null || data === '') return data;",
              "  var gene = String(data);",
              "  var href = 'https://www.genecards.org/card/' + encodeURIComponent(gene);",
              "  var label = $('<div>').text(gene).html();",
              "  return '<a href=\"' + href + '\" target=\"_blank\" rel=\"noopener noreferrer\">' + label + '</a>';",
              "}"
            )
          )
        )
      ),
      class = "compact stripe hover"
    )
    
    if (length(numeric_columns) > 0) {
      table <- DT::formatRound(
        table,
        columns = numeric_columns,
        digits = 3
      )
    }
    
    table
  }, server = TRUE)
  
  observeEvent(input$gene_statistics_table_rows_selected, {
    selected_row <- input$gene_statistics_table_rows_selected
    req(length(selected_row) == 1)
    
    current_statistics <- filtered_gene_statistics()
    selected_gene <- current_statistics[[statistics_gene_column]][selected_row]
    req(length(selected_gene) == 1, !is.na(selected_gene), nzchar(selected_gene))
    
    if (selected_gene %in% available_genes) {
      requested_gene(selected_gene)
      updateSelectizeInput(
        session,
        "gene",
        choices = available_genes,
        selected = selected_gene,
        server = TRUE
      )
      gene_status_message(
        paste("Selected", selected_gene, "from the statistics table.")
      )
    } else {
      showNotification(
        paste("This gene is not available:", selected_gene),
        type = "warning"
      )
    }
  }, ignoreInit = TRUE)
  
  session$onFlushed(function() {
    updateSelectizeInput(
      session,
      "gene",
      choices = available_genes,
      selected = character(),
      server = TRUE
    )
    updateSelectizeInput(
      session,
      "comparison_gene",
      choices = available_genes,
      selected = character(),
      server = TRUE
    )
  }, once = TRUE)
  
  observeEvent(input$gene, {
    if (
      length(input$gene) == 1 &&
      nzchar(input$gene) &&
      input$gene %in% available_genes
    ) {
      requested_gene(input$gene)
    }
  }, ignoreInit = TRUE)
  
  read_gene <- function(gene_symbol) {
    cache_key <- make_cache_key(gene_symbol)
    cached <- gene_cache$get(cache_key, missing = NULL)
    if (!is.null(cached)) return(cached)
    
    gene_key <- unname(gene_symbol_to_key[[gene_symbol]])
    if (is.null(gene_key) || is.na(gene_key)) {
      stop("No gene partition key was found for: ", gene_symbol)
    }
    
    gene_source <- bucket$path(
      file.path(
        s3_prefix,
        "counts_by_gene",
        paste0("gene=", gene_key)
      )
    )
    
    # open_dataset() deliberately supports one or many Parquet files in the
    # selected gene partition. Detection counts are optional during migration.
    gene_dataset <- arrow::open_dataset(
      gene_source,
      format = "parquet",
      partitioning = NULL,
      unify_schemas = FALSE
    )
    has_detection_counts <- "number_of_cells_expressing" %in%
      gene_dataset$schema$names
    
    if (has_detection_counts) {
      gene_counts <- gene_dataset |>
        select(sample_id, summed_counts, number_of_cells_expressing) |>
        collect() |>
        transmute(
          sample_id = as.character(sample_id),
          summed_counts = as.numeric(summed_counts),
          number_of_cells_expressing = as.numeric(number_of_cells_expressing)
        )
    } else {
      gene_counts <- gene_dataset |>
        select(sample_id, summed_counts) |>
        collect() |>
        transmute(
          sample_id = as.character(sample_id),
          summed_counts = as.numeric(summed_counts),
          number_of_cells_expressing = NA_real_
        )
    }
    
    if (nrow(gene_counts) == 0) {
      stop("No count data were found for gene: ", gene_symbol)
    }
    
    joined <- metadata |>
      inner_join(gene_counts, by = "sample_id") |>
      filter(
        is.finite(summed_counts),
        is.finite(CPM_scaling_factor),
        CPM_scaling_factor >= 0
      )
    
    if (!"number_of_cells" %in% names(joined)) {
      joined$number_of_cells <- NA_real_
    }
    
    joined <- joined |>
      mutate(
        number_of_cells = suppressWarnings(as.numeric(number_of_cells)),
        number_of_cells_expressing = suppressWarnings(
          as.numeric(number_of_cells_expressing)
        ),
        valid_detection_counts =
          is.finite(number_of_cells) & number_of_cells > 0 &
          is.finite(number_of_cells_expressing) &
          number_of_cells_expressing >= 0 &
          number_of_cells_expressing <= number_of_cells,
        fraction_expressing = if_else(
          valid_detection_counts,
          number_of_cells_expressing / number_of_cells,
          NA_real_
        ),
        CPM = summed_counts * CPM_scaling_factor
      )
    
    if (nrow(joined) == 0) {
      stop("Gene data were found, but no rows matched valid metadata and CPM scaling.")
    }
    
    gene_cache$set(cache_key, joined)
    joined
  }
  
  observeEvent(input$get_gene_data, {
    gene_to_load <- requested_gene()
    req(gene_to_load %in% available_genes)
    gene_status_message(paste("Retrieving", gene_to_load, "..."))
    
    tryCatch(
      {
        loaded <- withProgress(
          message = paste("Retrieving", gene_to_load),
          value = 0.25,
          {
            result <- read_gene(gene_to_load)
            incProgress(0.65)
            result
          }
        )
        gene_data(loaded)
        loaded_gene(gene_to_load)
        updateSelectizeInput(
          session,
          "gene",
          choices = available_genes,
          selected = gene_to_load,
          server = TRUE
        )
        gene_status_message(
          paste0(
            "Loaded ", gene_to_load, " (",
            comma(nrow(loaded)), " observations)."
          )
        )
      },
      error = function(e) {
        gene_status_message(
          paste("Unable to load", gene_to_load, ":", conditionMessage(e))
        )
        showNotification(
          paste("Unable to load gene data:", conditionMessage(e)),
          type = "error",
          duration = NULL
        )
      }
    )
  })
  
  output$gene_status <- renderText(gene_status_message())
  
  observeEvent(loaded_gene(), {
    gene_symbol <- loaded_gene()
    abc_atlas_url(NULL)
    
    req(
      !is.null(gene_data()),
      length(gene_symbol) == 1,
      !is.na(gene_symbol),
      nzchar(gene_symbol)
    )
    
    matched_url <- unname(
      abc_atlas_url_lookup[gene_symbol]
    )
    
    if (
      length(matched_url) == 1L &&
      !is.na(matched_url) &&
      nzchar(matched_url)
    ) {
      abc_atlas_url(matched_url)
    }
  }, ignoreInit = TRUE)
  
  output$plot_in_abc_atlas_button <- renderUI({
    url <- abc_atlas_url()
    req(
      !is.null(gene_data()),
      !is.null(loaded_gene()),
      length(url) == 1,
      !is.na(url),
      nzchar(url)
    )
    
    tags$a(
      href = url,
      target = "_blank",
      rel = "noopener noreferrer",
      class = "btn btn-outline-primary w-100",
      icon("arrow-up-right-from-square"),
      "Plot in ABC Atlas"
    )
  })
  
  observeEvent(input$comparison_gene, {
    selected_gene <- input$comparison_gene
    if (identical(input$plot_type, "correlation")) {
      validate(
        need(
          !is.null(input$comparison_gene) &&
            length(input$comparison_gene) == 1 &&
            nzchar(input$comparison_gene),
          "Select a comparison gene."
        ),
        need(
          !is.null(comparison_gene_data()) &&
            !is.null(loaded_comparison_gene()),
          paste(
            "Click 'Get comparison gene data'",
            "before generating the correlation plot."
          )
        ),
        need(
          identical(
            input$comparison_gene,
            loaded_comparison_gene()
          ),
          paste(
            "The selected comparison gene has not been retrieved.",
            "Click 'Get comparison gene data'."
          )
        ),
        need(
          !identical(
            loaded_gene(),
            loaded_comparison_gene()
          ),
          "Choose two different genes for the correlation plot."
        )
      )
    }
  }, ignoreInit = TRUE)
  
  observeEvent(input$get_comparison_gene_data, {
    comparison_gene <- input$comparison_gene
    req(
      length(comparison_gene) == 1,
      nzchar(comparison_gene),
      comparison_gene %in% available_genes
    )
    validate(need(
      !identical(comparison_gene, loaded_gene()),
      "Choose a comparison gene different from the primary gene."
    ))
    comparison_gene_status_message(
      paste("Retrieving", comparison_gene, "...")
    )
    
    tryCatch(
      {
        loaded <- withProgress(
          message = paste("Retrieving", comparison_gene),
          value = 0.25,
          {
            result <- read_gene(comparison_gene)
            incProgress(0.65)
            result
          }
        )
        comparison_gene_data(loaded)
        loaded_comparison_gene(comparison_gene)
        comparison_gene_status_message(
          paste0(
            "Loaded ", comparison_gene, " (",
            comma(nrow(loaded)), " observations)."
          )
        )
      },
      error = function(e) {
        comparison_gene_data(NULL)
        loaded_comparison_gene(NULL)
        comparison_gene_status_message(
          paste(
            "Unable to load ", comparison_gene, ": ",
            conditionMessage(e)
          )
        )
        showNotification(
          paste(
            "Unable to load comparison gene data:",
            conditionMessage(e)
          ),
          type = "error",
          duration = NULL
        )
      }
    )
  })
  
  output$comparison_gene_status <- renderText(
    comparison_gene_status_message()
  )
  
  automatic_expression_range <- reactive({
    data <- gene_data()
    req(data)
    
    filtered <- apply_filter_specification(
      data,
      active_filter_specification()
    )
    if (isTRUE(input$omit_zero_values)) {
      filtered <- filtered |>
        filter(summed_counts > 0)
    }
    
    values <- if (
      input$plot_type %in% c("trajectory", "heatmap", "violin") &&
      identical(input$trajectory_measure, "fraction")
    ) {
      filtered$fraction_expressing
    } else if (isTRUE(input$log_scale)) {
      log1p(filtered$CPM)
    } else {
      filtered$CPM
    }
    values <- values[is.finite(values)]
    validate(need(
      length(values) > 0,
      "No finite expression values remain after filtering."
    ))
    
    c(
      minimum = 0,
      maximum = if (
        identical(input$plot_type, "trajectory") &&
        identical(input$trajectory_measure, "fraction")
      ) {
        1
      } else {
        max(values, na.rm = TRUE)
      }
    )
  })
  
  observeEvent(automatic_expression_range(), {
    req(isTRUE(input$automatic_expression_limits))
    limits <- automatic_expression_range()
    maximum <- max(limits[["maximum"]], .Machine$double.eps)
    updateNumericInput(
      session,
      "expression_minimum",
      value = limits[["minimum"]],
      min = 0,
      max = maximum
    )
    updateNumericInput(
      session,
      "expression_maximum",
      value = maximum,
      min = 0,
      max = maximum
    )
  }, ignoreInit = TRUE)
  
  output$stacked_filter_controls <- renderUI({
    selected_fields <- input$filter_fields
    if (is.null(selected_fields) || length(selected_fields) == 0) {
      return(
        div(
          class = "small text-white-50",
          "No metadata filters selected."
        )
      )
    }
    
    selected_fields <- selected_fields[
      selected_fields %in% filter_fields
    ]
    
    tagList(lapply(selected_fields, function(field) {
      input_id <- filter_input_id(field)
      restored_filters <- isolate(restored_filter_values())
      restored_value <- if (
        !is.null(restored_filters) &&
        field %in% restored_filters$fields
      ) {
        restored_filters$values[[field]]
      } else {
        NULL
      }
      current_value <- if (!is.null(restored_value)) {
        restored_value
      } else {
        isolate(input[[input_id]])
      }
      
      if (field %in% numeric_filter_fields) {
        numeric_values <- suppressWarnings(
          as.numeric(metadata[[field]])
        )
        numeric_values <- numeric_values[is.finite(numeric_values)]
        
        if (length(numeric_values) == 0) {
          return(NULL)
        }
        
        minimum <- min(numeric_values)
        maximum <- max(numeric_values)
        selected_range <- if (
          length(current_value) == 2 &&
          all(is.finite(as.numeric(current_value)))
        ) {
          pmax(
            minimum,
            pmin(maximum, as.numeric(current_value))
          )
        } else {
          c(minimum, maximum)
        }
        
        step_size <- if (
          all(abs(numeric_values - round(numeric_values)) < 1e-9)
        ) {
          1
        } else {
          max((maximum - minimum) / 100, .Machine$double.eps)
        }
        
        div(
          class = "stacked-filter-control",
          sliderInput(
            input_id,
            field,
            min = minimum,
            max = maximum,
            value = selected_range,
            step = step_size,
            sep = ""
          )
        )
      } else {
        choices <- field_levels(metadata, field)
        selected_values <- as.character(current_value)
        selected_values <- selected_values[
          selected_values %in% choices
        ]
        
        div(
          class = "stacked-filter-control",
          selectizeInput(
            input_id,
            field,
            choices = choices,
            selected = selected_values,
            multiple = TRUE,
            options = list(
              placeholder = "Include values",
              maxOptions = 1000,
              closeAfterSelect = TRUE,
              plugins = list("remove_button")
            )
          )
        )
      }
    }))
  })
  
  active_filter_specification <- reactive({
    selected_fields <- input$filter_fields
    if (is.null(selected_fields) || length(selected_fields) == 0) {
      return(list())
    }
    
    selected_fields <- selected_fields[
      selected_fields %in% filter_fields
    ]
    
    restored_filters <- restored_filter_values()
    filters <- lapply(selected_fields, function(field) {
      restored_value <- if (
        !is.null(restored_filters) &&
        field %in% restored_filters$fields
      ) {
        restored_filters$values[[field]]
      } else {
        NULL
      }
      value <- if (!is.null(restored_value)) {
        restored_value
      } else {
        input[[filter_input_id(field)]]
      }
      
      if (field %in% numeric_filter_fields) {
        if (
          is.null(value) ||
          length(value) != 2 ||
          !all(is.finite(as.numeric(value)))
        ) {
          return(NULL)
        }
        list(
          field = field,
          type = "numeric",
          value = as.numeric(value)
        )
      } else {
        value <- as.character(value)
        value <- value[!is.na(value) & nzchar(value)]
        if (length(value) == 0) {
          return(NULL)
        }
        list(
          field = field,
          type = "categorical",
          value = value
        )
      }
    })
    
    Filter(Negate(is.null), filters)
  })
  
  observe({
    req(isTRUE(restore_plot_pending()))
    saved <- pending_bookmark_state()
    req(!is.null(saved), !is.null(loaded_gene()))
    
    restored_filters <- restored_filter_values()
    filters_ready <- TRUE
    if (!is.null(restored_filters) && length(restored_filters$fields) > 0) {
      filters_ready <- all(vapply(
        restored_filters$fields,
        function(field) {
          expected <- restored_filters$values[[field]]
          actual <- input[[filter_input_id(field)]]
          if (is.null(actual)) return(FALSE)
          if (field %in% numeric_filter_fields) {
            isTRUE(all.equal(
              as.numeric(actual),
              as.numeric(expected),
              tolerance = 1e-8,
              check.attributes = FALSE
            ))
          } else {
            identical(
              sort(as.character(actual)),
              sort(as.character(expected))
            )
          }
        },
        logical(1)
      ))
    }
    req(filters_ready)
    
    restored_filter_values(NULL)
    restore_plot_pending(FALSE)
    bookmark_plot_trigger(isolate(bookmark_plot_trigger()) + 1L)
    bslib::nav_select(
      id = "main_tabs",
      selected = "Plot gene of interest",
      session = session
    )
    pending_bookmark_state(NULL)
  }, priority = -100)
  
  control_metadata <- reactive({
    filtered <- apply_filter_specification(
      metadata,
      active_filter_specification()
    )
    
    if (nrow(filtered) == 0) {
      metadata[0, , drop = FALSE]
    } else {
      filtered
    }
  })
  
  observeEvent(control_metadata(), {
    data <- control_metadata()
    facet_choices <- fields_within_limit(
      data,
      unique(c(cell_type_fields, region_field)),
      50
    )
    color_choices <- fields_within_limit(
      data,
      plot_fields,
      50
    )
    
    facet_choices <- c(
      "Show all data together" = "none",
      facet_choices
    )
    color_choices <- c(
      "All data" = "none",
      color_choices
    )
    
    current_facet <- isolate(input$facet_variable)
    if (!current_facet %in% unname(facet_choices)) {
      current_facet <- if (default_facet %in% unname(facet_choices)) {
        default_facet
      } else {
        "none"
      }
    }
    
    current_color <- isolate(input$color_variable)
    if (!current_color %in% unname(color_choices)) {
      current_color <- if (default_color %in% unname(color_choices)) {
        default_color
      } else {
        "none"
      }
    }
    
    updateSelectInput(
      session,
      "facet_variable",
      choices = facet_choices,
      selected = current_facet
    )
    updateSelectInput(
      session,
      "color_variable",
      choices = color_choices,
      selected = current_color
    )
    current_correlation_color <- isolate(input$correlation_color_variable)
    if (!current_correlation_color %in% unname(color_choices)) {
      current_correlation_color <- if (
        default_color %in% unname(color_choices)
      ) {
        default_color
      } else {
        "none"
      }
    }
    updateSelectInput(
      session,
      "correlation_color_variable",
      choices = color_choices,
      selected = current_correlation_color
    )
  }, ignoreInit = FALSE)
  
  observeEvent(
    list(input$plot_type, input$x_variable, control_metadata()),
    {
      req(input$plot_type)
      data <- control_metadata()
      
      dimension_choices <- if (input$plot_type %in% c("heatmap", "dot")) {
        plot_fields
      } else {
        fields_within_limit(data, plot_fields, 50)
      }
      if (length(dimension_choices) < 2) dimension_choices <- plot_fields
      
      plot_type_changed <- !identical(
        isolate(previous_dimension_plot_type()),
        input$plot_type
      )
      
      if (
        plot_type_changed &&
        input$plot_type %in% c("heatmap", "dot") &&
        default_x %in% dimension_choices
      ) {
        current_x <- default_x
      } else {
        current_x <- isolate(input$x_variable)
        if (
          is.null(current_x) ||
          length(current_x) == 0 ||
          !current_x %in% dimension_choices
        ) {
          current_x <- if (default_x %in% dimension_choices) {
            default_x
          } else {
            dimension_choices[[1]]
          }
        }
      }
      
      updateSelectInput(
        session,
        "x_variable",
        choices = dimension_choices,
        selected = current_x
      )
      
      second_choices <- setdiff(dimension_choices, current_x)
      
      if (
        plot_type_changed &&
        input$plot_type %in% c("heatmap", "dot") &&
        default_second %in% second_choices
      ) {
        current_second <- default_second
      } else {
        current_second <- isolate(input$second_dimension)
        if (
          is.null(current_second) ||
          length(current_second) == 0 ||
          !current_second %in% second_choices
        ) {
          current_second <- if (default_second %in% second_choices) {
            default_second
          } else if (region_field %in% second_choices) {
            region_field
          } else {
            second_choices[[1]]
          }
        }
      }
      
      updateSelectInput(
        session,
        "second_dimension",
        choices = second_choices,
        selected = current_second
      )
      
      previous_dimension_plot_type(input$plot_type)
    },
    ignoreInit = FALSE
  )
  
  
  onRestored(function(state) {
    saved <- pending_bookmark_state()
    if (is.null(saved)) {
      return()
    }
    
    if (
      is.null(saved$schema_version) ||
      !identical(as.integer(saved$schema_version), bookmark_schema_version)
    ) {
      showNotification(
        paste(
          "This shared view was created with an older application version.",
          "Available settings will be restored where possible."
        ),
        type = "warning",
        duration = 8
      )
    }
    
    updateCheckboxInput(
      session, "omit_zero_values", value = isTRUE(saved$omit_zero_values)
    )
    updateSelectInput(
      session, "plot_type", selected = as.character(saved$plot_type)[1]
    )
    if (!is.null(saved$trajectory_measure)) {
      updateSelectInput(
        session,
        "trajectory_measure",
        selected = as.character(saved$trajectory_measure)[1]
      )
    }
    updateSelectInput(
      session,
      "progression_variable",
      selected = as.character(saved$progression_variable)[1]
    )
    updateSelectInput(
      session, "facet_variable", selected = as.character(saved$facet_variable)[1]
    )
    updateSelectInput(
      session, "color_variable", selected = as.character(saved$color_variable)[1]
    )
    updateSelectInput(
      session, "smoother", selected = as.character(saved$smoother)[1]
    )
    updateSelectInput(
      session, "x_variable", selected = as.character(saved$x_variable)[1]
    )
    updateSelectInput(
      session,
      "second_dimension",
      selected = as.character(saved$second_dimension)[1]
    )
    updateSelectInput(
      session,
      "correlation_color_variable",
      selected = as.character(saved$correlation_color_variable)[1]
    )
    updateCheckboxInput(
      session,
      "correlation_automatic_limits",
      value = isTRUE(saved$correlation_automatic_limits)
    )
    updateNumericInput(
      session,
      "correlation_axis_minimum",
      value = as.numeric(saved$correlation_axis_minimum)[1]
    )
    updateNumericInput(
      session,
      "correlation_axis_maximum",
      value = as.numeric(saved$correlation_axis_maximum)[1]
    )
    updateCheckboxInput(
      session, "show_orthogonal_fit", value = isTRUE(saved$show_orthogonal_fit)
    )
    if (!is.null(saved$correlation_size_minimum)) {
      updateNumericInput(
        session,
        "correlation_size_minimum",
        value = as.numeric(saved$correlation_size_minimum)[1]
      )
    }
    if (!is.null(saved$correlation_size_maximum)) {
      updateNumericInput(
        session,
        "correlation_size_maximum",
        value = as.numeric(saved$correlation_size_maximum)[1]
      )
    }
    updateCheckboxInput(
      session, "log_scale", value = isTRUE(saved$log_scale)
    )
    updateCheckboxInput(
      session,
      "automatic_expression_limits",
      value = isTRUE(saved$automatic_expression_limits)
    )
    updateNumericInput(
      session,
      "expression_minimum",
      value = as.numeric(saved$expression_minimum)[1]
    )
    updateNumericInput(
      session,
      "expression_maximum",
      value = as.numeric(saved$expression_maximum)[1]
    )
    if (!is.null(saved$dot_size_minimum)) {
      updateNumericInput(
        session,
        "dot_size_minimum",
        value = as.numeric(saved$dot_size_minimum)[1]
      )
    }
    if (!is.null(saved$dot_size_maximum)) {
      updateNumericInput(
        session,
        "dot_size_maximum",
        value = as.numeric(saved$dot_size_maximum)[1]
      )
    }
    updateCheckboxInput(
      session, "show_heatmap_counts", value = isTRUE(saved$show_heatmap_counts)
    )
    updateCheckboxInput(
      session, "show_points", value = isTRUE(saved$show_points)
    )
    if (!is.null(saved$violin_display)) {
      updateSelectInput(session, "violin_display", selected = saved$violin_display)
    }
    if (!is.null(saved$violin_facets_per_page)) {
      updateSelectInput(
        session,
        "violin_facets_per_page",
        selected = as.character(saved$violin_facets_per_page)[1]
      )
    }
    if (!is.null(saved$violin_page)) {
      updateSelectInput(
        session,
        "violin_page",
        selected = as.character(saved$violin_page)[1]
      )
    }
    if (!is.null(saved$show_violin_points)) {
      updateCheckboxInput(
        session,
        "show_violin_points",
        value = isTRUE(saved$show_violin_points)
      )
    }
    if (!is.null(saved$facets_per_row)) {
      updateSelectInput(
        session,
        "facets_per_row",
        selected = as.character(saved$facets_per_row)[1]
      )
    }
    if (!is.null(saved$facet_row_height)) {
      updateSelectInput(
        session,
        "facet_row_height",
        selected = as.character(saved$facet_row_height)[1]
      )
    }
    
    saved_filter_fields <- as.character(saved$filter_fields)
    saved_filter_fields <- saved_filter_fields[
      saved_filter_fields %in% filter_fields
    ]
    
    updateSelectizeInput(
      session,
      "filter_fields",
      choices = filter_fields,
      selected = saved_filter_fields,
      server = TRUE
    )
    
    session$onFlushed(function() {
      if (length(saved_filter_fields) > 0) {
        for (field in saved_filter_fields) {
          saved_value <- saved$stacked_filter_values[[field]]
          if (is.null(saved_value)) {
            next
          }
          input_id <- filter_input_id(field)
          if (field %in% numeric_filter_fields) {
            updateSliderInput(
              session,
              input_id,
              value = as.numeric(saved_value)
            )
          } else {
            updateSelectizeInput(
              session,
              input_id,
              choices = field_levels(metadata, field),
              selected = as.character(saved_value),
              server = TRUE
            )
          }
        }
      }
      
      session$onFlushed(function() {
        primary_gene <- as.character(saved$gene)[1]
        if (
          is.na(primary_gene) ||
          !nzchar(primary_gene) ||
          !primary_gene %in% available_genes
        ) {
          showNotification(
            "The primary gene in this shared view is no longer available.",
            type = "error",
            duration = NULL
          )
          pending_bookmark_state(NULL)
          restored_filter_values(NULL)
          restore_plot_pending(FALSE)
          return()
        }
        
        tryCatch(
          {
            primary_data <- read_gene(primary_gene)
            gene_data(primary_data)
            loaded_gene(primary_gene)
            requested_gene(primary_gene)
            updateSelectizeInput(
              session,
              "gene",
              choices = available_genes,
              selected = primary_gene,
              server = TRUE
            )
            gene_status_message(paste0(
              "Loaded ", primary_gene, " (",
              comma(nrow(primary_data)), " observations)."
            ))
            
            restored_plot_type <- as.character(saved$plot_type)[1]
            if (identical(restored_plot_type, "correlation")) {
              comparison_gene <- as.character(saved$comparison_gene)[1]
              if (
                is.na(comparison_gene) ||
                !nzchar(comparison_gene) ||
                !comparison_gene %in% available_genes ||
                identical(comparison_gene, primary_gene)
              ) {
                stop(
                  "The comparison gene in this shared view is unavailable or invalid."
                )
              }
              comparison_data <- read_gene(comparison_gene)
              comparison_gene_data(comparison_data)
              loaded_comparison_gene(comparison_gene)
              updateSelectizeInput(
                session,
                "comparison_gene",
                choices = available_genes,
                selected = comparison_gene,
                server = TRUE
              )
              comparison_gene_status_message(paste0(
                "Loaded ", comparison_gene, " (",
                comma(nrow(comparison_data)), " observations)."
              ))
            }
            
            restore_plot_pending(TRUE)
          },
          error = function(e) {
            showNotification(
              paste(
                "Unable to restore this shared view:",
                conditionMessage(e)
              ),
              type = "error",
              duration = NULL
            )
            pending_bookmark_state(NULL)
            restored_filter_values(NULL)
            restore_plot_pending(FALSE)
          }
        )
      }, once = TRUE)
    }, once = TRUE)
  })
  
  observeEvent(input$reset_defaults, {
    updateSelectizeInput(
      session,
      "filter_fields",
      selected = character(),
      server = TRUE
    )
    updateCheckboxInput(session, "omit_zero_values", value = FALSE)
    updateCheckboxInput(session, "log_scale", value = TRUE)
    updateCheckboxInput(
      session,
      "automatic_expression_limits",
      value = TRUE
    )
    updateNumericInput(session, "dot_size_minimum", value = 1)
    updateNumericInput(session, "dot_size_maximum", value = 9)
    updateCheckboxInput(
      session,
      "show_heatmap_counts",
      value = TRUE
    )
    updateSelectInput(session, "plot_type", selected = "trajectory")
    updateSelectInput(session, "trajectory_measure", selected = "expression")
    updateSelectInput(session, "progression_variable", selected = age_field)
    updateSelectInput(session, "facet_variable", selected = default_facet)
    updateSelectInput(session, "color_variable", selected = default_color)
    updateSelectInput(
      session,
      "correlation_color_variable",
      selected = default_color
    )
    updateCheckboxInput(
      session,
      "correlation_automatic_limits",
      value = TRUE
    )
    updateCheckboxInput(
      session,
      "show_orthogonal_fit",
      value = TRUE
    )
    updateNumericInput(session, "correlation_size_minimum", value = 1.5)
    updateNumericInput(session, "correlation_size_maximum", value = 7)
    updateSelectInput(session, "smoother", selected = "loess")
    updateCheckboxInput(session, "show_points", value = TRUE)
    updateSelectInput(session, "violin_display", selected = "violin")
    updateSelectInput(session, "violin_facets_per_page", selected = "12")
    updateSelectInput(session, "violin_page", selected = "1")
    updateCheckboxInput(session, "show_violin_points", value = TRUE)
    updateSelectInput(session, "facets_per_row", selected = "3")
    updateSelectInput(session, "facet_row_height", selected = "300")
    updateSelectInput(session, "x_variable", selected = default_x)
    updateSelectInput(session, "second_dimension", selected = default_second)
  })
  
  
  plot_settings <- reactiveVal(NULL)
  
  capture_plot_settings <- function() {
    req(gene_data(), loaded_gene(), input$plot_type)
    
    if (identical(input$plot_type, "correlation")) {
      req(comparison_gene_data(), loaded_comparison_gene())
      
      validate(
        need(
          !identical(
            loaded_gene(),
            loaded_comparison_gene()
          ),
          "Choose two different genes for the correlation plot."
        )
      )
    }
    
    list(
      plot_type = input$plot_type,
      trajectory_measure = input$trajectory_measure,
      filters = active_filter_specification(),
      omit_zero_values = isTRUE(input$omit_zero_values),
      log_scale = isTRUE(input$log_scale),
      automatic_expression_limits = isTRUE(
        input$automatic_expression_limits
      ),
      expression_minimum = input$expression_minimum,
      expression_maximum = input$expression_maximum,
      dot_size_minimum = input$dot_size_minimum,
      dot_size_maximum = input$dot_size_maximum,
      show_heatmap_counts = isTRUE(
        input$show_heatmap_counts
      ),
      progression_variable = input$progression_variable,
      facet_variable = input$facet_variable,
      color_variable = input$color_variable,
      correlation_color_variable =
        input$correlation_color_variable,
      correlation_automatic_limits = isTRUE(
        input$correlation_automatic_limits
      ),
      correlation_axis_minimum =
        input$correlation_axis_minimum,
      correlation_axis_maximum =
        input$correlation_axis_maximum,
      show_orthogonal_fit = isTRUE(
        input$show_orthogonal_fit
      ),
      correlation_size_minimum = input$correlation_size_minimum,
      correlation_size_maximum = input$correlation_size_maximum,
      comparison_gene = loaded_comparison_gene(),
      smoother = input$smoother,
      show_points = isTRUE(input$show_points),
      violin_display = if (
        is.null(input$violin_display) || length(input$violin_display) == 0
      ) "violin" else as.character(input$violin_display[[1]]),
      violin_facets_per_page = if (
        is.null(input$violin_facets_per_page) ||
        length(input$violin_facets_per_page) == 0
      ) "12" else as.character(input$violin_facets_per_page[[1]]),
      violin_page = {
        page_value <- suppressWarnings(as.integer(input$violin_page))
        if (length(page_value) == 0 || !is.finite(page_value[[1]])) {
          1L
        } else {
          page_value[[1]]
        }
      },
      show_violin_points = if (is.null(input$show_violin_points)) {
        TRUE
      } else {
        isTRUE(input$show_violin_points)
      },
      facets_per_row = as.integer(input$facets_per_row),
      facet_row_height = as.integer(input$facet_row_height),
      x_variable = input$x_variable,
      second_dimension = input$second_dimension
    )
  }
  
  observeEvent(
    list(
      input$make_plot,
      bookmark_plot_trigger()
    ),
    {
      plot_settings(capture_plot_settings())
      
      session$onFlushed(
        function() {
          bslib::nav_select(
            id = "main_tabs",
            selected = "Plot gene of interest",
            session = session
          )
        },
        once = TRUE
      )
    },
    ignoreInit = TRUE,
    priority = 100
  )
  
  output$violin_page_control <- renderUI({
    req(
      input$plot_type %in% c("trajectory", "violin"),
      gene_data()
    )
    if (
      identical(input$plot_type, "trajectory") &&
      identical(input$facet_variable, "none")
    ) {
      return(NULL)
    }
    data <- apply_filter_specification(
      gene_data(),
      active_filter_specification()
    )
    if (isTRUE(input$omit_zero_values)) {
      data <- data |> filter(summed_counts > 0)
    }
    facet_field <- if (identical(input$plot_type, "trajectory")) {
      input$facet_variable
    } else {
      input$second_dimension
    }
    if (is.null(facet_field) || !facet_field %in% names(data)) return(NULL)
    number_of_facets <- dplyr::n_distinct(data[[facet_field]], na.rm = TRUE)
    per_page_value <- input$violin_facets_per_page
    if (
      is.null(per_page_value) ||
      length(per_page_value) == 0 ||
      is.na(per_page_value[[1]]) ||
      !nzchar(as.character(per_page_value[[1]]))
    ) {
      per_page_value <- "12"
    } else {
      per_page_value <- as.character(per_page_value[[1]])
    }
    per_page <- if (identical(per_page_value, "all")) {
      max(1L, number_of_facets)
    } else {
      parsed_per_page <- suppressWarnings(as.integer(per_page_value))
      if (length(parsed_per_page) == 0 || !is.finite(parsed_per_page)) {
        parsed_per_page <- 12L
      }
      max(1L, parsed_per_page)
    }
    page_count <- max(1L, ceiling(number_of_facets / per_page))
    current <- suppressWarnings(as.integer(input$violin_page))
    if (
      length(current) == 0 ||
      !is.finite(current[[1]]) ||
      current[[1]] < 1L ||
      current[[1]] > page_count
    ) {
      current <- 1L
    } else {
      current <- current[[1]]
    }
    selectInput(
      "violin_page",
      "Page Number",
      choices = setNames(seq_len(page_count), paste0("Page ", seq_len(page_count), " of ", page_count)),
      selected = current
    )
  })
  
  output$plot_size_css <- renderUI({
    settings <- plot_settings()
    if (is.null(settings)) return(NULL)
    uses_facet_layout <- identical(settings$plot_type, "violin") ||
      (identical(settings$plot_type, "trajectory") &&
         !identical(settings$facet_variable, "none"))
    if (!uses_facet_layout) {
      return(tags$style(HTML(
        ".plot-wrapper { overflow: hidden !important; } #expression_plot { height: 100% !important; }"
      )))
    }
    data <- filtered_data()
    facet_field <- if (identical(settings$plot_type, "trajectory")) {
      settings$facet_variable
    } else {
      settings$second_dimension
    }
    number_of_facets <- dplyr::n_distinct(data[[facet_field]], na.rm = TRUE)
    if (settings$plot_type %in% c("trajectory", "violin")) {
      per_page_value <- settings$violin_facets_per_page
      if (is.null(per_page_value) || length(per_page_value) == 0) {
        per_page_value <- "12"
      } else {
        per_page_value <- as.character(per_page_value[[1]])
      }
      per_page <- if (identical(per_page_value, "all")) {
        max(1L, number_of_facets)
      } else {
        parsed_per_page <- suppressWarnings(as.integer(per_page_value))
        if (length(parsed_per_page) == 0 || !is.finite(parsed_per_page[[1]])) {
          parsed_per_page <- 12L
        }
        max(1L, parsed_per_page[[1]])
      }
      number_of_facets <- min(number_of_facets, per_page)
    }
    facets_per_row <- max(1L, as.integer(settings$facets_per_row))
    facet_height <- max(150L, as.integer(settings$facet_row_height))
    number_of_rows <- max(1L, ceiling(number_of_facets / facets_per_row))
    plot_height <- number_of_rows * facet_height
    tags$style(HTML(sprintf(
      ".plot-wrapper { height: calc(90vh - 233px) !important; overflow-y: auto !important; overflow-x: hidden !important; } .plot-wrapper .shiny-spinner-output-container, .plot-wrapper .load-container { height: auto !important; min-height: %dpx !important; overflow: visible !important; } #expression_plot { height: %dpx !important; min-height: %dpx !important; }",
      plot_height, plot_height, plot_height
    )))
  })
  
  filtered_data <- reactive({
    settings <- plot_settings()
    data <- gene_data()
    req(data)
    
    data <- apply_filter_specification(
      data,
      settings$filters
    )
    
    if (settings$omit_zero_values) {
      data <- data |> filter(summed_counts > 0)
    }
    
    validate(need(
      nrow(data) >= 2,
      "Fewer than two observations remain after filtering."
    ))
    
    use_fraction_measure <-
      settings$plot_type %in% c("trajectory", "heatmap", "violin") &&
      identical(settings$trajectory_measure, "fraction")
    
    if (use_fraction_measure) {
      validate(need(
        "fraction_expressing" %in% names(data) &&
          any(is.finite(data$fraction_expressing)),
        "Fraction-expressing data are unavailable for this gene."
      ))
    }
    
    data |>
      mutate(
        plotted_expression = if (use_fraction_measure) {
          fraction_expressing
        } else if (settings$log_scale) {
          log1p(CPM)
        } else {
          CPM
        }
      )
  })
  
  has_defensible_order <- function(data, field) {
    if (!field %in% names(data)) return(FALSE)
    values <- data[[field]]
    is.numeric(values) || identical(field, age_field) ||
      !is.null(value_set_for(field, values))
  }
  
  matrix_graph_cache <- cachem::cache_mem(max_size = 64 * 1024^2, max_age = Inf)
  matrix_graph_cache_key <- function(x_field, y_field, xo, yo, xi, yi) {
    source <- paste(x_field, y_field, as.integer(xo), as.integer(yo),
                    paste(xi, yi, sep = ":", collapse = "|"), sep = "::")
    bytes <- as.integer(charToRaw(enc2utf8(source))); mod <- 2147483629
    h1 <- 0; h2 <- 0
    for (byte in bytes) { h1 <- (h1 * 131 + byte + 1) %% mod; h2 <- (h2 * 137 + byte + 1) %% mod }
    paste0("matrix_", format(h1, scientific=FALSE, trim=TRUE), "_",
           format(h2, scientific=FALSE, trim=TRUE), "_", length(xi))
  }
  build_sparse_rook_graph <- function(xi, yi, xo, yo) {
    n <- length(xi); lookup <- seq_len(n); names(lookup) <- paste(xi, yi, sep=":")
    from <- integer(); to <- integer()
    for (i in seq_len(n)) {
      keys <- character()
      if (xo) keys <- c(keys, paste(xi[i]-1L,yi[i],sep=":"), paste(xi[i]+1L,yi[i],sep=":"))
      if (yo) keys <- c(keys, paste(xi[i],yi[i]-1L,sep=":"), paste(xi[i],yi[i]+1L,sep=":"))
      nbr <- unname(lookup[keys]); nbr <- nbr[!is.na(nbr)]
      if (length(nbr)) { from <- c(from, rep.int(i,length(nbr))); to <- c(to,nbr) }
    }
    if (!length(from)) return(NULL)
    connected <- tabulate(from, nbins=n) > 0
    if (sum(connected) < 3) return(NULL)
    map <- integer(n); map[connected] <- seq_len(sum(connected)); keep <- connected[from] & connected[to]
    from <- map[from[keep]]; to <- map[to[keep]]; totals <- tabulate(from, nbins=sum(connected))
    w <- 1/totals[from]
    list(connected=connected, from=from, to=to, edge_weights=w, weight_total=sum(w))
  }
  moran_rook_permutation <- function(data, x_field, y_field, value_field="expression",
                                     permutations=999L, seed=19050L) {
    m <- data |> transmute(x=as.character(.data[[x_field]]), y=as.character(.data[[y_field]]),
                           value=as.numeric(.data[[value_field]])) |>
      filter(!is.na(x),nzchar(x),!is.na(y),nzchar(y),is.finite(value))
    if (nrow(m)<3 || !is.finite(var(m$value)) || var(m$value)==0)
      return(list(available=FALSE,reason="Matrix autocorrelation requires at least three occupied cells with variation in the selected metric."))
    xo <- has_defensible_order(data,x_field); yo <- has_defensible_order(data,y_field)
    if (!xo && !yo) return(list(available=FALSE,reason="Matrix autocorrelation was not calculated because neither axis has a defined numeric or value-set order."))
    m$xi <- match(m$x,field_levels(data,x_field)); m$yi <- match(m$y,field_levels(data,y_field))
    key <- matrix_graph_cache_key(x_field,y_field,xo,yo,m$xi,m$yi)
    graph <- matrix_graph_cache$get(key,missing=NULL)
    if (is.null(graph)) { graph <- build_sparse_rook_graph(m$xi,m$yi,xo,yo); if(!is.null(graph)) matrix_graph_cache$set(key,graph) }
    if (is.null(graph)) return(list(available=FALSE,reason="Matrix autocorrelation was not calculated because occupied cells do not form a sufficient rook-neighbor graph."))
    z <- m$value[graph$connected]; z <- z-mean(z); den <- sum(z^2); n <- length(z); norm <- n/graph$weight_total
    calc <- function(v) norm*sum(graph$edge_weights*v[graph$from]*v[graph$to])/den
    obs <- calc(z); set.seed(seed)
    perms <- vapply(seq_len(permutations),function(i) calc(sample(z,n,FALSE)),numeric(1))
    list(available=TRUE,statistic=obs,p_value=(1+sum(perms>=obs))/(permutations+1),
         occupied_cells=n,permutations=permutations)
  }
  
  trajectory_progression_values <- function(data, field) {
    if (is.numeric(data[[field]])) return(as.numeric(data[[field]]))
    if (!has_defensible_order(data, field)) return(rep(NA_real_, nrow(data)))
    match(as.character(data[[field]]), field_levels(data, field))
  }
  
  summarized_data <- reactive({
    settings <- plot_settings()
    data <- filtered_data()
    req(settings$x_variable, settings$second_dimension)
    validate(need(
      settings$x_variable != settings$second_dimension,
      "Choose two different plot dimensions."
    ))
    
    data <- factor_field(data, settings$x_variable)
    data <- factor_field(data, settings$second_dimension)
    
    out <- data |>
      filter(
        !is.na(.data[[settings$x_variable]]),
        !is.na(.data[[settings$second_dimension]]),
        is.finite(plotted_expression)
      ) |>
      group_by(across(all_of(c(
        settings$x_variable,
        settings$second_dimension
      )))) |>
      summarise(
        expression = mean(plotted_expression, na.rm = TRUE),
        n_observations = dplyr::n(),
        detection_counts_available = any(valid_detection_counts),
        expressing_cells = if (any(valid_detection_counts)) {
          sum(number_of_cells_expressing[valid_detection_counts], na.rm = TRUE)
        } else NA_real_,
        total_cells = if (any(valid_detection_counts)) {
          sum(number_of_cells[valid_detection_counts], na.rm = TRUE)
        } else NA_real_,
        fraction_expressing = if (
          isTRUE(detection_counts_available) &&
          is.finite(total_cells) && total_cells > 0
        ) expressing_cells / total_cells else NA_real_,
        .groups = "drop"
      )
    
    validate(need(nrow(out) > 0, "No valid groups remain for this plot."))
    out
  })
  
  matrix_autocorrelation <- reactive({
    settings <- plot_settings(); req(settings$plot_type %in% c("heatmap","dot"))
    moran_rook_permutation(summarized_data(),settings$x_variable,settings$second_dimension,
                           permutations=999L,seed=19050L)
  })
  
  dot_fraction_autocorrelation <- reactive({
    settings <- plot_settings()
    req(identical(settings$plot_type, "dot"))
    data <- summarized_data()
    if (
      !"fraction_expressing" %in% names(data) ||
      !any(is.finite(data$fraction_expressing))
    ) {
      return(list(
        available = FALSE,
        reason = "Fraction-expressing Moran's I is unavailable for this gene."
      ))
    }
    moran_rook_permutation(
      data,
      settings$x_variable,
      settings$second_dimension,
      value_field = "fraction_expressing",
      permutations = 999L,
      seed = 19050L
    )
  })
  
  trajectory_statistics <- reactive({
    settings <- plot_settings(); req(identical(settings$plot_type,"trajectory"))
    data <- filtered_data(); field <- settings$progression_variable
    if (!has_defensible_order(data,field)) return(list(available=FALSE,reason="Trajectory statistics were not calculated because the selected progression axis has no defined numeric or value-set order."))
    data$.progression_stat <- trajectory_progression_values(data,field)
    facet_enabled <- !identical(settings$facet_variable,"none")
    data$.facet_stat <- if (facet_enabled) as.character(data[[settings$facet_variable]]) else "All data"
    groups <- split(data,data$.facet_stat,drop=TRUE)
    rows <- lapply(names(groups),function(name){ d<-groups[[name]]; d<-d[is.finite(d$.progression_stat)&is.finite(d$plotted_expression),,drop=FALSE]
    if(nrow(d)<5 || n_distinct(d$.progression_stat)<3 || n_distinct(d$plotted_expression)<2) return(NULL)
    test<-suppressWarnings(cor.test(d$.progression_stat,d$plotted_expression,method="spearman",exact=FALSE))
    data.frame(facet=name,rho=unname(test$estimate),p=test$p.value,n=nrow(d),stringsAsFactors=FALSE) })
    result <- bind_rows(rows)
    if(!nrow(result)) return(list(available=FALSE,reason="No trajectory facet had at least 5 observations, 3 progression levels, and variation in the selected metric."))
    result$adjusted_p <- if(nrow(result)>1) p.adjust(result$p,method="BH") else result$p
    result <- result[order(result$adjusted_p,-abs(result$rho),result$facet),,drop=FALSE]
    list(available=TRUE,results=result,faceted=facet_enabled)
  })
  
  violin_statistics <- reactive({
    settings <- plot_settings()
    req(identical(settings$plot_type, "violin"))
    data <- filtered_data()
    group_field <- settings$x_variable
    facet_field <- settings$second_dimension
    req(group_field, facet_field)
    
    data <- data |>
      filter(
        !is.na(.data[[group_field]]),
        !is.na(.data[[facet_field]]),
        is.finite(plotted_expression)
      )
    groups <- split(
      data,
      as.character(data[[facet_field]]),
      drop = TRUE
    )
    rows <- lapply(names(groups), function(facet_name) {
      facet_data <- groups[[facet_name]]
      group_values <- as.character(facet_data[[group_field]])
      group_counts <- table(group_values)
      eligible_groups <- names(group_counts[group_counts >= 2])
      facet_data <- facet_data[group_values %in% eligible_groups, , drop = FALSE]
      if (
        nrow(facet_data) < 4 ||
        dplyr::n_distinct(facet_data[[group_field]]) < 2 ||
        dplyr::n_distinct(facet_data$plotted_expression) < 2
      ) {
        return(NULL)
      }
      test <- tryCatch(
        stats::kruskal.test(
          facet_data$plotted_expression,
          as.factor(facet_data[[group_field]])
        ),
        error = function(e) NULL
      )
      if (is.null(test) || !is.finite(test$p.value)) return(NULL)
      data.frame(
        facet = facet_name,
        statistic = unname(test$statistic),
        degrees_freedom = unname(test$parameter),
        p = test$p.value,
        n = nrow(facet_data),
        groups = dplyr::n_distinct(facet_data[[group_field]]),
        stringsAsFactors = FALSE
      )
    })
    result <- dplyr::bind_rows(rows)
    if (!nrow(result)) {
      return(list(
        available = FALSE,
        reason = paste(
          "No violin facet had at least two groups with two observations",
          "per group and variation in the selected metric."
        )
      ))
    }
    result$adjusted_p <- stats::p.adjust(result$p, method = "BH")
    result <- result[
      order(result$adjusted_p, -result$statistic, result$facet),
      ,
      drop = FALSE
    ]
    list(
      available = TRUE,
      results = result,
      group_field = group_field,
      facet_field = facet_field
    )
  })
  
  output$plot_statistics_note <- renderText({
    settings <- plot_settings()
    if(settings$plot_type %in% c("heatmap","dot")) {
      r <- matrix_autocorrelation()
      if (!isTRUE(r$available)) return(r$reason)
      format_moran <- function(result, label) {
        paste0(
          label,
          " Moran's I = ",
          formatC(result$statistic, digits = 3, format = "f"),
          ", one-sided permutation p = ",
          format.pval(result$p_value, digits = 3, eps = 1e-3),
          " (", result$occupied_cells, " occupied cells; ",
          result$permutations, " permutations)"
        )
      }
      if (identical(settings$plot_type, "dot")) {
        fraction_result <- dot_fraction_autocorrelation()
        fraction_text <- if (isTRUE(fraction_result$available)) {
          format_moran(fraction_result, "Fraction-expressing")
        } else {
          fraction_result$reason
        }
        return(paste0(
          "Matrix autocorrelation: ",
          format_moran(r, "Expression"),
          "; ",
          fraction_text,
          "."
        ))
      }
      metric_description <- if (
        identical(settings$trajectory_measure, "fraction")
      ) {
        "Fraction-expressing"
      } else {
        "Expression"
      }
      return(paste0(
        "Evidence of ", tolower(metric_description),
        " autocorrelation across the current matrix arrangement: ",
        format_moran(r, metric_description),
        "."
      ))
    }
    if (identical(settings$plot_type, "violin")) {
      r <- violin_statistics()
      if (!isTRUE(r$available)) return(r$reason)
      tab <- r$results
      fmt <- function(row) paste0(
        row$facet,
        " (Kruskal-Wallis chi-squared = ",
        formatC(row$statistic, digits = 2, format = "f"),
        ", df = ",
        formatC(row$degrees_freedom, digits = 0, format = "f"),
        ", BH-adjusted p = ",
        format.pval(row$adjusted_p, digits = 3, eps = 1e-300),
        ")"
      )
      significant <- tab[tab$adjusted_p < 0.05, , drop = FALSE]
      if (!nrow(significant)) {
        return(paste0(
          "No ", r$facet_field,
          " facet showed a significant difference across ",
          r$group_field,
          " groups (BH-adjusted p < 0.05; ",
          nrow(tab), " eligible facets tested)."
        ))
      }
      first <- fmt(significant[1, , drop = FALSE])
      others <- significant[-1, , drop = FALSE]
      shown <- head(others, 4)
      text <- paste0(
        "Strongest group difference among ",
        r$group_field,
        " values: ",
        first,
        "."
      )
      if (nrow(shown)) {
        text <- paste0(
          text,
          " Other significant ",
          r$facet_field,
          " facets: ",
          paste(
            vapply(
              seq_len(nrow(shown)),
              function(i) fmt(shown[i, , drop = FALSE]),
              character(1)
            ),
            collapse = "; "
          ),
          "."
        )
      }
      remaining <- nrow(others) - nrow(shown)
      if (remaining > 0) {
        text <- paste0(
          text,
          " ", remaining,
          " additional significant facets."
        )
      }
      return(text)
    }
    if(!identical(settings$plot_type,"trajectory")) return("")
    r<-trajectory_statistics(); if(!isTRUE(r$available)) return(r$reason); tab<-r$results
    fmt <- function(row) paste0(row$facet," (Spearman rho = ",formatC(row$rho,digits=2,format="f"),
                                ", ",if(r$faceted) "BH-adjusted " else "","p = ",format.pval(row$adjusted_p,digits=3,eps=1e-300),")")
    sig<-tab[tab$adjusted_p<0.05,,drop=FALSE]
    if(!nrow(sig)) return(paste0("No ",if(r$faceted) "facet showed" else "overall data showed",
                                 " a significant monotonic trajectory (",if(r$faceted) "BH-adjusted " else "","p < 0.05; ",nrow(tab)," eligible ",if(r$faceted) "facets" else "test"," tested)."))
    if(!r$faceted) return(paste0("Overall monotonic trajectory: ",fmt(sig[1,,drop=FALSE]),"."))
    first<-fmt(sig[1,,drop=FALSE]); others<-sig[-1,,drop=FALSE]; shown<-head(others,4)
    text<-paste0("Strongest monotonic trajectory: ",first,".")
    if(nrow(shown)) text<-paste0(text," Other significant facets: ",paste(vapply(seq_len(nrow(shown)),function(i) fmt(shown[i,,drop=FALSE]),character(1)),collapse="; "),".")
    remaining<-nrow(others)-nrow(shown); if(remaining>0) text<-paste0(text," ",remaining," additional significant facets.")
    text
  })
  
  orthogonal_fit <- function(x, y) {
    complete <- is.finite(x) & is.finite(y)
    x <- x[complete]
    y <- y[complete]
    if (length(x) < 2) {
      return(NULL)
    }
    covariance <- stats::cov(cbind(x, y))
    decomposition <- eigen(covariance, symmetric = TRUE)
    direction <- decomposition$vectors[, which.max(decomposition$values)]
    center <- c(mean(x), mean(y))
    if (abs(direction[[1]]) < sqrt(.Machine$double.eps)) {
      return(list(type = "vertical", intercept = center[[1]]))
    }
    slope <- direction[[2]] / direction[[1]]
    list(
      type = "line",
      slope = slope,
      intercept = center[[2]] - slope * center[[1]]
    )
  }
  
  correlation_analysis <- reactive({
    settings <- plot_settings()
    req(identical(settings$plot_type, "correlation"))
    primary <- gene_data()
    comparison <- comparison_gene_data()
    req(primary, comparison)
    
    comparison_values <- comparison |>
      transmute(
        sample_id = as.character(sample_id),
        comparison_CPM = as.numeric(CPM),
        comparison_fraction_expressing = as.numeric(fraction_expressing)
      )
    
    data <- primary |>
      inner_join(comparison_values, by = "sample_id")
    data <- apply_filter_specification(data, settings$filters)
    if (settings$omit_zero_values) {
      data <- data |>
        filter(CPM > 0, comparison_CPM > 0)
    }
    use_fraction_measure <- identical(
      settings$trajectory_measure,
      "fraction"
    )
    if (use_fraction_measure) {
      validate(need(
        any(is.finite(data$fraction_expressing)) &&
          any(is.finite(data$comparison_fraction_expressing)),
        "Fraction-expressing data are unavailable for one or both genes."
      ))
    }
    data <- data |>
      mutate(
        primary_expression = if (use_fraction_measure) {
          fraction_expressing
        } else {
          log1p(CPM)
        },
        comparison_expression = if (use_fraction_measure) {
          comparison_fraction_expressing
        } else {
          log1p(comparison_CPM)
        }
      ) |>
      filter(
        is.finite(primary_expression),
        is.finite(comparison_expression)
      )
    
    validate(
      need(nrow(data) >= 3, "At least three matched observations are required."),
      need(
        dplyr::n_distinct(data$primary_expression) >= 2,
        "The primary gene has no variation in the filtered observations."
      ),
      need(
        dplyr::n_distinct(data$comparison_expression) >= 2,
        "The comparison gene has no variation in the filtered observations."
      )
    )
    
    correlation_method <- if (use_fraction_measure) {
      "spearman"
    } else {
      "pearson"
    }
    test <- suppressWarnings(stats::cor.test(
      data$primary_expression,
      data$comparison_expression,
      method = correlation_method,
      exact = if (use_fraction_measure) FALSE else NULL
    ))
    
    list(
      data = data,
      estimate = unname(test$estimate),
      p_value = test$p.value,
      method = correlation_method,
      use_fraction_measure = use_fraction_measure
    )
  })
  
  displayed_facet_levels <- function(data, field, settings) {
    levels <- field_levels(data, field)
    levels <- levels[
      levels %in% as.character(data[[field]])
    ]
    per_page_value <- settings$violin_facets_per_page
    if (is.null(per_page_value) || length(per_page_value) == 0) {
      per_page_value <- "12"
    } else {
      per_page_value <- as.character(per_page_value[[1]])
    }
    per_page <- if (identical(per_page_value, "all")) {
      max(1L, length(levels))
    } else {
      parsed <- suppressWarnings(as.integer(per_page_value))
      if (length(parsed) == 0 || !is.finite(parsed[[1]])) parsed <- 12L
      max(1L, parsed[[1]])
    }
    page_count <- max(1L, ceiling(length(levels) / per_page))
    page <- suppressWarnings(as.integer(settings$violin_page))
    if (length(page) == 0 || !is.finite(page[[1]])) page <- 1L
    page <- min(max(page[[1]], 1L), page_count)
    first <- (page - 1L) * per_page + 1L
    last <- min(page * per_page, length(levels))
    if (length(levels) == 0 || first > last) character() else levels[seq.int(first, last)]
  }
  
  output$expression_plot <- renderPlot({
    settings <- plot_settings()
    data <- filtered_data()
    
    y_label <- if (
      settings$plot_type %in% c("trajectory", "heatmap", "violin") &&
      identical(settings$trajectory_measure, "fraction")
    ) {
      "Fraction of cells expressing"
    } else if (settings$log_scale) {
      "ln(CPM + 1)"
    } else {
      "Counts per million"
    }
    
    automatic_limits <- range(
      data$plotted_expression,
      finite = TRUE
    )
    automatic_limits[[1]] <- 0
    expression_limits <- if (settings$automatic_expression_limits) {
      automatic_limits
    } else {
      c(
        suppressWarnings(as.numeric(settings$expression_minimum)),
        suppressWarnings(as.numeric(settings$expression_maximum))
      )
    }
    validate(need(
      length(expression_limits) == 2 &&
        all(is.finite(expression_limits)) &&
        expression_limits[[1]] < expression_limits[[2]],
      "The expression minimum must be smaller than the expression maximum."
    ))
    
    if (settings$plot_type == "correlation") {
      analysis <- correlation_analysis()
      plot_data <- analysis$data
      color_field <- settings$correlation_color_variable
      color_enabled <- !identical(color_field, "none")
      
      if (color_enabled) {
        plot_data <- factor_field(plot_data, color_field)
        plot_data <- plot_data |>
          filter(!is.na(.data[[color_field]]))
      } else {
        plot_data$correlation_color_group <- "All data"
        color_field <- "correlation_color_group"
      }
      
      observed_limits <- range(
        c(
          plot_data$primary_expression,
          plot_data$comparison_expression
        ),
        finite = TRUE
      )
      observed_limits[[1]] <- min(0, observed_limits[[1]])
      if (isTRUE(analysis$use_fraction_measure)) {
        observed_limits <- c(0, 1)
      }
      correlation_limits <- if (settings$correlation_automatic_limits) {
        observed_limits
      } else {
        c(
          suppressWarnings(as.numeric(settings$correlation_axis_minimum)),
          suppressWarnings(as.numeric(settings$correlation_axis_maximum))
        )
      }
      validate(need(
        length(correlation_limits) == 2 &&
          all(is.finite(correlation_limits)) &&
          correlation_limits[[1]] < correlation_limits[[2]],
        "The correlation-axis minimum must be smaller than the maximum."
      ))
      
      correlation_size_limits <- c(
        suppressWarnings(as.numeric(settings$correlation_size_minimum)),
        suppressWarnings(as.numeric(settings$correlation_size_maximum))
      )
      validate(need(
        length(correlation_size_limits) == 2 &&
          all(is.finite(correlation_size_limits)) &&
          correlation_size_limits[[1]] >= 0 &&
          correlation_size_limits[[1]] < correlation_size_limits[[2]],
        "Min. dot size must be nonnegative and smaller than Max. dot size."
      ))
      
      point_mapping <- aes(
        x = primary_expression,
        y = comparison_expression,
        color = .data[[color_field]]
      )
      use_cell_size <- "number_of_cells" %in% names(plot_data) &&
        any(is.finite(suppressWarnings(as.numeric(plot_data$number_of_cells))))
      if (use_cell_size) {
        plot_data$point_size <- log10(
          pmax(as.numeric(plot_data$number_of_cells), 1)
        )
        point_mapping <- aes(
          x = primary_expression,
          y = comparison_expression,
          color = .data[[color_field]],
          size = point_size
        )
      }
      
      plot <- ggplot(plot_data, point_mapping) +
        geom_point(alpha = 0.68) +
        coord_fixed(
          ratio = 1,
          xlim = correlation_limits,
          ylim = correlation_limits,
          expand = FALSE,
          clip = "on"
        ) +
        labs(
          x = if (isTRUE(analysis$use_fraction_measure)) {
            paste0(loaded_gene(), " fraction expressing")
          } else {
            paste0("ln(", loaded_gene(), " CPM + 1)")
          },
          y = if (isTRUE(analysis$use_fraction_measure)) {
            paste0(loaded_comparison_gene(), " fraction expressing")
          } else {
            paste0("ln(", loaded_comparison_gene(), " CPM + 1)")
          },
          color = if (color_enabled) color_field else NULL,
          size = if (use_cell_size) "log10(number of cells)" else NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          aspect.ratio = 1,
          legend.position = if (color_enabled || use_cell_size) {
            "bottom"
          } else {
            "none"
          }
        )
      
      if (settings$show_orthogonal_fit) {
        fit <- orthogonal_fit(
          plot_data$primary_expression,
          plot_data$comparison_expression
        )
        if (!is.null(fit)) {
          if (identical(fit$type, "vertical")) {
            plot <- plot + geom_vline(
              xintercept = fit$intercept,
              color = "#111827",
              linewidth = 0.9,
              linetype = "solid"
            )
          } else {
            plot <- plot + geom_abline(
              slope = fit$slope,
              intercept = fit$intercept,
              color = "#111827",
              linewidth = 0.9,
              linetype = "solid"
            )
          }
        }
      }
      
      if (color_enabled) {
        plot <- plot + manual_color_scale(plot_data, color_field, "color")
      } else {
        plot <- plot + scale_color_manual(
          values = c("All data" = "#214E68"),
          guide = "none"
        )
      }
      if (use_cell_size) {
        plot <- plot + scale_size_continuous(
          range = correlation_size_limits
        )
      } else {
        plot <- plot + guides(size = "none")
      }
      plot
      
    } else if (settings$plot_type == "trajectory") {
      facet_enabled <- !identical(settings$facet_variable, "none")
      color_enabled <- !identical(settings$color_variable, "none")
      
      if (facet_enabled) {
        data <- factor_field(data, settings$facet_variable)
      }
      if (color_enabled) {
        data <- factor_field(data, settings$color_variable)
      } else {
        data$plot_color_group <- "All data"
      }
      
      progression <- progression_positions(
        data,
        settings$progression_variable,
        maximum_step = 3
      )
      
      data$progression_value <- unname(
        progression$lookup[
          as.character(data[[settings$progression_variable]])
        ]
      )
      
      plot_data <- data |>
        filter(
          is.finite(progression_value),
          is.finite(plotted_expression)
        )
      
      if (facet_enabled) {
        plot_data <- plot_data |>
          filter(!is.na(.data[[settings$facet_variable]]))
        displayed_facets <- displayed_facet_levels(
          plot_data,
          settings$facet_variable,
          settings
        )
        plot_data <- plot_data |>
          filter(
            as.character(.data[[settings$facet_variable]]) %in%
              displayed_facets
          )
        plot_data[[settings$facet_variable]] <- factor(
          as.character(plot_data[[settings$facet_variable]]),
          levels = displayed_facets
        )
      }
      if (color_enabled) {
        plot_data <- plot_data |>
          filter(!is.na(.data[[settings$color_variable]]))
      }
      
      validate(
        need(nrow(plot_data) >= 2, "Too few observations remain."),
        need(
          n_distinct(plot_data$progression_value) >= 2,
          "At least two progression values are required."
        )
      )
      
      color_field <- if (color_enabled) {
        settings$color_variable
      } else {
        "plot_color_group"
      }
      
      plot <- ggplot(
        plot_data,
        aes(
          x = progression_value,
          y = plotted_expression,
          color = .data[[color_field]],
          group = .data[[color_field]]
        )
      )
      
      if (settings$show_points) {
        plot <- plot + geom_point(alpha = 0.42, size = 1.0)
      }
      if (settings$smoother == "loess") {
        plot <- plot + geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = FALSE,
          linewidth = 0.9,
          na.rm = TRUE
        )
      } else if (settings$smoother == "lm") {
        plot <- plot + geom_smooth(
          method = "lm",
          formula = y ~ x,
          se = FALSE,
          linewidth = 0.9,
          na.rm = TRUE
        )
      }
      
      if (color_enabled) {
        plot <- plot + manual_color_scale(
          plot_data,
          settings$color_variable,
          "color"
        )
      } else {
        plot <- plot + scale_color_manual(
          values = c("All data" = "#214E68"),
          guide = "none"
        )
      }
      
      plot <- plot +
        scale_x_continuous(
          breaks = progression$breaks,
          labels = progression$labels,
          limits = progression$limits
        ) +
        scale_y_continuous(
          limits = expression_limits,
          oob = scales::squish
        ) +
        labs(
          x = settings$progression_variable,
          y = y_label,
          color = if (color_enabled) settings$color_variable else NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text = element_text(face = "bold", size = 9),
          legend.position = if (color_enabled) "bottom" else "none"
        )
      
      if (facet_enabled) {
        plot <- plot + facet_wrap(
          vars(.data[[settings$facet_variable]]),
          scales = "fixed",
          ncol = settings$facets_per_row,
          axes = "all_x",
          axis.labels = "all_x",
          drop = TRUE
        )
      }
      
      plot
      
    } else if (settings$plot_type == "heatmap") {
      plot_data <- summarized_data()
      midpoint <- mean(expression_limits)
      plot_data <- plot_data |>
        mutate(
          count_text_color = if_else(
            expression >= midpoint,
            "white",
            "black"
          )
        )
      
      plot <- ggplot(
        plot_data,
        aes(
          x = .data[[settings$x_variable]],
          y = .data[[settings$second_dimension]],
          fill = expression
        )
      ) +
        geom_tile(color = "white", linewidth = 0.2) +
        scale_fill_viridis_c(
          option = "C",
          name = y_label,
          limits = expression_limits,
          oob = scales::squish
        ) +
        labs(x = settings$x_variable, y = settings$second_dimension) +
        theme_minimal(base_size = 11) +
        theme(
          panel.grid = element_blank(),
          axis.text.x = element_text(angle = 55, hjust = 1)
        )
      
      if (settings$show_heatmap_counts) {
        plot <- plot + geom_text(
          aes(
            label = n_observations,
            color = count_text_color
          ),
          size = 3,
          show.legend = FALSE
        ) +
          scale_color_identity()
      }
      
      plot
      
    } else if (settings$plot_type == "dot") {
      plot_data <- summarized_data() |>
        mutate(
          bounded_expression = scales::squish(expression, range = expression_limits),
          size_expression = pmax(
            if (settings$log_scale) bounded_expression else log1p(bounded_expression),
            0.5
          )
        )
      
      dot_size_limits <- c(
        suppressWarnings(as.numeric(settings$dot_size_minimum)),
        suppressWarnings(as.numeric(settings$dot_size_maximum))
      )
      validate(need(
        length(dot_size_limits) == 2 &&
          all(is.finite(dot_size_limits)) &&
          dot_size_limits[[1]] >= 0 &&
          dot_size_limits[[1]] < dot_size_limits[[2]],
        "Min. dot size must be nonnegative and smaller than Max. dot size."
      ))
      
      use_fraction_expressing <-
        "fraction_expressing" %in% names(plot_data) &&
        any(is.finite(plot_data$fraction_expressing))
      
      if (use_fraction_expressing) {
        size_break_fractions <- c(0, 0.25, 0.5, 0.75, 1)
        plot_data <- plot_data |>
          mutate(
            dot_display_size = if_else(
              is.finite(fraction_expressing),
              sqrt(
                dot_size_limits[[1]]^2 +
                  scales::squish(fraction_expressing, c(0, 1)) *
                  (dot_size_limits[[2]]^2 - dot_size_limits[[1]]^2)
              ),
              NA_real_
            )
          )
        size_break_values <- sqrt(
          dot_size_limits[[1]]^2 +
            size_break_fractions *
            (dot_size_limits[[2]]^2 - dot_size_limits[[1]]^2)
        )
        plot <- ggplot(
          plot_data,
          aes(
            x = .data[[settings$x_variable]],
            y = .data[[settings$second_dimension]],
            color = expression,
            size = dot_display_size
          )
        ) +
          geom_point(alpha = 0.9, na.rm = TRUE) +
          scale_size_identity(
            name = "Cells expressing",
            breaks = size_break_values,
            labels = scales::percent(size_break_fractions),
            guide = "legend"
          )
      } else {
        plot <- ggplot(
          plot_data,
          aes(
            x = .data[[settings$x_variable]],
            y = .data[[settings$second_dimension]],
            color = expression,
            size = size_expression
          )
        ) +
          geom_point(alpha = 0.9) +
          scale_size_continuous(
            name = "ln(CPM + 1)",
            range = dot_size_limits,
            oob = scales::squish
          ) +
          labs(caption = paste(
            "Fraction-expressing data are unavailable;",
            "dot size represents expression."
          ))
      }
      
      plot +
        scale_color_viridis_c(
          option = "C",
          name = y_label,
          limits = expression_limits,
          oob = scales::squish
        ) +
        labs(x = settings$x_variable, y = settings$second_dimension) +
        theme_minimal(base_size = 11) +
        theme(axis.text.x = element_text(angle = 55, hjust = 1))
      
    } else {
      data <- factor_field(data, settings$x_variable)
      data <- factor_field(data, settings$second_dimension)
      all_plot_data <- data |>
        filter(
          !is.na(.data[[settings$x_variable]]),
          !is.na(.data[[settings$second_dimension]]),
          is.finite(plotted_expression)
        ) |>
        group_by(across(all_of(c(
          settings$x_variable,
          settings$second_dimension
        )))) |>
        filter(n() >= 2) |>
        ungroup()
      
      validate(need(
        nrow(all_plot_data) >= 2,
        "No violin group has at least two observations."
      ))
      
      facet_levels <- field_levels(
        all_plot_data,
        settings$second_dimension
      )
      facet_levels <- facet_levels[
        facet_levels %in% as.character(all_plot_data[[settings$second_dimension]])
      ]
      per_page_value <- settings$violin_facets_per_page
      if (
        is.null(per_page_value) ||
        length(per_page_value) == 0 ||
        is.na(per_page_value[[1]])
      ) {
        per_page_value <- "12"
      } else {
        per_page_value <- as.character(per_page_value[[1]])
      }
      facets_per_page <- if (identical(per_page_value, "all")) {
        max(1L, length(facet_levels))
      } else {
        parsed_per_page <- suppressWarnings(as.integer(per_page_value))
        if (length(parsed_per_page) == 0 || !is.finite(parsed_per_page[[1]])) {
          parsed_per_page <- 12L
        }
        max(1L, parsed_per_page[[1]])
      }
      page_count <- max(1L, ceiling(length(facet_levels) / facets_per_page))
      selected_page <- suppressWarnings(as.integer(settings$violin_page))
      if (length(selected_page) == 0 || !is.finite(selected_page[[1]])) {
        selected_page <- 1L
      } else {
        selected_page <- selected_page[[1]]
      }
      selected_page <- min(max(selected_page, 1L), page_count)
      first_facet <- (selected_page - 1L) * facets_per_page + 1L
      last_facet <- min(selected_page * facets_per_page, length(facet_levels))
      displayed_facets <- facet_levels[seq.int(first_facet, last_facet)]
      
      plot_data <- all_plot_data |>
        filter(
          as.character(.data[[settings$second_dimension]]) %in% displayed_facets
        )
      plot_data[[settings$second_dimension]] <- factor(
        as.character(plot_data[[settings$second_dimension]]),
        levels = displayed_facets
      )
      
      plot <- ggplot(
        plot_data,
        aes(
          x = .data[[settings$x_variable]],
          y = plotted_expression,
          fill = .data[[settings$x_variable]]
        )
      )
      
      if (identical(settings$violin_display, "boxplot")) {
        plot <- plot + geom_boxplot(
          width = 0.65,
          outlier.shape = NA,
          color = "#214E68",
          alpha = 0.75,
          na.rm = TRUE
        )
      } else {
        plot <- plot + geom_violin(
          scale = "width",
          trim = TRUE,
          color = "#214E68",
          alpha = 0.75,
          na.rm = TRUE
        )
      }
      
      if (isTRUE(settings$show_violin_points)) {
        plot <- plot + geom_jitter(
          width = 0.10,
          color = "black",
          alpha = 0.20,
          size = 0.35
        )
      }
      
      plot +
        stat_summary(
          fun.data = mean_se,
          geom = "errorbar",
          width = 0.18,
          color = "black"
        ) +
        stat_summary(
          fun = mean,
          geom = "point",
          shape = 21,
          fill = "white",
          color = "black",
          size = 2
        ) +
        manual_color_scale(plot_data, settings$x_variable, "fill") +
        facet_wrap(
          vars(.data[[settings$second_dimension]]),
          scales = "fixed",
          ncol = settings$facets_per_row,
          axes = "all_x",
          axis.labels = "all_x",
          drop = TRUE
        ) +
        scale_y_continuous(
          limits = expression_limits,
          oob = scales::squish
        ) +
        labs(x = settings$x_variable, y = y_label) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 9)
        )
    }
  }, res = 96, execOnResize = FALSE)
  
  output$plot_title <- renderText({
    if (is.null(loaded_gene())) {
      return("Retrieve a gene, then generate a plot")
    }
    settings <- plot_settings()
    zero_note <- if (settings$omit_zero_values) {
      "; zero-count observations removed"
    } else {
      ""
    }
    
    if (settings$plot_type == "correlation") {
      analysis <- correlation_analysis()
      paste0(
        loaded_gene(),
        " and ",
        loaded_comparison_gene(),
        if (isTRUE(analysis$use_fraction_measure)) {
          ": fraction-expressing Spearman rho = "
        } else {
          ": Pearson R = "
        },
        formatC(analysis$estimate, digits = 3, format = "f"),
        ", p-value = ",
        format.pval(analysis$p_value, digits = 3, eps = 1e-300),
        zero_note
      )
    } else if (settings$plot_type == "trajectory") {
      facet_note <- if (identical(settings$facet_variable, "none")) {
        ""
      } else {
        paste0(" by ", settings$facet_variable)
      }
      color_note <- if (identical(settings$color_variable, "none")) {
        ""
      } else {
        paste0(" (colored by ", settings$color_variable, ")")
      }
      paste0(
        loaded_gene(),
        if (identical(settings$trajectory_measure, "fraction")) {
          " fraction expressing across "
        } else {
          " across "
        },
        settings$progression_variable,
        facet_note,
        color_note,
        zero_note
      )
    } else {
      paste0(
        loaded_gene(),
        if (
          settings$plot_type %in% c("heatmap", "violin") &&
          identical(settings$trajectory_measure, "fraction")
        ) {
          " fraction expressing by "
        } else {
          " expression by "
        },
        settings$x_variable,
        " and ",
        settings$second_dimension,
        zero_note
      )
    }
  })
}

shinyApp(
  ui = ui,
  server = server,
  enableBookmarking = "url"
)
