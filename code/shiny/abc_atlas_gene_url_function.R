# Create a gene-specific ABC Atlas URL while preserving the remaining atlas state.
#
# Requirements:
#   - Node.js available on PATH
#   - Node package "kiwi-schema" installed
#   - abc_atlas_schema_kiwi.txt available beside app.R, unless schema_path is supplied

create_abc_atlas_gene_urls <- function(
    desired_genes,
    starting_gene="RELN",
    starting_url="https://knowledge.brain-map.org/abcatlas#AQIBSEZWWDM5QkI2S1QyVExCNlc0QQACREFUU0VBQURNVUxUSVRFU1REM0UxMUUAAwAEAQECgijKOIIW9hUDg26OwYRMMZ4ABQEBUkVMTgAABgEAAlJFTE4ABX%2BUC4yCsj3EBsrbAQcCI0ZGRkZGRgADAAZHRU5FAAcACAEJAAAACEU2Qzg3MzM3NkIyQThBQTA2MwAJVU1TVlhURElBWlRBRktHRTQzVAAKAAsBbm9uZQACbm9uZQADAQQBAAIjMDAwMDAwAAPIAQAFAQECIzAwMDAwMAADyAEAAAwADQF%2FAAAAAn8AAAAADgF%2FAAAAAn5WVVUDdGhlbWUAAAABUU5aTUdPWlhINEVQQ0lVM1YzQQACUlJSNkZOM1dMUjJMOUlWV1hRVAADAAQBAAKMuCC6jSQ8hAOOLokfj87qCAQyTlFUSUU3VEFNUDhQUUFITzRQAAWJT5zAhgPp14lOnMCLmAGzBgAHAAAFAQFQQVg2AAAGAQECQlpZN1hCRzM0VVJBVDFJVTNDNwAACDJDS09RSUsyTjlJRTdUNFdJOEYACUZIVVVLMTlBNVhKN0ZDOVlZVDQACgALAW5vbmUAAm5vbmUAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAMAA0BfwAAAAJ%2FAAAAAA4BfwAAAAJ%2BVlVVA3RoZW1lAAAAAgMDAAA%3D",
    schema_path = "abc_atlas_schema_kiwi.txt",
    node_path = Sys.which("node")
    
) {
  desired_genes <- as.character(desired_genes)
  
  if (
    length(starting_url) != 1L ||
    length(starting_gene) != 1L ||
    is.na(starting_url) ||
    is.na(starting_gene) ||
    !nzchar(starting_url) ||
    !nzchar(starting_gene)
  ) {
    stop(
      "starting_url and starting_gene must each be one non-empty value.",
      call. = FALSE
    )
  }
  
  if (
    length(desired_genes) == 0L ||
    anyNA(desired_genes) ||
    any(!nzchar(desired_genes))
  ) {
    stop(
      "desired_genes must contain at least one non-empty, non-missing gene symbol.",
      call. = FALSE
    )
  }
  
  if (!file.exists(schema_path)) {
    stop("ABC Atlas Kiwi schema not found: ", schema_path, call. = FALSE)
  }
  
  if (!nzchar(node_path)) {
    stop("Node.js is required but 'node' was not found on PATH.", call. = FALSE)
  }
  
  hash_position <- regexpr("#", starting_url, fixed = TRUE)
  
  if (hash_position < 1L) {
    stop(
      "The starting URL does not contain an encoded state after '#'.",
      call. = FALSE
    )
  }
  
  url_prefix <- substr(starting_url, 1L, hash_position)
  encoded_state <- substr(
    starting_url,
    hash_position + 1L,
    nchar(starting_url)
  )
  
  if (!nzchar(encoded_state)) {
    stop("The starting URL contains an empty encoded state.", call. = FALSE)
  }
  
  work_dir <- tempfile("abc_atlas_gene_urls_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  
  encoded_input_path <- file.path(work_dir, "state.txt")
  genes_input_path <- file.path(work_dir, "genes.json")
  urls_output_path <- file.path(work_dir, "urls.json")
  node_script_path <- file.path(work_dir, "convert_genes.js")
  
  writeLines(encoded_state, encoded_input_path, useBytes = TRUE)
  
  jsonlite::write_json(
    desired_genes,
    genes_input_path,
    auto_unbox = FALSE
  )
  
  node_script <- c(
    "const kiwi = require('kiwi-schema');",
    "const fs = require('fs');",
    "",
    "const [schemaPath, inputPath, genesPath, outputPath, startingGene, urlPrefix] = process.argv.slice(2);",
    "",
    "function base64ToBytes(base64) {",
    "  return Uint8Array.from(Buffer.from(decodeURIComponent(base64.trim()), 'base64'));",
    "}",
    "",
    "function bytesToBase64(bytes) {",
    "  return encodeURIComponent(Buffer.from(bytes).toString('base64'));",
    "}",
    "",
    "function replaceGene(payload, fromGene, toGene) {",
    "  let replacements = 0;",
    "",
    "  for (const frame of payload.frames || []) {",
    "    for (const gene of frame.genes || []) {",
    "      if (gene.symbol === fromGene) {",
    "        gene.symbol = toGene;",
    "        replacements += 1;",
    "      }",
    "    }",
    "",
    "    if (frame.colorBy && frame.colorBy.value === fromGene) {",
    "      frame.colorBy.value = toGene;",
    "      replacements += 1;",
    "    }",
    "",
    "    for (const filter of frame.quantitativeFilters || []) {",
    "      if (filter.symbol === fromGene) {",
    "        filter.symbol = toGene;",
    "        replacements += 1;",
    "      }",
    "    }",
    "  }",
    "",
    "  return replacements;",
    "}",
    "",
    "const schema = kiwi.compileSchema(",
    "  kiwi.parseSchema(fs.readFileSync(schemaPath, 'utf8'))",
    ");",
    "",
    "const encodedInput = fs.readFileSync(inputPath, 'utf8');",
    "const desiredGenes = JSON.parse(fs.readFileSync(genesPath, 'utf8'));",
    "const basePayload = schema.decodeExplorePageInitPayload(base64ToBytes(encodedInput));",
    "",
    "const testPayload = structuredClone(basePayload);",
    "if (replaceGene(testPayload, startingGene, startingGene) === 0) {",
    "  throw new Error(`Starting gene '${startingGene}' was not found in the decoded atlas state.`);",
    "}",
    "",
    "const urls = desiredGenes.map((desiredGene) => {",
    "  const payload = structuredClone(basePayload);",
    "  replaceGene(payload, startingGene, desiredGene);",
    "  const encodedOutput = bytesToBase64(",
    "    schema.encodeExplorePageInitPayload(payload)",
    "  );",
    "  return urlPrefix + encodedOutput;",
    "});",
    "",
    "fs.writeFileSync(outputPath, JSON.stringify(urls));"
  )
  
  writeLines(node_script, node_script_path, useBytes = TRUE)
  
  node_output <- system2(
    command = node_path,
    args = c(
      node_script_path,
      normalizePath(schema_path, mustWork = TRUE),
      encoded_input_path,
      genes_input_path,
      urls_output_path,
      starting_gene,
      url_prefix
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  status <- attr(node_output, "status")
  
  if (!is.null(status) && status != 0L) {
    stop(
      paste(
        c(
          "ABC Atlas batch state conversion failed.",
          "Confirm that Node.js and the 'kiwi-schema' Node package are installed.",
          node_output
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }
  
  if (!file.exists(urls_output_path)) {
    stop("The ABC Atlas batch encoder did not create an output file.", call. = FALSE)
  }
  
  converted_urls <- jsonlite::read_json(
    urls_output_path,
    simplifyVector = TRUE
  )
  
  converted_urls <- as.character(converted_urls)
  
  if (
    length(converted_urls) != length(desired_genes) ||
    anyNA(converted_urls) ||
    any(!nzchar(converted_urls))
  ) {
    stop("The ABC Atlas batch encoder returned invalid output.", call. = FALSE)
  }
  
  names(converted_urls) <- desired_genes
  converted_urls
}

# Backward-compatible wrapper for code that still requests one gene at a time.
create_abc_atlas_gene_url <- function(
    desired_gene,
    starting_gene = "RELN",
    starting_url,
    schema_path = "abc_atlas_schema_kiwi.txt",
    node_path = Sys.which("node")
) {
  unname(create_abc_atlas_gene_urls(
    desired_genes = desired_gene,
    starting_gene = starting_gene,
    starting_url = starting_url,
    schema_path = schema_path,
    node_path = node_path
  ))
}
