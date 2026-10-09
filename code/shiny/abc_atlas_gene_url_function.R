# Create a gene-specific ABC Atlas URL while preserving the remaining atlas state.
#
# Requirements:
#   - Node.js available on PATH
#   - Node package "kiwi-schema" installed
#   - abc_atlas_schema_kiwi.txt available beside app.R, unless schema_path is supplied

create_abc_atlas_gene_url <- function(
    desired_gene,
    starting_gene="Reln",
    starting_url="https://knowledge.brain-map.org/abcatlas#AQEBMTEyV0hYS0E0UkowSTZZVlVOSgACV1MxRkJIRU9OQk9KWjBVU0xRSQADAAQBAQKBQEPMgWLHIgOEfGCChPyeFQAFAQFSZWxuAAAGAQACUmVsbgAFgIDUB4LQb8oGuHcHAiNGRkZGRkYAAwAGR0VORQAHAAgBCQAAAAg1OTJGRTk2NTdDRkY2MTExNTMACUxaSzAxNlFFUlZWSVZXRUJEMTIACgALAW5vbmUAAm5vbmUAAwEEAQACIzAwMDAwMAADyAEABQEBAiMwMDAwMDAAA8gBAAAMAA0BfwAAAAJ%2FAAAAAA4BfwAAAAJ%2BVlVVA3RoZW1lAAAAAgEDAAA%3D",
    schema_path = "abc_atlas_schema_kiwi.txt",
    node_path = Sys.which("node")
    
) {
  inputs <- c(starting_url, starting_gene, desired_gene)

  if (
    any(lengths(list(starting_url, starting_gene, desired_gene)) != 1L) ||
    any(is.na(inputs)) ||
    any(!nzchar(inputs))
  ) {
    stop(
      "starting_url, starting_gene, and desired_gene must each be one non-empty value.",
      call. = FALSE
    )
  }

  if (!file.exists(schema_path)) {
    stop(
      "ABC Atlas Kiwi schema not found: ",
      schema_path,
      call. = FALSE
    )
  }

  if (!nzchar(node_path)) {
    stop(
      "Node.js is required but 'node' was not found on PATH.",
      call. = FALSE
    )
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
    stop(
      "The starting URL contains an empty encoded state.",
      call. = FALSE
    )
  }

  work_dir <- tempfile("abc_atlas_gene_url_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(
    unlink(work_dir, recursive = TRUE, force = TRUE),
    add = TRUE
  )

  encoded_input_path <- file.path(work_dir, "state.txt")
  encoded_output_path <- file.path(work_dir, "state_updated.txt")
  node_script_path <- file.path(work_dir, "convert_gene.js")

  writeLines(
    encoded_state,
    encoded_input_path,
    useBytes = TRUE
  )

  node_script <- c(
    "const kiwi = require('kiwi-schema');",
    "const fs = require('fs');",
    "",
    "const [schemaPath, inputPath, outputPath, startingGene, desiredGene] = process.argv.slice(2);",
    "",
    "function base64ToBytes(base64) {",
    "  const binary = atob(decodeURIComponent(base64.trim()));",
    "  return Uint8Array.from(binary, character => character.codePointAt(0));",
    "}",
    "",
    "function bytesToBase64(bytes) {",
    "  const chunkSize = 0x8000;",
    "  let binary = '';",
    "  for (let start = 0; start < bytes.length; start += chunkSize) {",
    "    binary += String.fromCodePoint(...bytes.subarray(start, start + chunkSize));",
    "  }",
    "  return encodeURIComponent(btoa(binary));",
    "}",
    "",
    "const schema = kiwi.compileSchema(",
    "  kiwi.parseSchema(fs.readFileSync(schemaPath, 'utf8'))",
    ");",
    "",
    "const encodedInput = fs.readFileSync(inputPath, 'utf8');",
    "const payload = schema.decodeExplorePageInitPayload(base64ToBytes(encodedInput));",
    "",
    "let replacements = 0;",
    "",
    "for (const frame of payload.frames || []) {",
    "  for (const gene of frame.genes || []) {",
    "    if (gene.symbol === startingGene) {",
    "      gene.symbol = desiredGene;",
    "      replacements += 1;",
    "    }",
    "  }",
    "",
    "  if (frame.colorBy && frame.colorBy.value === startingGene) {",
    "    frame.colorBy.value = desiredGene;",
    "    replacements += 1;",
    "  }",
    "",
    "  for (const filter of frame.quantitativeFilters || []) {",
    "    if (filter.symbol === startingGene) {",
    "      filter.symbol = desiredGene;",
    "      replacements += 1;",
    "    }",
    "  }",
    "}",
    "",
    "if (replacements === 0) {",
    "  throw new Error(`Starting gene '${startingGene}' was not found in the decoded atlas state.`);",
    "}",
    "",
    "const encodedOutput = bytesToBase64(",
    "  schema.encodeExplorePageInitPayload(payload)",
    ");",
    "",
    "fs.writeFileSync(outputPath, encodedOutput);"
  )

  writeLines(
    node_script,
    node_script_path,
    useBytes = TRUE
  )

  node_output <- system2(
    command = node_path,
    args = c(
      node_script_path,
      normalizePath(schema_path, mustWork = TRUE),
      encoded_input_path,
      encoded_output_path,
      starting_gene,
      desired_gene
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  status <- attr(node_output, "status")

  if (!is.null(status) && status != 0L) {
    stop(
      paste(
        c(
          "ABC Atlas state conversion failed.",
          "Confirm that Node.js and the 'kiwi-schema' Node package are installed.",
          node_output
        ),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  converted_state <- paste0(
    readLines(encoded_output_path, warn = FALSE),
    collapse = ""
  )

  if (!nzchar(converted_state)) {
    stop(
      "The ABC Atlas encoder returned an empty state.",
      call. = FALSE
    )
  }

  paste0(url_prefix, converted_state)
}
