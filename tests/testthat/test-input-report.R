response_json <- '{
  "gene_list": ["GENEA", "HLA-C.7"],
  "rsId_genes_map": {"rs1": ["GENEA", "HLA-C.7"]},
  "panther_gene_info": {},
  "gene_panther_mapping": {},
  "overrepresentation_results": [],
  "reference_gene_list": null,
  "input_report": {
    "summary": {
      "n_input": 3, "n_annotated": 1, "n_no_record": 1, "n_no_gene": 1,
      "n_records_not_matched_to_input": 0, "n_identifiers_submitted": 2,
      "n_identifiers_mapped_by_panther": 2, "n_identifiers_not_mapped_by_panther": 0,
      "n_identifiers_not_listed_by_panther": 0
    },
    "variants": [
      {"input": "rs1", "variant_key": "rs1", "status": "annotated", "n_genes": 2},
      {"input": "rs2", "variant_key": "rs2", "status": "no_gene", "n_genes": 0},
      {"input": "rs404", "variant_key": null, "status": "no_record", "n_genes": 0}
    ],
    "gene_identifiers": [
      {"identifier": "GENEA", "mapped_by_panther": "yes", "panther_listed_as": "GENEA", "panther_match": "exact"},
      {"identifier": "HLA-C.7", "mapped_by_panther": "yes", "panther_listed_as": "HLA-C", "panther_match": "version_suffix_removed"}
    ]
  }
}'

parse_response <- function(text) jsonlite::fromJSON(text, simplifyVector = FALSE)

test_that("input report becomes a summary and data frames", {
  report <- .build_snpway_input_report(parse_response(response_json))

  expect_equal(report$summary$n_input, 3)
  expect_equal(report$variants$input, c("rs1", "rs2", "rs404"))
  expect_equal(report$variants$status, c("annotated", "no_gene", "no_record"))
  expect_equal(report$variants$variant_key, c("rs1", "rs2", NA))
  expect_type(report$variants$n_genes, "integer")
  expect_equal(report$gene_identifiers$panther_listed_as, c("GENEA", "HLA-C"))
})

test_that("null gene identifier report stays NULL", {
  # As /workflow/gene_mappings sends it: JSON nulls.
  text <- '{"input_report": {
    "summary": {"n_input": 1, "n_identifiers_mapped_by_panther": null},
    "variants": [{"input": "rs1", "variant_key": "rs1", "status": "annotated", "n_genes": 1}],
    "gene_identifiers": null
  }}'

  report <- .build_snpway_input_report(parse_response(text))

  expect_null(report$gene_identifiers)
  expect_true(is.na(report$summary$n_identifiers_mapped_by_panther))
})

test_that("no input report from the backend gives NULL", {
  response <- parse_response(response_json)
  response$input_report <- NULL

  expect_null(.build_snpway_input_report(response))
})

test_that("empty variant list gives a zero-row data frame", {
  response <- parse_response(response_json)
  response$input_report$variants <- list()

  report <- .build_snpway_input_report(response)

  expect_equal(nrow(report$variants), 0)
  expect_equal(names(report$variants), c("input", "variant_key", "status", "n_genes"))
})

test_that("workflow query keeps its existing elements and adds input_report", {
  local_mocked_bindings(.post_snpway_request = function(...) parse_response(response_json))

  result <- snpwayOverrepresentationWorkflowQuery(rsid_list = c("rs1", "rs2", "rs404"))

  expect_equal(
    names(result),
    c("mapping", "panther", "overrepresentation", "csv", "input_report")
  )
  expect_equal(result$input_report$variants$status[[3]], "no_record")
})

test_that("gene mappings query adds input_report", {
  local_mocked_bindings(.post_snpway_request = function(...) parse_response(response_json))

  result <- snpwayGeneMappingsQuery(rsid_list = "rs1")

  expect_equal(names(result), c("mapping", "panther", "input_report"))
})
