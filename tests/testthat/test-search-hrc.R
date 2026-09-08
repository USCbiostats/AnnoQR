# annoq-site#78. api-v2 accepts search_hrc on the nine SNP search/count/download
# endpoints and rejects it on /snpAttributes. It is omitted when FALSE rather
# than sent as "false", so existing calls stay byte-identical on the wire.

# Capture the query list httr would send, without any network access.
capture_query <- function(expr) {
  captured <- NULL
  mock <- function(url, query = NULL, ...) {
    captured <<- query
    structure(list(status_code = 200L), class = "response")
  }

  testthat::with_mocked_bindings(
    GET = mock,
    POST = mock,
    .package = "httr",
    try(force(expr), silent = TRUE)
  )
  captured
}

searches <- list(
  function(...) regionQuery("18", 1, 100, ...),
  function(...) rsidsQuery(c("rs1", "rs2"), ...),
  function(...) geneQuery("ZMYND11", ...)
)

counts <- list(
  function(...) countRegionQuery("18", 1, 100, ...),
  function(...) countRsidsQuery(c("rs1", "rs2"), ...),
  function(...) countGeneQuery("ZMYND11", ...)
)

test_that("searches send search_hrc when set", {
  for (fn in searches) {
    expect_equal(capture_query(fn(search_hrc = TRUE))[["search_hrc"]], "true")
  }
})

test_that("searches omit search_hrc by default", {
  for (fn in searches) {
    expect_null(capture_query(fn())[["search_hrc"]])
  }
})

# The fetch_all path switches to POST /<mode>/download. Both paths build from the
# same params list, but a filter that silently stops applying to large result sets
# is exactly the bug worth a test.
test_that("the download path carries search_hrc", {
  for (fn in searches) {
    expect_equal(capture_query(fn(fetch_all = TRUE, search_hrc = TRUE))[["search_hrc"]], "true")
  }
})

test_that("the download path omits search_hrc by default", {
  for (fn in searches) {
    expect_null(capture_query(fn(fetch_all = TRUE))[["search_hrc"]])
  }
})

test_that("counts send search_hrc when set", {
  for (fn in counts) {
    expect_equal(capture_query(fn(search_hrc = TRUE))[["search_hrc"]], "true")
  }
})

test_that("counts omit search_hrc by default", {
  for (fn in counts) {
    expect_null(capture_query(fn())[["search_hrc"]])
  }
})

# /snpAttributes does not accept the parameter, and the SNPWay wrappers call
# snpway.annoq.org rather than api-v2 -- HRC support there is
# Annoq_Overrepr_Workflow#9, not this change.
test_that("excluded functions do not gain the argument", {
  excluded <- c(
    "snpAttributesQuery",
    "snpwayGeneMappingsQuery",
    "snpwayOverrepresentationWorkflowQuery"
  )
  for (name in excluded) {
    expect_false("search_hrc" %in% names(formals(get(name, asNamespace("AnnoQR")))))
  }
})
