# The base URL must be overridable: only api-v2-dev.topmed.annoq.org carries
# search_hrc, so a hardcoded default makes the HRC filter untestable.

test_that("the base URL defaults to production", {
  expect_equal(annoq_api_url(), "https://api-v2.annoq.org")
})

test_that("annoq_api_url sets the URL for the session", {
  original <- annoq_api_url()
  on.exit(annoq_api_url(original), add = TRUE)

  annoq_api_url("https://api-v2-dev.topmed.annoq.org")
  expect_equal(annoq_api_url(), "https://api-v2-dev.topmed.annoq.org")
})

test_that("setting the AnnoQ URL does not disturb the SNPWay default", {
  original <- annoq_api_url()
  on.exit(annoq_api_url(original), add = TRUE)

  annoq_api_url("https://example.invalid")
  expect_equal(AnnoQR:::SNPWAY_BASE_URL_DEFAULT, "http://snpway.annoq.org")
})
