test_that("type-specific LDAK flanks override the legacy common fallback", {
  config <- write_network_test_config()
  config$ldak$flank_bp <- 500L
  config$ldak$enhancer_flank_bp <- 1000L
  config$ldak$promoter_flank_bp <- 2000L
  config$ldak$result_dirs <- list(enhancer = "enhancer", promoter = "promoter_2k")

  expect_equal(hotnetR2:::ldak_flank_bp(config, "enhancer"), 1000L)
  expect_equal(hotnetR2:::ldak_flank_bp(config, "promoter"), 2000L)
  expect_equal(
    basename(hotnetR2:::ldak_job_specs(config)$promoter$output),
    "promoter_2k"
  )

  enhancer_lines <- hotnetR2:::make_ldak_job_lines(
    config, "enhancer",
    hotnetR2:::ldak_job_specs(config)$enhancer$annotation,
    hotnetR2:::ldak_job_specs(config)$enhancer$output
  )
  promoter_lines <- hotnetR2:::make_ldak_job_lines(
    config, "promoter",
    hotnetR2:::ldak_job_specs(config)$promoter$annotation,
    hotnetR2:::ldak_job_specs(config)$promoter$output
  )
  expect_true(any(grepl("--gene-buffer 1000", enhancer_lines, fixed = TRUE)))
  expect_true(any(grepl("--gene-buffer 2000", promoter_lines, fixed = TRUE)))
})

test_that("legacy common LDAK flank remains supported", {
  config <- write_network_test_config()
  expect_equal(hotnetR2:::ldak_flank_bp(config, "enhancer"), 0L)
  expect_equal(hotnetR2:::ldak_flank_bp(config, "promoter"), 0L)
})
