test_that("schema v1 migrates to equivalent local and HPC scientific settings", {
  source <- system.file("extdata", "v2.yaml", package = "hotnetR2", mustWork = TRUE)
  root <- withr::local_tempdir(pattern = "schema-v2 migration ")
  file <- file.path(root, "analysis.yaml")
  expect_true(file.copy(source, file))

  migrated <- migrate_analysis_config(file, generate = TRUE)
  expect_true(file.exists(file.path(root, "analysis_schema_v1.yaml")))
  raw <- yaml::read_yaml(migrated)
  expect_identical(as.integer(raw$schema_version), 2L)
  expect_setequal(names(raw$execution$profiles), c("local", "hpc"))

  local <- read_analysis_config(migrated, profile = "local")
  hpc <- read_analysis_config(migrated, profile = "hpc")
  expect_identical(local$execution_profile, "local")
  expect_identical(hpc$execution_profile, "hpc")
  expect_identical(local$ldak$pruning, hpc$ldak$pruning)
  expect_identical(local$networks, hpc$networks)
  expect_identical(local$reference$ldak_exe, "analysis1/ldak6.2.mac")
  expect_identical(hpc$reference$ldak_exe, "software/ldak7.01.linux")
  expect_match(analysis_paths(local)[["ldak_scripts"]], "/scripts/local$")
  expect_match(analysis_paths(hpc)[["ldak_scripts"]], "/scripts/hpc$")

  expected <- c(
    "run_analysis_local.sh", "run_analysis_hpc.sh", "run_analysis.R",
    "LDAK/scripts/local/run_ldak_enhancer.sh",
    "LDAK/scripts/hpc/run_ldak_enhancer.sh",
    "hHotnet/scripts/local/run_hhotnet.sh",
    "hHotnet/scripts/hpc/run_hhotnet.sh"
  )
  expect_true(all(file.exists(file.path(root, expected))))
  testthat::skip_on_os("windows")
  expect_identical(system2("bash", c("-n", shQuote(file.path(root, "run_analysis_local.sh")))), 0L)
  expect_identical(system2("bash", c("-n", shQuote(file.path(root, "run_analysis_hpc.sh")))), 0L)
  expect_identical(system2("bash", c("-n", shQuote(file.path(root, "LDAK/scripts/local/run_ldak_enhancer.sh")))), 0L)
  expect_identical(system2("bash", c("-n", shQuote(file.path(root, "hHotnet/scripts/hpc/run_hhotnet.sh")))), 0L)
})

test_that("schema v2 validates profiles, advanced settings, and minimum versions", {
  source <- system.file("extdata", "v2.yaml", package = "hotnetR2", mustWork = TRUE)
  config <- yaml::read_yaml(source)
  config <- hotnetR2:::legacy_to_schema_v2(config)
  root <- withr::local_tempdir()
  file <- file.path(root, "analysis.yaml")
  yaml::write_yaml(config, file)

  expect_silent(read_analysis_config(file, profile = "local"))
  expect_silent(read_analysis_config(file, profile = "hpc"))
  expect_error(read_analysis_config(file, profile = "missing"), "Unknown execution profile")

  config$unexpected <- TRUE
  yaml::write_yaml(config, file)
  expect_error(read_analysis_config(file), "Unknown schema-v2 field")

  config$unexpected <- NULL
  config$requirements$hotnetR2 <- ">= 999.0.0"
  yaml::write_yaml(config, file)
  expect_error(read_analysis_config(file), "installed version")
  expect_true(hotnetR2:::version_meets_requirement("7.1", ">= 6.3"))
  expect_false(hotnetR2:::version_meets_requirement("6.2", ">= 6.3"))
})

test_that("execution preflight accepts newer LDAK versions with exact reporting", {
  testthat::skip_on_os("windows")
  source <- system.file("extdata", "v2.yaml", package = "hotnetR2", mustWork = TRUE)
  config <- hotnetR2:::legacy_to_schema_v2(yaml::read_yaml(source))
  root <- withr::local_tempdir()
  executable <- file.path(root, "ldak-test")
  writeLines(c("#!/usr/bin/env bash", "echo 'LDAK 7.1'"), executable)
  Sys.chmod(executable, "0755")
  config$execution$profiles$local$platform$os <-
    hotnetR2:::normalise_platform_os(Sys.info()[["sysname"]])
  config$execution$profiles$local$platform$architecture <-
    hotnetR2:::normalise_architecture(R.version$arch)
  config$execution$profiles$local$ldak$executable <- executable
  config$execution$profiles$local$ldak$version <- ">= 6.3"
  file <- file.path(root, "analysis.yaml")
  yaml::write_yaml(config, file)

  report <- check_execution_profile(read_analysis_config(file, profile = "local"))
  expect_true(all(report$status == "OK"))
  expect_identical(report$observed[report$check == "ldak_version"], "7.1")
})

test_that("project paths accept relative symlinks and leading-slash absolute paths", {
  root <- withr::local_tempdir()
  dir.create(file.path(root, "software"))
  target <- file.path(root, "actual-tool")
  writeLines("tool", target)
  expect_true(file.symlink(target, file.path(root, "software", "tool")))
  expect_identical(
    hotnetR2:::resolve_config_path("software/tool", root),
    normalizePath(target)
  )
  expect_identical(hotnetR2:::resolve_config_path(target, root), normalizePath(target))
  expect_error(hotnetR2:::resolve_config_path("~/tool", root), "beginning with /")
  expect_error(hotnetR2:::resolve_config_path("$TOOL", root), "shell expansion")
})

test_that("shared h-HotNet compatibility installation is manifest verified", {
  project <- withr::local_tempdir(pattern = "compat project ")
  expect_identical(check_project_compat(project)$status, "absent")
  first <- initialize_project(project)
  expect_identical(first$status, "ready")
  expect_identical(initialize_project(project)$status, "ready")

  sitecustomize <- file.path(project, "software", "hHotnet", "compat", "sitecustomize.py")
  write("# modified", sitecustomize, append = TRUE)
  expect_identical(check_project_compat(project)$status, "modified_or_incomplete")
  expect_error(initialize_project(project), "update = TRUE")
  expect_identical(initialize_project(project, update = TRUE)$status, "ready")
  backups <- list.files(file.path(project, "software", "hHotnet"),
                        pattern = "^compat-backup-", full.names = TRUE)
  expect_length(backups, 1L)
})
