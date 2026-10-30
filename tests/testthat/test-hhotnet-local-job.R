test_that("local Python h-HotNet runner uses the configured source and compatibility shim", {
  config <- write_network_test_config(
    available = "[network_1, network_2, network_3, network_4]",
    selected = "[network_1, network_2, network_3, network_4]"
  )
  config$hhotnet$execution_mode <- "local_python"
  config$hhotnet$local_dir <- "/tmp/hierarchical-hotnet"
  config$hhotnet$python_executable <- "python3"
  config$hhotnet$compat_dir <- "compat"
  config$hhotnet$apptainer_image <- "unused.sif"
  lines <- make_hhotnet_job_lines(config, networks = "network_4")
  expect_true(any(grepl("DEFAULT_MODE='local_python'", lines, fixed = TRUE)))
  expect_true(any(grepl("${HHOTNET_DIR}/src/construct_hierarchy.py", lines, fixed = TRUE)))
  expect_true(any(grepl("DEFAULT_NETWORKS=(network_4)", lines, fixed = TRUE)))
  expect_true(any(grepl("Network is already running", lines, fixed = TRUE)))
  script <- tempfile(fileext = ".sh")
  writeLines(lines, script)
  expect_equal(system2("bash", c("-n", script)), 0)
})

test_that("generated runner archives stale results when rebuilding upstream files", {
  testthat::skip_on_os("windows")
  config <- write_network_test_config(
    available = "[network_1]", selected = "[network_1]"
  )
  root <- analysis_paths(config)[["root"]]
  source <- file.path(root, "fake-hhotnet")
  compat <- file.path(root, "compat")
  dir.create(file.path(source, "src"), recursive = TRUE)
  dir.create(compat, recursive = TRUE)
  config$hhotnet$execution_mode <- "local_python"
  config$hhotnet$local_dir <- source
  config$hhotnet$python_executable <- Sys.which("python3")
  config$hhotnet$compat_dir <- compat
  config$hhotnet$score_name <- "nodes_all"
  config$hhotnet$num_permutations <- 1L

  fake <- c(
    "import os, sys",
    "name = os.path.basename(sys.argv[0])",
    "def value(flag): return sys.argv[sys.argv.index(flag) + 1]",
    "def emit(path, text):",
    "    os.makedirs(os.path.dirname(path), exist_ok=True)",
    "    with open(path, 'w') as handle: handle.write(text + '\\n')",
    "if name == 'construct_similarity_matrix.py':",
    "    emit(value('-o'), 'matrix')",
    "    emit(value('-bof'), 'beta')",
    "elif name == 'find_permutation_bins.py': emit(value('-o'), 'bins')",
    "elif name == 'permute_scores.py': emit(value('-o'), 'permuted scores')",
    "elif name == 'construct_hierarchy.py':",
    "    emit(value('-helf'), 'edge hierarchy')",
    "    emit(value('-higf'), 'index hierarchy')",
    "elif name == 'process_hierarchies.py':",
    "    emit(value('-cf'), 'new cluster')",
    "    emit(value('-pf'), 'new plot')",
    "    print('PROCESS_HIERARCHIES')"
  )
  scripts <- c(
    "construct_similarity_matrix.py", "find_permutation_bins.py",
    "permute_scores.py", "construct_hierarchy.py", "process_hierarchies.py"
  )
  for (script_name in scripts) {
    writeLines(fake, file.path(source, "src", script_name))
  }

  paths <- analysis_paths(config)
  dir.create(paths[["hhotnet_data"]], recursive = TRUE)
  writeLines("1\t2", file.path(paths[["hhotnet_data"]], "network_1_edge_list.tsv"))
  writeLines(c("1\tA", "2\tB"),
             file.path(paths[["hhotnet_data"]], "network_1_index_gene.tsv"))
  writeLines(c("A\t2", "B\t1"),
             file.path(paths[["hhotnet_data"]], "nodes_all.tsv"))
  old_results <- file.path(paths[["hhotnet_results"]], "network_1")
  dir.create(old_results, recursive = TRUE)
  writeLines("old cluster", file.path(old_results, "clusters_network_1_nodes_all.tsv"))
  writeLines("old plot", file.path(old_results, "sizes_network_1_nodes_all.pdf"))

  runner <- generate_hhotnet_job(config, "network_1", dry_run = FALSE)
  first <- processx::run("bash", runner, error_on_status = TRUE)
  expect_identical(first$status, 0L)
  expect_identical(
    readLines(file.path(old_results, "clusters_network_1_nodes_all.tsv")),
    "new cluster"
  )
  archives <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  archives <- archives[grepl("^superseded-hhotnet-", basename(archives))]
  expect_length(archives, 1L)
  expect_identical(
    readLines(file.path(archives, "1-network_1", "clusters_network_1_nodes_all.tsv")),
    "old cluster"
  )
  log <- file.path(paths[["logs"]], "network_1_nodes_all.log")
  expect_equal(sum(grepl("PROCESS_HIERARCHIES", readLines(log), fixed = TRUE)), 1L)

  second <- processx::run("bash", runner, error_on_status = TRUE)
  expect_identical(second$status, 0L)
  status_file <- file.path(paths[["hhotnet_summary"]], "status", "network_1.tsv")
  expect_equal(
    sum(grepl("PROCESS_HIERARCHIES", readLines(log), fixed = TRUE)), 1L,
    info = paste(readLines(status_file), collapse = "\n")
  )
  status <- readr::read_tsv(
    status_file,
    show_col_types = FALSE
  )
  expect_true(any(status$stage == "results" & status$status == "archived"))
  expect_true(any(status$stage == "process_hierarchies" & status$status == "reused"))

  hierarchy <- file.path(
    paths[["hhotnet_intermediate"]],
    "network_1_nodes_all", "hierarchy_edge_list_0.tsv"
  )
  expect_true(Sys.setFileTime(hierarchy, Sys.time() + 5))
  third <- processx::run("bash", runner, error_on_status = TRUE)
  expect_identical(third$status, 0L)
  expect_equal(sum(grepl("PROCESS_HIERARCHIES", readLines(log), fixed = TRUE)), 2L)
  status <- readr::read_tsv(status_file, show_col_types = FALSE)
  expect_true(any(
    status$stage == "results" & status$status == "archived" &
      grepl("results predate upstream files", status$detail, fixed = TRUE)
  ))
})
