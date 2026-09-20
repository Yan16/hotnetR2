test_that("manual mode shell ends at similarity and skips every permutation command", {
  config <- write_network_test_config()
  config$hhotnet$analysis_mode <- "manual_delta"
  config$hhotnet$execution_mode <- "local_python"
  config$hhotnet$local_dir <- "/tmp/hhn"
  config$hhotnet$apptainer_image <- "unused.sif"
  config$hhotnet$score_name <- "nodes_all"
  lines <- make_hhotnet_job_lines(config)
  expect_true(any(grepl("run_hhotnet construct_similarity_matrix.py", lines, fixed = TRUE)))
  expect_false(any(grepl("permute_scores.py|find_permutation_bins.py|construct_hierarchy.py|process_hierarchies.py", lines)))
  expect_true(any(grepl("release_lock", tail(lines, 3), fixed = TRUE)))
  script <- tempfile(fileext = ".sh")
  writeLines(lines, script)
  expect_equal(system2("bash", c("-n", script)), 0L)
  expect_error(extract_hhotnet_delta_clusters(config, c(0.1, NA)), "positive finite")
  expect_error(extract_hhotnet_delta_clusters(config, 0), "positive finite")
  expect_error(extract_hhotnet_delta_clusters(config, 0.1, min_cluster_size = 0), "positive integer")
})

test_that("delta exports contain only within-cluster edges and support singleton-only cuts", {
  config <- write_network_test_config(selected = "[network_1]")
  config$hhotnet$execution_mode <- "local_python"
  config$hhotnet$local_dir <- tempdir()
  config$hhotnet$score_name <- "nodes_all"
  paths <- analysis_paths(config)
  data <- paths[["hhotnet_data"]]
  intermediate <- file.path(paths[["hhotnet_intermediate"]], "network_1")
  dir.create(data, recursive = TRUE)
  dir.create(intermediate, recursive = TRUE)
  writeLines("fixture matrix", file.path(intermediate, "similarity_matrix.h5"))
  writeLines(c("1\tA", "2\tB", "3\tC"), file.path(data, "network_1_index_gene.tsv"))
  writeLines(c("A\t3", "B\t2", "C\t1"), file.path(data, "nodes_all.tsv"))
  readr::write_tsv(tibble::tibble(gene = c("A", "B", "C"), node_type = "promoter_gene",
    node_sources = "JEME", score = 3:1, LRT_P_Perm = c(.001, .01, .1)), file.path(data, "nodes_all_info.tsv"))
  readr::write_tsv(tibble::tibble(from = c("A", "B"), to = c("B", "C"), source = "JEME"),
                   file.path(data, "network_1_edge_list_full.tsv"))
  hierarchy_calls <- 0L
  local_mocked_bindings(delta_python = function(config, script, args) {
    arg <- function(flag) args[[match(flag, args) + 1L]]
    if (basename(script) == "construct_hierarchy.py") {
      hierarchy_calls <<- hierarchy_calls + 1L
      writeLines("1\t2\t0.2", arg("-helf"))
      writeLines(c("1\tA", "2\tB", "3\tC"), arg("-higf"))
    } else {
      writeLines(c("# Significance: not evaluated", if (as.numeric(arg("--delta")) <= .2) "A\tB" else c("A", "B"), "C"), arg("--output"))
    }
  })
  manifest <- extract_hhotnet_delta_clusters(config, c(.1, .5))
  expect_equal(manifest$nodes, c(2L, 0L))
  expect_equal(manifest$edges, c(1L, 0L))
  expect_equal(manifest$largest_cluster, c(2L, 1L))
  expect_equal(hierarchy_calls, 1L)
  extract_hhotnet_delta_clusters(config, .2)
  expect_equal(hierarchy_calls, 1L)
  cx <- jsonlite::read_json(manifest$cx2_file[[1]])
  expect_length(cx2_aspect(cx, "nodes"), 2L)
  expect_length(cx2_aspect(cx, "edges"), 1L)
  empty <- jsonlite::read_json(manifest$cx2_file[[2]])
  expect_length(cx2_aspect(empty, "nodes"), 0L)
})
