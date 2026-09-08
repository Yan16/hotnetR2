test_that("h-HotNet cluster parser preserves headers and every cluster line", {
  file <- tempfile(fileext = ".tsv")
  writeLines(c(
    "# Observed cut height: 0.25",
    "# p-value: 0.04",
    "# Clusters:",
    "A\tB",
    "C"
  ), file)

  result <- read_hhotnet_cluster_result(file, "network_1")
  expect_equal(nrow(result$headers), 3L)
  expect_equal(result$headers$numeric_value[[1]], 0.25)
  expect_equal(result$clusters$cluster_size, c(2L, 1L))
  expect_equal(result$membership$node, c("A", "B", "C"))
  expect_true(result$clusters$is_multi_node[[1]])
  expect_true(result$clusters$is_singleton[[2]])
})

test_that("h-HotNet parser rejects nodes assigned to multiple clusters", {
  file <- tempfile(fileext = ".tsv")
  writeLines(c("A\tB", "A"), file)
  expect_error(read_hhotnet_cluster_result(file, "network_1"), "more than one")
})

test_that("split h-HotNet tables contain only the requested cluster", {
  config <- write_network_test_config()
  result <- list(
    clusters = data.frame(cluster=c(1L,2L),is_multi_node=c(TRUE,TRUE)),
    cluster_nodes = data.frame(network="network_1",cluster=c(1L,1L,2L),node=c("A","B","C")),
    cluster_edges = data.frame(network="network_1",cluster=c(1L,2L),from=c("A","C"),to=c("B","D"))
  )
  manifest_file <- write_hhotnet_cluster_tables(result,config,"network_1")
  manifest <- readr::read_tsv(manifest_file,show_col_types=FALSE)
  first_nodes <- readr::read_tsv(manifest$node_file[[1L]],show_col_types=FALSE)
  second_nodes <- readr::read_tsv(manifest$node_file[[2L]],show_col_types=FALSE)
  expect_equal(first_nodes$cluster,c(1L,1L))
  expect_equal(first_nodes$node,c("A","B"))
  expect_equal(second_nodes$cluster,2L)
  expect_equal(second_nodes$node,"C")
})

test_that("graph construction uses named endpoint columns rather than leading numeric IDs", {
  nodes <- data.frame(name = c("A", "B"), label = c("A", "B"), stringsAsFactors = FALSE)
  edges <- data.frame(id1 = 10L, id2 = 20L, from = "A", to = "B", source = "JEME")
  graph <- hhotnet_igraph_from_tables(nodes, edges)
  expect_equal(igraph::vcount(graph), 2L)
  expect_equal(igraph::ecount(graph), 1L)
})

test_that("permutation summaries parse observed and permuted hierarchy files", {
  config <- write_network_test_config()
  score_dir <- file.path(dirname(config$config_file), "hHotnet", "intermediate", "network_1_nodes_all")
  dir.create(score_dir, recursive = TRUE)
  for (permutation in 0:1) {
    writeLines(c("1\t2\t0.5", "2\t3\t0.1"), file.path(score_dir, paste0("hierarchy_edge_list_", permutation, ".tsv")))
    writeLines(c("1\tA", "2\tB", "3\tC"), file.path(score_dir, paste0("hierarchy_index_gene_", permutation, ".tsv")))
  }
  summary <- summarize_hhotnet_permutation_outputs(config, "network_1")
  expect_equal(summary$permutation, c(0L, 1L))
  expect_equal(summary$hierarchy_kind, c("observed", "permuted"))
  expect_true(all(summary$valid_index))
  expect_true(all(summary$positive_hierarchy_ids))
  expect_true(all(summary$edge_endpoints_in_gene_index))
})

test_that("CX2 declares every exported edge column", {
  skip_if_not_installed("jsonlite")
  nodes <- data.frame(
    name = c("A", "B"), label = c("A", "B"), cy_shape = c("ellipse", "diamond"),
    fill_color = c("#FFFFFF", "#000000"), node_size = c(20, 30), stringsAsFactors = FALSE
  )
  edges <- data.frame(
    id1 = 1L, id2 = 2L, from = "A", to = "B", source = "JEME", weight = 0.5,
    edge_color = "#4C78A8", edge_width = 1.5, edge_line_style = "SOLID",
    stringsAsFactors = FALSE
  )
  file <- tempfile(fileext = ".cx2")
  write_hhotnet_cx2(nodes, edges, file, "test", "network_1", TRUE)
  cx2 <- jsonlite::read_json(file, simplifyVector = FALSE)
  declarations <- cx2_aspect(cx2, "attributeDeclarations")
  expect_true(all(names(edges) %in% names(declarations[[1]]$edges)))
  expect_equal(length(cx2_aspect(cx2, "nodes")), 2L)
  expect_equal(length(cx2_aspect(cx2, "edges")), 1L)
})
