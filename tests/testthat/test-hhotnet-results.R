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

test_that("GraphML metadata uses analysis and network names", {
  nodes <- data.frame(name = c("A", "B"), label = c("A", "B"), stringsAsFactors = FALSE)
  edges <- data.frame(from = "A", to = "B", source = "JEME", stringsAsFactors = FALSE)
  graph <- hhotnet_igraph_from_tables(nodes, edges)
  graph <- name_hhotnet_graph(graph, "analysis_E092_Gastric_STRINGdb", "network_4")
  file <- tempfile(fileext = ".graphml")
  igraph::write_graph(graph, file, format = "graphml")
  imported <- igraph::read_graph(file, format = "graphml")

  expect_identical(igraph::graph_attr(imported, "collection_name"),
                   "analysis_E092_Gastric_STRINGdb")
  expect_identical(igraph::graph_attr(imported, "network_name"), "network_4")
  expect_identical(igraph::graph_attr(imported, "name"), "network_4")
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

test_that("CX2 declares attributes and uses mappings without bypasses", {
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
  network_attributes <- cx2_aspect(cx2, "networkAttributes")[[1L]]
  visual <- cx2_aspect(cx2, "visualProperties")[[1L]]
  expect_identical(network_attributes$name, "network_1")
  expect_identical(network_attributes$collection_name, "test")
  expect_identical(network_attributes$network_name, "network_1")
  expect_identical(network_attributes$analysis, "test")
  expect_identical(network_attributes$network, "network_1")
  expect_true(all(names(edges) %in% names(declarations[[1]]$edges)))
  expect_equal(length(cx2_aspect(cx2, "nodes")), 2L)
  expect_equal(length(cx2_aspect(cx2, "edges")), 1L)
  expect_null(cx2_aspect(cx2, "nodeBypasses"))
  expect_null(cx2_aspect(cx2, "edgeBypasses"))
  expect_setequal(
    names(visual$nodeMapping),
    c("NODE_LABEL", "NODE_SHAPE", "NODE_BACKGROUND_COLOR", "NODE_WIDTH", "NODE_HEIGHT")
  )
  expect_setequal(
    names(visual$edgeMapping),
    c("EDGE_LINE_COLOR", "EDGE_WIDTH", "EDGE_LINE_STYLE")
  )
  expect_true(all(vapply(visual$nodeMapping, function(x) identical(x$type, "PASSTHROUGH"), logical(1))))
  expect_true(all(vapply(visual$edgeMapping, function(x) identical(x$type, "PASSTHROUGH"), logical(1))))
  expect_identical(visual$nodeMapping$NODE_SHAPE$definition$attribute, "cy_shape")
  expect_identical(visual$nodeMapping$NODE_BACKGROUND_COLOR$definition$attribute, "fill_color")
  expect_identical(visual$nodeMapping$NODE_WIDTH$definition$attribute, "node_size")
  expect_identical(visual$edgeMapping$EDGE_LINE_COLOR$definition$attribute, "edge_color")
  expect_identical(visual$edgeMapping$EDGE_WIDTH$definition$attribute, "edge_width")
  expect_identical(visual$edgeMapping$EDGE_LINE_STYLE$definition$attribute, "edge_line_style")
})

test_that("Cytoscape import assigns analysis collection and short network names", {
  config <- write_network_test_config(selected = "[network_1]")
  paths <- analysis_paths(config)
  dir.create(paths[["hhotnet_cytoscape"]], recursive = TRUE, showWarnings = FALSE)
  file <- file.path(
    paths[["hhotnet_cytoscape"]],
    paste0(hhotnet_export_prefix(config, "network_1", TRUE), ".cx2")
  )
  writeLines("[]", file)
  observed <- NULL
  testthat::local_mocked_bindings(
    cyrest_import_named_network = function(file, collection_name, network_name, base.url) {
      observed <<- list(file = file, collection_name = collection_name,
                        network_name = network_name, base.url = base.url)
      list(collection_suid = 101L, network_suid = 202L)
    }
  )

  result <- import_hhotnet_cytoscape(config, "network_1", format = "cx2")
  expect_identical(observed$collection_name, "test")
  expect_identical(observed$network_name, "network_1")
  expect_identical(observed$file, normalizePath(file))
  expect_identical(result$collection_name, "test")
  expect_identical(result$network_name, "network_1")
  expect_identical(result$collection_suid, "101")
  expect_identical(result$network_suid, "202")
})
