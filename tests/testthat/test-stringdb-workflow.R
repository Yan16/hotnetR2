test_that("STRINGdb configuration is additive and requires ARACNe to be disabled", {
  config <- write_network_test_config(
    available = "[network_1, network_2, network_3, network_4]",
    selected = "[network_1, network_2, network_3, network_4]"
  )
  config$aracne <- list(enabled = FALSE)
  config$interaction_network <- list(
    type = "STRINGdb",
    stringdb = list(
      version = "12.0", species = 9606L, score_threshold = 700,
      cache_dir = ".cache/STRINGdb",
      mapping_policy = "legacy_shortest_gene_per_string_id"
    )
  )
  expect_silent(validate_analysis_config(config))
  expect_identical(interaction_network_type(config), "stringdb")
  expect_equal(
    basename(stringdb_cache_paths(config)),
    c(
      "9606.protein.aliases.v12.0.txt.gz",
      "9606.protein.info.v12.0.txt.gz",
      "9606.protein.links.v12.0.txt.gz"
    )
  )

  invalid <- config
  invalid$aracne$enabled <- TRUE
  expect_error(
    validate_analysis_config(invalid),
    "Disable ARACNe"
  )
})

test_that("disabled augmentation cannot advertise Networks 3 or 4", {
  config <- write_network_test_config(
    available = "[network_1, network_2]",
    selected = "[network_1, network_2]"
  )
  config$aracne <- list(enabled = FALSE)
  config$interaction_network <- list(type = "none")
  expect_silent(validate_analysis_config(config))
  config$networks$available <- c("network_1", "network_3")
  config$networks$selected <- c("network_1", "network_3")
  expect_error(validate_analysis_config(config), "require interaction_network")
})

test_that("STRINGdb cache files are declared analysis inputs", {
  config <- write_network_test_config(
    available = "[network_1, network_2, network_3, network_4]",
    selected = "[network_3, network_4]"
  )
  config$aracne <- list(enabled = FALSE)
  config$interaction_network <- list(
    type = "STRINGdb",
    stringdb = list(cache_dir = "string-cache")
  )
  expect_silent(validate_analysis_config(config))
  paths <- hotnetR2:::input_paths(config)
  expect_false("aracne" %in% names(paths))
  expect_true(all(c(
    "stringdb_aliases", "stringdb_info", "stringdb_links"
  ) %in% names(paths)))
})

test_that("cached STRINGdb client maps and reads links without an API request", {
  cache <- withr::local_tempdir()
  paths <- c(
    stringdb_aliases = file.path(cache, "aliases.txt.gz"),
    stringdb_info = file.path(cache, "info.txt.gz"),
    stringdb_links = file.path(cache, "links.txt.gz")
  )
  write_gzip <- function(path, lines) {
    connection <- gzfile(path, open = "wt")
    on.exit(close(connection), add = TRUE)
    writeLines(lines, connection)
  }
  write_gzip(paths[["stringdb_info"]], c(
    "#string_protein_id\tpreferred_name",
    "9606.p1\tGENEA",
    "9606.p2\tGENEB"
  ))
  write_gzip(paths[["stringdb_aliases"]], c(
    "#string_protein_id\talias\tsource",
    "9606.p1\tALIASA\ttest",
    "9606.p2\tALIASB\ttest"
  ))
  write_gzip(paths[["stringdb_links"]], c(
    "protein1 protein2 combined_score",
    "9606.p1 9606.p2 700"
  ))

  client <- hotnetR2:::new_cached_stringdb_client(paths, 700)
  mapped <- client$map(
    data.frame(gene = c("genea", "ALIASB", "missing")), "gene",
    removeUnmappedRows = FALSE
  )
  expect_equal(mapped$STRING_id, c("9606.p1", "9606.p2", NA_character_))
  expect_equal(
    client$get_interactions(c("9606.p1", "9606.p2"))$combined_score,
    700
  )
})

test_that("STRINGdb augmentation is deterministic, restricted, and audited", {
  vertices <- data.frame(
    name = c("A", "LONGNAME", "C", "D", "chr1_10_20"),
    node_type = c(
      "promoter_gene", "promoter_gene", "promoter_gene",
      "promoter_gene", "enhancer"
    ),
    score = c(4, 3, 2, 1, 2), stringsAsFactors = FALSE
  )
  base_edges <- data.frame(
    from = "A", to = "C", source = "JEME", weight = NA_real_,
    stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(
    base_edges, vertices = vertices, directed = FALSE
  )
  observed <- new.env(parent = emptyenv())
  fake <- list(
    map = function(data, column, removeUnmappedRows) {
      observed$mapped_input <- data[[column]]
      data.frame(
        gene = c("A", "LONGNAME", "C", "D"),
        STRING_id = c("p1", "p1", "p2", "p3"),
        stringsAsFactors = FALSE
      )
    },
    get_interactions = function(ids) {
      observed$interaction_ids <- ids
      data.frame(
        from = c("p1", "p1", "p3", "p3", "p2"),
        to = c("p2", "p3", "p1", "p3", "p3"),
        combined_score = c(700, 800, 790, 900, 699),
        stringsAsFactors = FALSE
      )
    }
  )
  result <- hotnetR2:::prepare_stringdb_augmentation(
    graph, fake, score_threshold = 700
  )

  expect_false("chr1_10_20" %in% observed$mapped_input)
  expect_setequal(observed$interaction_ids, c("p1", "p2", "p3"))
  expect_equal(result$conflicts$selected_gene, "A")
  expect_equal(result$summary$interactions_at_threshold, 4L)
  expect_equal(result$summary$self_loops_removed, 1L)
  expect_equal(result$summary$duplicate_gene_pairs_removed, 1L)
  expect_equal(result$summary$overlaps_existing_base, 1L)
  expect_equal(result$summary$string_edges_added, 1L)
  expect_true(any(
    result$edges$Regulator == "A" & result$edges$Target == "C" &
      result$edges$combined_score == 700 & result$edges$overlaps_base
  ))
  graph_edges <- igraph::as_data_frame(result$graph, what = "edges")
  expect_equal(sum(graph_edges$source == "JEME"), 1L)
  expect_equal(sum(graph_edges$source == "STRINGdb"), 1L)
  expect_true(any(
    pmin(graph_edges$from, graph_edges$to) == "A" &
      pmax(graph_edges$from, graph_edges$to) == "D"
  ))
  expect_false(any(graph_edges$from == graph_edges$to))
})
