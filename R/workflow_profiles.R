compat_layer_version <- function() "1.0.0"

compat_manifest_name <- function() ".hotnetR2-compat-manifest.tsv"

normalise_platform_os <- function(x) {
  x <- tolower(as.character(x))
  if (x %in% c("darwin", "macos", "mac")) return("macos")
  if (x %in% c("linux", "gnu/linux")) return("linux")
  x
}

normalise_architecture <- function(x) {
  x <- tolower(as.character(x))
  if (x %in% c("aarch64", "arm64")) return("arm64")
  if (x %in% c("amd64", "x86_64")) return("x86_64")
  x
}

probe_executable_version <- function(executable, timeout = 10000L) {
  result <- tryCatch(
    processx::run(
      executable, "--version", stdout = "|", stderr_to_stdout = TRUE,
      error_on_status = FALSE, timeout = timeout
    ),
    error = function(e) e
  )
  if (inherits(result, "error")) {
    return(list(version = NA_character_, output = conditionMessage(result), status = NA_integer_))
  }
  output <- paste(as.character(result$stdout %||% ""), collapse = "\n")
  match <- regexpr("[0-9]+(?:[.][0-9]+)+", output, perl = TRUE)
  version <- if (match[[1L]] < 0L) NA_character_ else regmatches(output, match)
  status <- as.integer(result$status %||% NA_integer_)
  detail <- trimws(output)
  if (!nzchar(detail)) {
    detail <- if (is.na(status)) "no version output" else
      paste0("exit status ", status, "; no version output")
  }
  list(version = version, output = detail, status = status)
}

#' Check the selected execution profile on the current host
#'
#' Validates OS/architecture, executable permission, and the configured LDAK
#' version policy. A profile can use `version_check: warn` when its LDAK build
#' is known to print an unreliable version banner; the raw output remains in
#' the report and provenance, but it does not block execution.
#'
#' @param config Configuration returned by [read_analysis_config()].
#' @param probe_versions Whether to execute harmless `--version` probes.
#' @return A data frame of preflight checks.
#' @export
check_execution_profile <- function(config, probe_versions = TRUE) {
  profile <- config$execution_profile %||% "legacy"
  runtime <- config$selected_execution %||% list()
  expected_os <- normalise_platform_os(runtime$platform$os %||%
    normalise_platform_os(Sys.info()[["sysname"]]))
  expected_arch <- normalise_architecture(runtime$platform$architecture %||%
    normalise_architecture(R.version$arch))
  observed_os <- normalise_platform_os(Sys.info()[["sysname"]])
  observed_arch <- normalise_architecture(R.version$arch)
  ldak <- resolve_config_path(config$reference$ldak_exe, config$project_root)
  ldak_exists <- file.exists(ldak)
  ldak_executable <- ldak_exists && file.access(ldak, mode = 1L) == 0L
  rows <- list(
    data.frame(
      check = "operating_system", status = if (identical(expected_os, observed_os)) "OK" else "FAIL",
      expected = expected_os, observed = observed_os, detail = profile,
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "architecture", status = if (identical(expected_arch, observed_arch)) "OK" else "FAIL",
      expected = expected_arch, observed = observed_arch, detail = profile,
      stringsAsFactors = FALSE
    ),
    data.frame(
      check = "ldak_executable", status = if (ldak_executable) "OK" else "FAIL",
      expected = ldak,
      observed = if (!ldak_exists) "missing" else if (ldak_executable) "executable" else "not executable",
      detail = Sys.readlink(ldak), stringsAsFactors = FALSE
    )
  )
  if (identical(observed_os, "macos") && ldak_exists && nzchar(Sys.which("xattr"))) {
    quarantine <- processx::run(
      Sys.which("xattr"), c("-p", "com.apple.quarantine", ldak),
      stdout = "|", stderr_to_stdout = TRUE, error_on_status = FALSE
    )
    quarantined <- identical(as.integer(quarantine$status), 0L)
    rows[[length(rows) + 1L]] <- data.frame(
      check = "ldak_quarantine", status = if (quarantined) "FAIL" else "OK",
      expected = "absent",
      observed = if (quarantined) trimws(quarantine$stdout) else "absent",
      detail = if (quarantined) {
        paste0(
          "After independently verifying the download, remove with: xattr -d ",
          "com.apple.quarantine ", shQuote(ldak)
        )
      } else {
        ""
      },
      stringsAsFactors = FALSE
    )
  }
  requirement <- runtime$ldak$version %||% NULL
  if (!is.null(requirement) && isTRUE(probe_versions) && ldak_executable) {
    probe <- probe_executable_version(ldak)
    policy <- runtime$ldak$version_check %||% "enforce"
    advisory <- identical(policy, "warn")
    ok <- !advisory && !is.na(probe$version) &&
      version_meets_requirement(probe$version, requirement, "LDAK version requirement")
    observed <- if (is.na(probe$version)) "unavailable/unreliable" else probe$version
    detail <- if (advisory) {
      paste0(
        "Version comparison skipped because this profile uses version_check: warn. ",
        "Verify the executable source and SHA-256. Raw --version output: ",
        probe$output
      )
    } else {
      probe$output
    }
    rows[[length(rows) + 1L]] <- data.frame(
      check = "ldak_version", status = if (advisory) "WARN" else if (ok) "OK" else "FAIL",
      expected = requirement, observed = observed,
      detail = detail, stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

compat_source_files <- function() {
  source <- system.file("compat", package = "hotnetR2", mustWork = TRUE)
  files <- list.files(source, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                      no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[!grepl("(^|/)__pycache__/|[.]pyc$", files)]
  stats::setNames(files, substring(files, nchar(source) + 2L))
}

compat_manifest <- function(root, version = compat_layer_version()) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE,
                      no.. = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- files[basename(files) != compat_manifest_name()]
  files <- files[!grepl("(^|/)__pycache__/|[.]pyc$", files)]
  relative <- substring(files, nchar(root) + 2L)
  data.frame(
    compat_version = rep(version, length(files)),
    path = relative,
    sha256 = vapply(files, digest::digest, character(1), file = TRUE, algo = "sha256"),
    stringsAsFactors = FALSE
  )
}

#' Check the shared h-HotNet compatibility installation
#'
#' @param project_root Project root containing `software/hHotnet/compat`.
#' @param requirement Minimum compatible asset version using `>= VERSION`.
#' @return A one-row data frame describing the installation.
#' @export
check_project_compat <- function(project_root = ".", requirement = ">= 1.0.0") {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  path <- file.path(project_root, "software", "hHotnet", "compat")
  manifest_file <- file.path(path, compat_manifest_name())
  status <- "absent"
  version <- NA_character_
  detail <- "run hotnetR2::initialize_project(project_root)"
  if (dir.exists(path) && !file.exists(manifest_file)) {
    status <- "unmanifested"
    detail <- "compatibility directory exists without a hotnetR2 manifest"
  } else if (file.exists(manifest_file)) {
    manifest <- tryCatch(
      readr::read_tsv(manifest_file, show_col_types = FALSE),
      error = function(e) NULL
    )
    if (is.null(manifest) || !all(c("compat_version", "path", "sha256") %in% names(manifest)) ||
        !nrow(manifest)) {
      status <- "invalid_manifest"
      detail <- "compatibility manifest is empty or malformed"
    } else {
      version <- as.character(manifest$compat_version[[1L]])
      targets <- file.path(path, manifest$path)
      present <- file.exists(targets)
      hashes <- rep(NA_character_, length(targets))
      hashes[present] <- vapply(targets[present], digest::digest, character(1),
                                file = TRUE, algo = "sha256")
      intact <- present & hashes == manifest$sha256
      if (!all(intact)) {
        status <- "modified_or_incomplete"
        detail <- paste(manifest$path[!intact], collapse = ", ")
      } else if (!version_meets_requirement(version, requirement,
                                            "h-HotNet compatibility requirement")) {
        status <- "too_old"
        detail <- paste0("found ", version, "; requires ", requirement)
      } else {
        status <- "ready"
        detail <- paste0("version ", version)
      }
    }
  }
  data.frame(
    path = path, status = status, version = version, detail = detail,
    stringsAsFactors = FALSE
  )
}

#' Install project-wide h-HotNet compatibility files
#'
#' This is an explicit, idempotent project-initialization step. Package loading
#' and `R CMD INSTALL` never modify a project automatically.
#'
#' @param project_root Project root in which to create
#'   `software/hHotnet/compat`.
#' @param update If `TRUE`, back up and replace a non-current installation.
#' @return The result of `check_project_compat()`, invisibly.
#' @export
initialize_project <- function(project_root = ".", update = FALSE) {
  project_root <- normalizePath(project_root, mustWork = TRUE)
  destination <- file.path(project_root, "software", "hHotnet", "compat")
  current <- check_project_compat(project_root)
  if (identical(current$status, "ready")) return(invisible(current))
  if (dir.exists(destination)) {
    if (!isTRUE(update)) {
      stop(
        "Shared h-HotNet compatibility installation is ", current$status,
        ": ", destination, ". Re-run with update = TRUE to back it up and replace it.",
        call. = FALSE
      )
    }
    backup <- paste0(destination, "-backup-", format(Sys.time(), "%Y%m%d-%H%M%S"))
    if (!file.rename(destination, backup)) {
      stop("Could not back up compatibility directory: ", destination, call. = FALSE)
    }
  }
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  source_files <- compat_source_files()
  for (relative in names(source_files)) {
    target <- file.path(destination, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(source_files[[relative]], target, overwrite = FALSE,
                   copy.mode = TRUE, copy.date = TRUE)) {
      stop("Could not install compatibility file: ", relative, call. = FALSE)
    }
  }
  manifest <- compat_manifest(destination)
  readr::write_tsv(manifest, file.path(destination, compat_manifest_name()))
  result <- check_project_compat(project_root)
  if (!identical(result$status, "ready")) {
    stop("Installed compatibility layer did not validate: ", result$detail, call. = FALSE)
  }
  invisible(result)
}

legacy_to_schema_v2 <- function(config) {
  if (!identical(as.integer(config$schema_version), 1L)) {
    stop("Only schema_version: 1 can be migrated", call. = FALSE)
  }
  old_ldak_exe <- config$reference$ldak_exe
  old_hhotnet <- config$hhotnet
  old_ldak_threads <- config$ldak$threads %||% 8L
  old_slurm <- config$slurm %||% list()
  old_classification <- config$enhancer_classification %||% list(enabled = FALSE)
  if (isTRUE(old_classification$enabled)) {
    old_classification$promoter_flank_bp <- old_classification$promoter_flank_bp %||%
      old_classification$gene_node_flank_bp
    old_classification$promoter_window <- old_classification$promoter_window %||%
      "gene_body"
  }
  old_aliases <- config$reference$promoter_alias_file
  old_alias_columns <- config$reference$promoter_alias_columns
  old_score_cap <- config$scores$infinite_score_cap
  old_aracne_enabled <- config$aracne$enabled %||% FALSE

  config$schema_version <- 2L
  config$features <- list(
    enhancer_classification = isTRUE(old_classification$enabled),
    aracne = isTRUE(old_aracne_enabled)
  )
  old_classification$enabled <- NULL
  config$advanced_settings <- list(
    promoter_aliases = if (is.null(old_aliases)) NULL else list(
      file = old_aliases,
      columns = old_alias_columns %||% list(alias = "alias", gene_symbol = "gene_symbol")
    ),
    enhancer_classification = old_classification,
    score_infinity_cap = old_score_cap
  )
  config$requirements <- list(
    hotnetR2 = ">= 0.3.1",
    hhotnet_compat = ">= 1.0.0"
  )

  config$reference$ldak_exe <- NULL
  config$reference$promoter_alias_file <- NULL
  config$reference$promoter_alias_columns <- NULL
  config$enhancer_classification <- NULL
  config$aracne$enabled <- NULL
  if (!length(config$aracne)) config$aracne <- NULL
  config$scores$infinite_score_cap <- NULL
  config$ldak$threads <- NULL
  config$hhotnet$execution_mode <- NULL
  config$hhotnet$local_dir <- NULL
  config$hhotnet$python_executable <- NULL
  config$hhotnet$compat_dir <- NULL
  config$hhotnet$apptainer_image <- NULL
  config$hhotnet$cpus <- NULL
  if (!is.null(config$interaction_network$stringdb)) {
    config$interaction_network$stringdb$mapping_policy <- NULL
  }
  if (!is.null(config$regulatory$promoter_harmonization)) {
    config$regulatory$promoter_harmonization$jeme_strategy <- NULL
    config$regulatory$promoter_harmonization$hic_strategy <- NULL
  }

  threads <- old_hhotnet$cpus %||% 16L
  ldak_threads <- old_ldak_threads
  local_ldak <- if (identical(old_ldak_exe, "software/ldak")) {
    "software/ldak6.3.mac"
  } else {
    old_ldak_exe
  }
  config$execution <- list(
    default_profile = "local",
    profiles = list(
      local = list(
        platform = list(os = "macos", architecture = "arm64"),
        ldak = list(executable = local_ldak, version = ">= 6.3"),
        hhotnet = list(
          mode = "local_python",
          source_dir = "software/hierarchical-hotnet",
          python_executable = "python3"
        ),
        resources = list(ldak_threads = ldak_threads, hhotnet_cpus = threads),
        scheduler = list(type = "direct")
      ),
      hpc = list(
        platform = list(os = "linux", architecture = "x86_64"),
        environment = list(modules = c(
          "gcc/11.2.0", "bedtools/2.30.0", "openmpi", "r",
          "r-bundle-bioconductor"
        )),
        ldak = list(
          executable = "software/ldak7.01.linux", version = ">= 7.01",
          version_check = "warn"
        ),
        hhotnet = list(mode = "apptainer", image = "software/hhotnet.sif"),
        resources = list(ldak_threads = ldak_threads, hhotnet_cpus = threads),
        scheduler = list(
          type = "slurm",
          time = old_slurm$time %||% "12:00:00",
          memory = old_slurm$memory %||% "16G",
          hhotnet_memory = "128G",
          partition = old_slurm$partition %||% "general-compute",
          qos = old_slurm$qos %||% "general-compute"
        )
      )
    )
  )
  config$slurm <- NULL
  config
}

#' Migrate an analysis YAML from schema v1 to schema v2
#'
#' @param file Existing schema-v1 `analysis.yaml`.
#' @param output Destination YAML; defaults to replacing `file` after backup.
#' @param backup Whether to retain `analysis_schema_v1.yaml` beside the source.
#' @param generate Whether to generate launchers and both profiles' job scripts.
#' @return Normalized path to the schema-v2 YAML.
#' @export
migrate_analysis_config <- function(file, output = file, backup = TRUE,
                                    generate = TRUE) {
  file <- normalizePath(file, mustWork = TRUE)
  config <- yaml::read_yaml(file)
  migrated <- legacy_to_schema_v2(config)
  output <- normalizePath(output, mustWork = FALSE)
  if (isTRUE(backup)) {
    backup_file <- file.path(dirname(file), "analysis_schema_v1.yaml")
    if (!file.exists(backup_file) && !file.copy(file, backup_file, copy.mode = TRUE,
                                                copy.date = TRUE)) {
      stop("Could not create schema-v1 backup: ", backup_file, call. = FALSE)
    }
    analysis_root <- dirname(file)
    for (launcher in c("run_analysis.sh", "run_analysis.R", "run_analysis_linux.sh")) {
      source_launcher <- file.path(analysis_root, launcher)
      if (file.exists(source_launcher)) {
        destination_launcher <- file.path(
          analysis_root,
          paste0(tools::file_path_sans_ext(launcher), "_schema_v1.",
                 tools::file_ext(launcher))
        )
        if (!file.exists(destination_launcher) &&
            !file.rename(source_launcher, destination_launcher)) {
          stop("Could not preserve legacy launcher: ", source_launcher, call. = FALSE)
        }
      }
    }
    for (relative in c(file.path("LDAK", "scripts"),
                       file.path("hHotnet", "scripts"))) {
      script_root <- file.path(analysis_root, relative)
      legacy_files <- list.files(script_root, full.names = TRUE, recursive = FALSE)
      legacy_files <- legacy_files[file.info(legacy_files)$isdir %in% FALSE]
      if (length(legacy_files)) {
        archive <- file.path(script_root, "schema_v1")
        dir.create(archive, recursive = TRUE, showWarnings = FALSE)
        for (legacy_file in legacy_files) {
          destination <- file.path(archive, basename(legacy_file))
          if (!file.exists(destination) && !file.rename(legacy_file, destination)) {
            stop("Could not preserve legacy generated script: ", legacy_file,
                 call. = FALSE)
          }
        }
      }
    }
  }
  yaml::write_yaml(migrated, output)
  read_analysis_config(output, profile = "local")
  read_analysis_config(output, profile = "hpc")
  if (isTRUE(generate)) generate_analysis_scripts(output, overwrite = TRUE)
  output
}

launcher_lines <- function(profile, minimum_version, modules = character()) {
  is_hpc <- identical(profile, "hpc")
  profile_label <- if (is_hpc) "HPC" else "local"
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "",
    "analysis_root=\"$(cd \"$(dirname \"${BASH_SOURCE[0]}\")\" && pwd)\"",
    "config_file=\"$analysis_root/analysis.yaml\"",
    "dispatcher=\"$analysis_root/run_analysis.R\"",
    paste0("profile=", shell_quote(profile)),
    "stage=\"${1:-check}\"",
    paste0("minimum_version=\"${HOTNETR2_MIN_VERSION:-", minimum_version, "}\""),
    "r_script=\"${R_SCRIPT:-Rscript}\"",
    ""
  )
  if (is_hpc && length(modules)) {
    lines <- c(lines,
      "command -v module >/dev/null 2>&1 || { echo 'Environment Modules command is unavailable' >&2; exit 127; }",
      paste0("module load ", paste(vapply(modules, shell_quote, character(1)), collapse = " ")),
      ""
    )
  }
  lines <- c(lines,
    "[[ -f \"$config_file\" ]] || { echo \"Missing $config_file\" >&2; exit 2; }",
    "[[ -f \"$dispatcher\" ]] || { echo \"Missing $dispatcher\" >&2; exit 2; }",
    "command -v \"$r_script\" >/dev/null 2>&1 || { echo \"Rscript executable not found: $r_script\" >&2; exit 127; }",
    ""
  )
  if (is_hpc) {
    lines <- c(lines,
      "if [[ \"$stage\" == \"submit\" ]]; then",
      "  component=\"${2:-}\"",
      "  command -v sbatch >/dev/null 2>&1 || { echo 'sbatch is unavailable' >&2; exit 127; }",
      "  case \"$component\" in",
      "    ldak)",
      "      \"$r_script\" --vanilla \"$dispatcher\" \"$config_file\" \"$profile\" generate-ldak \"$minimum_version\"",
      "      sbatch \"$analysis_root/LDAK/scripts/hpc/run_ldak_enhancer.sh\"",
      "      sbatch \"$analysis_root/LDAK/scripts/hpc/run_ldak_promoter.sh\"",
      "      ;;",
      "    hhotnet|hotnet)",
      "      \"$r_script\" --vanilla \"$dispatcher\" \"$config_file\" \"$profile\" generate-hotnet \"$minimum_version\"",
      "      sbatch \"$analysis_root/hHotnet/scripts/hpc/run_hhotnet.sh\" \"${@:3}\"",
      "      ;;",
      "    *) echo 'Usage: run_analysis_hpc.sh submit {ldak|hhotnet} [network ...]' >&2; exit 2 ;;",
      "  esac",
      "  exit 0",
      "fi",
      "if [[ \"$stage\" == \"all\" || \"$stage\" == \"ldak\" || \"$stage\" == \"hotnet\" ]]; then",
      "  echo 'Heavy HPC stages require explicit submission: run_analysis_hpc.sh submit {ldak|hhotnet}' >&2",
      "  exit 2",
      "fi",
      ""
    )
  }
  c(lines,
    paste0("# Dispatch the ", profile_label, " profile."),
    "exec \"$r_script\" --vanilla \"$dispatcher\" \"$config_file\" \"$profile\" \"$stage\" \"$minimum_version\""
  )
}

#' Generate local and HPC analysis launchers and job scripts
#'
#' @param config_file Schema-v2 analysis YAML.
#' @param overwrite Whether to replace existing generated launchers.
#' @return Named paths to generated launchers and scripts.
#' @export
generate_analysis_scripts <- function(config_file = "analysis.yaml", overwrite = FALSE) {
  config_file <- normalizePath(config_file, mustWork = TRUE)
  raw <- yaml::read_yaml(config_file)
  if (!identical(as.integer(raw$schema_version), 2L)) {
    stop("Profile script generation requires schema_version: 2", call. = FALSE)
  }
  analysis_root <- dirname(config_file)
  requirement <- raw$requirements$hotnetR2 %||% ">= 0.3.1"
  minimum <- version_requirement_minimum(requirement, "requirements.hotnetR2")
  dispatcher <- file.path(analysis_root, "run_analysis.R")
  source_dispatcher <- system.file("scripts", "run_analysis.R", package = "hotnetR2",
                                   mustWork = TRUE)
  outputs <- list(dispatcher = dispatcher)
  if (file.exists(dispatcher) && !isTRUE(overwrite)) {
    stop("Launcher already exists: ", dispatcher, call. = FALSE)
  }
  if (!file.copy(source_dispatcher, dispatcher, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not install analysis dispatcher", call. = FALSE)
  }
  modules <- unlist(raw$execution$profiles$hpc$environment$modules %||% character(),
                    use.names = FALSE)
  for (profile in c("local", "hpc")) {
    launcher <- file.path(analysis_root, paste0("run_analysis_", profile, ".sh"))
    writeLines(launcher_lines(profile, minimum, if (profile == "hpc") modules else character()),
               launcher)
    Sys.chmod(launcher, "0755")
    outputs[[paste0("launcher_", profile)]] <- launcher
    config <- read_analysis_config(config_file, profile = profile)
    outputs[[paste0("ldak_", profile)]] <- generate_ldak_jobs(config, dry_run = FALSE)
    if (profile == "local") generate_ldak_local_runner(config, dry_run = FALSE)
    outputs[[paste0("hhotnet_", profile)]] <- generate_hhotnet_job(config, dry_run = FALSE)
  }
  outputs
}
