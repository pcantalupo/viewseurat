#!/usr/bin/env Rscript

# Poll the system for running ViewSeurat servers and report their ports.
# A ViewSeurat server is an R interpreter process listening on a TCP port whose
# command line invokes the viewseurat.R launcher.

pacman::p_load('optparse')

##################### Options ########################
opt_parser <- OptionParser(
  option_list=list(
    make_option("--all", default = FALSE, type = "logical", action = "store_true",
                help = "Show every R process listening on TCP, not just ViewSeurat [default: %default]"),
    make_option("--killall", default = FALSE, type = "logical", action = "store_true",
                help = "Send SIGTERM to every ViewSeurat server found (never affects non-ViewSeurat R processes) [default: %default]")
  )
)
opts <- parse_args(opt_parser)
####################################################

# Query lsof for R-interpreter processes listening on TCP, in machine-readable
# (-F) format. Fields: p = PID, c = command, n = network address. The /^R$/
# matches the R interpreter exactly (not "Reminders" etc.), and -a ANDs the
# selectors together (without it, lsof ORs them).
lsof_args <- c("-nP", "-a", "-c", "/^R$/", "-iTCP", "-sTCP:LISTEN", "-F", "pcn")
out <- suppressWarnings(system2("lsof", lsof_args, stdout = TRUE, stderr = FALSE))

# Parse the -F records into one row per listening socket.
servers <- data.frame(
  pid = character(0), addr = character(0), port = character(0),
  stringsAsFactors = FALSE
)
pid <- NA_character_
for (line in out) {
  tag <- substr(line, 1, 1)
  val <- substr(line, 2, nchar(line))
  if (tag == "p") {
    pid <- val
  } else if (tag == "n") {
    port <- sub(".*:", "", val)
    servers <- rbind(servers, data.frame(
      pid = pid, addr = val, port = port, stringsAsFactors = FALSE
    ))
  }
}

if (nrow(servers) == 0) {
  cat("No R servers are listening on a TCP port.\n")
  quit(save = "no", status = 0)
}

# Pull each PID's full command line so we can identify ViewSeurat and extract
# the --seurat argument it was launched with.
get_cmdline <- function(pid) {
  cl <- suppressWarnings(system2("ps", c("-o", "command=", "-p", pid),
                                 stdout = TRUE, stderr = FALSE))
  paste(cl, collapse = " ")
}
servers$cmdline <- vapply(servers$pid, get_cmdline, character(1))

servers$is_viewseurat <- grepl("viewseurat", servers$cmdline, ignore.case = TRUE)

extract_arg <- function(cmdline, flag) {
  m <- regmatches(cmdline, regexpr(paste0(flag, "[= ]+\\S+"), cmdline))
  if (length(m) == 0) {
    return(NA_character_)
  }
  sub(paste0(flag, "[= ]+"), "", m)
}
servers$seurat <- vapply(servers$cmdline, extract_arg, character(1), flag = "--seurat")

# The --seurat path is relative to the directory the server was launched from,
# so resolve it against each process's working directory (reported by lsof).
get_cwd <- function(pid) {
  o <- suppressWarnings(system2("lsof", c("-a", "-p", pid, "-d", "cwd", "-Fn"),
                                stdout = TRUE, stderr = FALSE))
  n <- o[startsWith(o, "n")]
  if (length(n) == 0) {
    return(NA_character_)
  }
  sub("^n", "", n[1])
}
resolve_path <- function(seurat, cwd) {
  if (is.na(seurat)) {
    return(NA_character_)
  }
  if (startsWith(seurat, "/")) {
    return(seurat)
  }
  if (is.na(cwd)) {
    return(seurat)
  }
  normalizePath(file.path(cwd, seurat), mustWork = FALSE)
}
servers$cwd <- vapply(servers$pid, get_cwd, character(1))
servers$seurat <- mapply(resolve_path, servers$seurat, servers$cwd)

if (!opts$all) {
  servers <- servers[servers$is_viewseurat, , drop = FALSE]
}

if (nrow(servers) == 0) {
  cat("No ViewSeurat servers are running. Use --all to list other R servers.\n")
  quit(save = "no", status = 0)
}

report <- servers[, c("pid", "port", "seurat")]
names(report) <- c("PID", "PORT", "SEURAT")
report[is.na(report)] <- ""
print(report, row.names = FALSE, right = FALSE)

if (opts$killall) {
  # Only ever target ViewSeurat-matched processes, even when --all listed
  # other R servers above. This is the safety guarantee: --killall cannot
  # touch an unrelated R session.
  targets <- servers[servers$is_viewseurat, , drop = FALSE]
  if (nrow(targets) == 0) {
    cat("\nNo ViewSeurat servers to kill.\n")
    quit(save = "no", status = 0)
  }
  cat("\nSending SIGTERM to", nrow(targets), "ViewSeurat server(s):\n")
  for (pid in targets$pid) {
    status <- suppressWarnings(system2("kill", pid, stdout = FALSE, stderr = FALSE))
    if (identical(status, 0L)) {
      cat("  killed PID", pid, "\n")
    } else {
      cat("  failed to kill PID", pid, "(exit", status, ")\n")
    }
  }
  quit(save = "no", status = 0)
}

