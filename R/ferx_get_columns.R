#' Inspect column headers of a NONMEM data file
#'
#' Reads only the header row of a NONMEM-format CSV file and prints the column
#' names grouped into required NONMEM columns, optional NONMEM columns, and
#' covariates / user-defined columns. Only the first line is read, so this is
#' fast even on large datasets.
#'
#' @param data One of:
#'   \itemize{
#'     \item A character string path to a NONMEM CSV file.
#'     \item A \code{ferx_fit} object -- uses \code{fit$data_path}.
#'     \item The list returned by \code{\link{ferx_example}} -- uses
#'       \code{x$data}.
#'   }
#'
#' @return A character vector of column names, returned invisibly.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#' ferx_get_columns(ex)        # pass the ferx_example() list directly
#' ferx_get_columns(ex$data)   # or pass the path
#'
#' \dontrun{
#' fit <- ferx_fit(ex$model, ex$data)
#' ferx_get_columns(fit)       # inspect the data used for a fit
#' }
#'
#' @family utilities
#' @export
ferx_get_columns <- function(data) {
  path <- .ferx_resolve_data_path(data)

  if (!file.exists(path)) {
    stop("Data file not found: '", path, "'")
  }

  con <- tryCatch(
    file(path, open = "r"),
    error = function(e) stop("Cannot open '", path, "': ", conditionMessage(e))
  )
  on.exit(close(con), add = TRUE)

  header_line <- tryCatch(
    readLines(con, n = 1L, warn = FALSE),
    error = function(e) stop("Failed to read '", path, "': ", conditionMessage(e))
  )

  if (length(header_line) == 0L || !nzchar(trimws(header_line))) {
    stop("Data file '", path, "' appears to be empty.")
  }

  cols <- tryCatch(
    names(utils::read.csv(
      text = header_line, nrows = 0L,
      stringsAsFactors = FALSE, check.names = FALSE
    )),
    error = function(e) stop(
      "Could not parse column header of '", path, "': ", conditionMessage(e)
    )
  )

  if (length(cols) == 0L) {
    stop("No columns found in header of '", path, "'.")
  }

  nonmem_required <- c("ID", "TIME", "DV", "EVID", "AMT", "CMT")
  nonmem_optional <- c("RATE", "MDV", "II", "SS", "CENS", "OCC")
  nonmem_all      <- c(nonmem_required, nonmem_optional)

  cols_upper <- toupper(cols)
  idx_req <- which(cols_upper %in% nonmem_required)
  idx_opt <- which(cols_upper %in% nonmem_optional)
  idx_cov <- which(!cols_upper %in% nonmem_all)

  n_total  <- length(cols)
  lbl_w    <- 22L
  indent   <- strrep(" ", lbl_w + 1L)
  con_w    <- getOption("width", 80L)

  fmt_group <- function(label, indices) {
    if (length(indices) == 0L) return(invisible(NULL))
    items   <- sprintf("[%d] %s", indices, cols[indices])
    lbl_str <- formatC(label, width = -lbl_w)
    cat(" ", lbl_str, sep = "")
    cur_pos <- lbl_w + 1L
    for (i in seq_along(items)) {
      sep   <- if (i < length(items)) "  " else ""
      token <- paste0(items[i], sep)
      if (cur_pos > lbl_w + 1L && cur_pos + nchar(items[i]) > con_w) {
        cat("\n", indent, sep = "")
        cur_pos <- lbl_w + 1L
      }
      cat(token)
      cur_pos <- cur_pos + nchar(token)
    }
    cat("\n")
  }

  cat(sprintf("Data file: %s\n", path))
  cat(sprintf("%d column%s\n\n", n_total, if (n_total == 1L) "" else "s"))
  fmt_group("Required NONMEM:", idx_req)
  fmt_group("Optional NONMEM:", idx_opt)
  fmt_group("Covariates / other:", idx_cov)
  cat("\n")

  invisible(cols)
}

.ferx_resolve_data_path <- function(data) {
  if (inherits(data, "ferx_fit")) {
    p <- data$data_path
    if (is.null(p) || length(p) != 1L || is.na(p) || !nzchar(p)) {
      stop("The ferx_fit object does not have a valid data_path.")
    }
    return(p)
  }
  if (is.list(data) && !is.null(names(data)) && "data" %in% names(data)) {
    p <- data[["data"]]
    if (!is.character(p) || length(p) != 1L || is.na(p) || !nzchar(p)) {
      stop("The list's '$data' element is not a valid file path string.")
    }
    return(p)
  }
  if (is.character(data) && length(data) == 1L && !is.na(data) && nzchar(data)) {
    return(data)
  }
  stop(
    "`data` must be a file path (character), a ferx_fit object, ",
    "or a ferx_example() list."
  )
}
