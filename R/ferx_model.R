#' Create a ferx_model object
#'
#' Constructs a \code{ferx_model} S3 object that bundles a \code{.ferx} model
#' file path with an optional data path. This is the entry point for pipe-based
#' workflows. Both the model file and data path are validated at construction
#' time (the data path may be omitted and supplied later to
#' \code{\link{ferx_fit}}).
#'
#' \strong{Argument order:} \code{data} comes first so a data object can flow
#' naturally into a pipeline:
#'
#' \preformatted{
#' ex$data |> ferx_model(ex$model) |> ferx_fit() |> summary()
#' }
#'
#' This is a change from earlier versions (where \code{model} was the first
#' argument). Old positional calls of the form \code{ferx_model("pk.ferx")} or
#' \code{ferx_model("pk.ferx", "data.csv")} are detected by the \code{.ferx}
#' extension on what is now the \code{data} slot and silently rewritten with a
#' deprecation warning; this auto-correction will be removed in a future
#' release. Calls that name \code{data} explicitly
#' (\code{ferx_model("pk.ferx", data = "data.csv")}) keep working unchanged
#' because R matches \code{data =} by name first and the remaining positional
#' argument falls into the \code{model} slot.
#'
#' All fit options (\code{method}, \code{covariance}, \code{threads},
#' \code{settings}, ...) can still be passed directly to \code{ferx_fit()} in
#' the pipe - the \code{ferx_model} object only carries the file paths. See
#' \code{\link{ferx_fit}} for the full list of options and post-fit outputs.
#'
#' \strong{Scaffold mode.} Passing \code{template =} (or \code{print = TRUE})
#' switches to scaffold mode: a new \code{.ferx} file is written to \code{path}
#' from a built-in skeleton and wrapped in a \code{ferx_model} object, so the
#' result pipes straight into \code{\link{ferx_fit}}. \code{print = TRUE} prints
#' the skeleton to the console and writes nothing (returns \code{NULL}). This
#' replaces the former \code{ferx_model_new()}.
#'
#' @param data Optional path to a NONMEM-format CSV data file. Can be omitted
#'   here and supplied later to \code{\link{ferx_fit}}. When omitted, it defaults
#'   to the dataset declared in the model file's \code{[data]} block
#'   (\code{path = ...}) if the model file has one.
#' @param model Path to an existing \code{.ferx} model file (wrap mode). The
#'   file must exist. Leave \code{NULL} in scaffold mode.
#' @param template Scaffold mode. One of \code{"1cpt_oral"} (default),
#'   \code{"1cpt_iv"}, \code{"2cpt_oral"}, \code{"2cpt_iv"}, or \code{"ode"}.
#'   Supplying it (or \code{print = TRUE}) triggers scaffold mode.
#' @param path Scaffold mode. Path for the new \code{.ferx} file. Required
#'   unless \code{print = TRUE}. Must not already exist unless
#'   \code{overwrite = TRUE}.
#' @param overwrite Scaffold mode. Logical. Overwrite \code{path} if it already
#'   exists? Defaults to \code{FALSE}.
#' @param edit Scaffold mode. Logical. Open the new file in an editor after
#'   writing? Defaults to \code{TRUE}.
#' @param print Scaffold mode. Logical. If \code{TRUE}, print the skeleton to
#'   the console instead of writing a file (no file created, no editor opened).
#'   Defaults to \code{FALSE}.
#'
#' @return An object of class \code{ferx_model} with fields \code{$model} and
#'   \code{$data}. In scaffold mode with \code{print = TRUE}, \code{NULL}
#'   invisibly.
#'
#' @examples
#' ex <- ferx_example("warfarin")
#'
#' # Inspect the object (prints model path, data path, and structure summary)
#' m <- ferx_model(ex$data, ex$model)
#' print(m)
#'
#' # Without data: supply at fit time via ferx_fit(data = ...)
#' m <- ferx_model(model = ex$model)
#'
#' # Scaffold a new model from a template (writes the file, returns a ferx_model)
#' tmp <- tempfile(fileext = ".ferx")
#' m <- ferx_model(template = "1cpt_oral", path = tmp, edit = FALSE)
#' print(m)
#'
#' # Peek at a skeleton without writing anything
#' ferx_model(template = "2cpt_iv", print = TRUE)
#'
#' \dontrun{
#' # ?? Minimal pipe ????????????????????????????????????????????????????????
#' # ferx_fit() picks up $data automatically from the ferx_model object.
#' fit <- ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_fit(method = "focei", covariance = TRUE) |>
#'   summary()
#'
#' # ?? Override fit options before fitting ?????????????????????????????????
#' # ferx_model_set_section() rewrites [fit_options] on disk and passes the
#' # ferx_model through so the pipe continues. When the model file lives
#' # inside the installed package (as for ferx_example()), the file is
#' # copied to tempdir() first so the bundled example is never mutated.
#' fit <- ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_model_set_section("fit_options", c(
#'     "  method     = focei",
#'     "  maxiter    = 500",
#'     "  covariance = true"
#'   )) |>
#'   ferx_fit()
#'
#' summary(fit)
#' ferx_model_inspect(fit)   # structure (no path needed post-fit)
#' fit$cor_matrix            # parameter correlation matrix
#' plot(fit)                 # OFV + gradient norm over iterations
#'
#' # ?? Read a section before fitting ?????????????????????????????????????????
#' lines <- ferx_model_get_section(ex$model, "parameters")
#' fit <- ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_fit(method = "focei")
#'
#' # ?? Validate initialisation before a long run ????????????????????????????
#' # ferx_check_init() runs 5 iterations and returns trace + diagnostics.
#' chk <- ferx_check_init(ex$model, ex$data, method = "focei")
#' chk$summary   # ofv_start, ofv_end, ofv_drop - confirm OFV is dropping
#' plot(chk$fit)             # visual check of first few iterations
#'
#' # ?? Multi-stage chain: SAEM ? FOCEI ?????????????????????????????????????
#' fit <- ex$data |>
#'   ferx_model(ex$model) |>
#'   ferx_fit(method = c("saem", "focei"), covariance = TRUE)
#'
#' # ?? Simulate and predict from fitted parameters ??????????????????????????
#' sim  <- ferx_simulate(ex$model, ex$data, n_sim = 100, seed = 42, fit = fit)
#' pred <- ferx_predict(ex$model, ex$data, fit = fit)
#'
#' # ?? Data can be overridden at fit time ???????????????????????????????????
#' # (substitute the path to your own dataset for "other_cohort.csv")
#' ferx_model(ex$data, ex$model) |>
#'   ferx_fit(data = "other_cohort.csv")
#' }
#'
#' @seealso \code{\link{ferx_model_set_section}}, \code{\link{ferx_model_get_section}},
#'   \code{\link{ferx_fit}}, \code{\link{ferx_check_init}}
#' @family model-editing
#' @export
ferx_model <- function(data = NULL, model = NULL, template = NULL,
                       path = NULL, overwrite = FALSE, edit = TRUE,
                       print = FALSE) {
  # ---- scaffold mode ------------------------------------------------------
  # Create a new .ferx file from a skeleton template, then wrap it in a
  # ferx_model object so the result pipes straight into ferx_fit(). Triggered
  # by `template =` (or `print = TRUE` to peek at a skeleton without writing).
  if (!is.null(template) || isTRUE(print)) {
    skeleton <- .ferx_model_skeleton(if (is.null(template)) "1cpt_oral" else template)
    if (isTRUE(print)) {
      cat(skeleton, sep = "\n")
      cat("\n")
      return(invisible(NULL))
    }
    if (is.null(path)) stop("'path' is required when print = FALSE")
    if (tools::file_ext(path) != "ferx") stop("'path' must end in .ferx")
    if (file.exists(path) && !overwrite) {
      stop(path, " already exists. Use overwrite = TRUE to replace it.")
    }
    # Validate data before any side effect so a bad path doesn't leave a
    # written file / spawned editor behind (matches wrap mode's ordering).
    if (!is.null(data) && !file.exists(data)) stop("Data file not found: ", data)
    writeLines(skeleton, path)
    message("Created ", path)
    if (isTRUE(edit)) utils::file.edit(path)
    return(structure(list(model = path, data = data), class = "ferx_model"))
  }

  # ---- wrap mode ----------------------------------------------------------
  # Backwards-compat shim: in earlier versions the signature was
  # `ferx_model(model, data = NULL)`. Detect old-style positional calls by
  # the .ferx extension on what is now `data` and silently rewrite them
  # with a deprecation warning.
  data_is_ferx <- !is.null(data) && is.character(data) &&
    length(data) == 1L && tolower(tools::file_ext(data)) == "ferx"

  if (data_is_ferx) {
    if (is.null(model)) {
      # `ferx_model("pk.ferx")` ? treat as `ferx_model(model = "pk.ferx")`.
      warning(
        "ferx_model() argument order changed: data is now the first ",
        "argument and model is the second. Call has been auto-corrected ",
        "for compatibility; pass `model =` by name to silence this ",
        "warning. The compatibility shim will be removed in a future release.",
        call. = FALSE
      )
      model <- data
      data  <- NULL
    } else if (is.character(model) && length(model) == 1L &&
               tolower(tools::file_ext(model)) != "ferx") {
      # `ferx_model("pk.ferx", "data.csv")` ? swap arguments.
      warning(
        "ferx_model() argument order changed: data is now the first ",
        "argument and model is the second. Positional call has been ",
        "auto-corrected for compatibility; pass arguments by name to ",
        "silence this warning. The compatibility shim will be removed ",
        "in a future release.",
        call. = FALSE
      )
      tmp   <- model
      model <- data
      data  <- tmp
    }
  }

  if (is.null(model)) {
    stop("'model' is required. Pass a path to a .ferx file.")
  }
  if (!file.exists(model)) stop("File not found: ", model)
  if (tolower(tools::file_ext(model)) != "ferx") stop("'model' must be a .ferx file")
  # When no data path is supplied, fall back to the dataset declared in the
  # model file's `[data]` block (#254) so print() and downstream pipes see it.
  # Resolve before the existence check so a bad declared path fails the same
  # way an explicit `data` argument would, instead of silently constructing a
  # ferx_model that carries a non-existent dataset.
  if (is.null(data)) data <- .ferx_model_data_path(model)
  if (!is.null(data) && !file.exists(data)) stop("Data file not found: ", data)
  structure(list(model = model, data = data), class = "ferx_model")
}

#' @export
print.ferx_model <- function(x, ...) {
  cat("ferx_model\n")
  cat("  Model: ", x$model, "\n", sep = "")
  cat("  Data:  ", if (is.null(x$data)) "<none>" else x$data, "\n", sep = "")
  s <- tryCatch(.ferx_parse_structure(x$model), error = function(e) NULL)
  if (!is.null(s)) {
    cat("  ---\n")
    .ferx_print_structure(s)
  }
  invisible(x)
}

# Return the skeleton lines for a built-in model template.
# Templates: "1cpt_oral" (default), "1cpt_iv", "2cpt_oral", "2cpt_iv", "ode".
# Errors on an unknown template name. Used by ferx_model()'s scaffold mode.
.ferx_model_skeleton <- function(template = "1cpt_oral") {
  templates <- list(
    `1cpt_oral` = c(
      "# One-compartment oral PK model",
      "",
      "[parameters]",
      "  theta TVCL(1.0, 0.001, 100.0)",
      "  theta TVV(10.0, 0.1, 1000.0)",
      "  theta TVKA(1.0, 0.01, 50.0)",
      "",
      "  omega ETA_CL ~ 0.09",
      "  omega ETA_V  ~ 0.09",
      "  omega ETA_KA ~ 0.25",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V  = TVV  * exp(ETA_V)",
      "  KA = TVKA * exp(ETA_KA)",
      "",
      "[structural_model]",
      "  pk one_cpt_oral(cl=CL, v=V, ka=KA)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[fit_options]",
      "  method     = foce",
      "  maxiter    = 300",
      "  covariance = true"
    ),
    `1cpt_iv` = c(
      "# One-compartment IV PK model (RATE column drives bolus vs infusion)",
      "",
      "[parameters]",
      "  theta TVCL(5.0, 0.1, 100.0)",
      "  theta TVV(20.0, 1.0, 500.0)",
      "",
      "  omega ETA_CL ~ 0.09",
      "  omega ETA_V  ~ 0.09",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V  = TVV  * exp(ETA_V)",
      "",
      "[structural_model]",
      "  pk one_cpt_iv(cl=CL, v=V)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[fit_options]",
      "  method     = foce",
      "  maxiter    = 300",
      "  covariance = true"
    ),
    `2cpt_oral` = c(
      "# Two-compartment oral PK model",
      "",
      "[parameters]",
      "  theta TVCL(5.0, 0.1, 100.0)",
      "  theta TVV1(50.0, 1.0, 500.0)",
      "  theta TVQ(10.0, 0.1, 100.0)",
      "  theta TVV2(100.0, 1.0, 1000.0)",
      "  theta TVKA(1.2, 0.01, 10.0)",
      "",
      "  omega ETA_CL ~ 0.10",
      "  omega ETA_V1 ~ 0.10",
      "  omega ETA_Q  ~ 0.05",
      "  omega ETA_V2 ~ 0.05",
      "  omega ETA_KA ~ 0.15",
      "",
      "  sigma PROP_ERR ~ 0.02",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V1 = TVV1 * exp(ETA_V1)",
      "  Q  = TVQ  * exp(ETA_Q)",
      "  V2 = TVV2 * exp(ETA_V2)",
      "  KA = TVKA * exp(ETA_KA)",
      "",
      "[structural_model]",
      "  pk two_cpt_oral(cl=CL, v1=V1, q=Q, v2=V2, ka=KA)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[fit_options]",
      "  method     = focei",
      "  maxiter    = 500",
      "  covariance = true"
    ),
    `2cpt_iv` = c(
      "# Two-compartment IV PK model (RATE column drives bolus vs infusion)",
      "",
      "[parameters]",
      "  theta TVCL(5.0, 0.1, 100.0)",
      "  theta TVV1(15.0, 1.0, 500.0)",
      "  theta TVQ(3.0, 0.01, 100.0)",
      "  theta TVV2(30.0, 1.0, 500.0)",
      "",
      "  omega ETA_CL ~ 0.10",
      "  omega ETA_V1 ~ 0.10",
      "  omega ETA_Q  ~ 0.10",
      "  omega ETA_V2 ~ 0.10",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  CL = TVCL * exp(ETA_CL)",
      "  V1 = TVV1 * exp(ETA_V1)",
      "  Q  = TVQ  * exp(ETA_Q)",
      "  V2 = TVV2 * exp(ETA_V2)",
      "",
      "[structural_model]",
      "  pk two_cpt_iv(cl=CL, v1=V1, q=Q, v2=V2)",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[fit_options]",
      "  method     = foce",
      "  maxiter    = 500",
      "  covariance = true"
    ),
    ode = c(
      "# ODE-based PK model",
      "",
      "[parameters]",
      "  theta TVPARAM(1.0, 0.001, 1000.0)",
      "",
      "  omega ETA_PARAM ~ 0.09",
      "",
      "  sigma PROP_ERR ~ 0.01",
      "",
      "[individual_parameters]",
      "  PARAM = TVPARAM * exp(ETA_PARAM)",
      "",
      "[structural_model]",
      "  ode(obs_cmt=central, states=[depot, central])",
      "",
      "[odes]",
      "  d/dt(depot)   = 0  # replace with your equations",
      "  d/dt(central) = 0",
      "",
      "[error_model]",
      "  DV ~ proportional(PROP_ERR)",
      "",
      "[fit_options]",
      "  method     = focei",
      "  maxiter    = 500",
      "  covariance = true"
    )
  )

  valid <- names(templates)
  if (!template %in% valid) {
    stop(
      "Unknown template '", template, "'. ",
      "Choose one of: ", paste(valid, collapse = ", ")
    )
  }

  templates[[template]]
}
