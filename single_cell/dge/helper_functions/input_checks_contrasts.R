# Helper functions for input checks for contrasts argument in DGE functions.  
#' @param metadata Sample metadata, with rownames containing column names of \code{counts}.
#' @param contrasts NULL or the set of contrasts to extract from a model. must be formatted as a named list.
#' @param model The model formula used for the DGE analysis, as a string. This is used to check that the contrasts are compatible with the model
.input_check_contrasts <- function(metadata, contrasts, model=NULL) {
  if (!is.list(contrasts)) stop("Error. Contrasts must be provided as a named list.")
  if (is.null(names(contrasts)) | any(names(contrasts) == "")) stop("Error. Contrasts must be a named list with non-empty names.")

  model_terms <-  stats::terms(stats::as.formula(model))

  lapply(names(contrasts), function(contrast_name) {
    contrast_str <- contrasts[[contrast_name]]
    if (!is.character(contrast_str) | length(contrast_str) != 1) {
      stop("Error. Contrast '", contrast_name, "' must be a single string.")
    }

    # check that all variables in the contrast are in the model
    vars <- all.vars(parse(text = contrast_str))
    bad <- setdiff(vars, model_terms)
    if (length(bad > 0)) {
      stop("Error. Unknown term(s) in contrast '", contrast_name, "': ",
        paste(bad, collapse = ", ")
      )
    }

    # check that it is parseable
    out <- tryCatch(
      limma::makeContrasts(contrasts = contrast_str, levels = design),
      error = function(e) {
        stop("Error. Contrast '", contrast_name, "' could not be parsed: ", e$message)
      }
    )

    if (is.null(out)) {
      stop("Error. Contrast '", contrast_name, "' could not be parsed. Please check contrast formatting.")
    }
    NULL
  })
}