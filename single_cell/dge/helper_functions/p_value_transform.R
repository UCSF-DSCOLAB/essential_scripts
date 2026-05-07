#' Produces rounded p-values or *-transformed p-values
#' @param pvals Numeric vector of adjusted or unadjusted p-values
#' @parma summarize Logical, TRUE by default, naming whether p-values should be summarized into bins of not significant, significant, more significant, ect.
#' with the exact details of these bins and the exact summary style determined by additional inputs below.
#' @param round Logical, TRUE by default, naming whether pvalues should be rounded to \code{round_digits} digits.
#' @param round_digits Integer giving the number of decimal points for rounding
#' @param sig_cut_l1,sig_cut_l2,sig_cut_l3 Numbers which set the boundaries of significance bins when \code{summarize = TRUE}.
#' In particular, \code{sig_cut_l1} set the primary cutoff where \code{pvals > sig_cut_l1} are deemed not significant.
#' @param non_sig_label String, "NS" by default, giving the label used for non-significant pvalues when \code{summarize = TRUE}
#' @param summarize_style One of '*', 'value', or 'value-l1' which toggles three distinct options for how pvals are transformed when \code{summarize = TRUE}:
#' \itemize{
#' \item \code{*} = 4 bins using all \code{sig_cut} inputs labeled as \code{non_sig_label}, \code{*}, \code{**}, and \code{***}
#' \item \code{value} = 2 bins using only \code{sig_cut_l1}, labeled as \code{non_sig_label} for non-significant pvals, and the pval itself otherwise (rounded if \code{round = TRUE})
#' \item \code{value-l1} = 3 bins set using \code{sig_cut_l1} and \code{sig_cut_l2}, labeled as \code{non_sig_label} for non-significant pvals, the pval itself for pvals between \code{sig_cut_l1} and \code{sig_cut_l2} (rounded if \code{round = TRUE}), and '\code{value_style_p_string} < \code{sig_cut_l2}' for pvals smaller than \code{sig_cut_l2}.
#' }
#' @param value_style_p_string String, 'padj' by default, giving a name to call the p-values when \code{summarize = TRUE} and \vode{summarize_style} is \code{'value'} or \code{'value-l1'}.
#' Note that when p-values are not multiple hypothesis corrected this should be changed to 'p' or similar to not give the wrong impression.
#' @param value_style_symbol String, '=' by default, giving the spacer to use between \code{value_style_p_string} and the p-value descriptor, when \code{summarize = TRUE} and \vode{summarize_style} is \code{'value'} or \code{'value-l1'}.
#' @param value_style_compact Logical determining whether to include spaces between \code{value_style_p_string}, \code{value_style_symbol}, and the p-value descriptor, when \code{summarize = TRUE} and \vode{summarize_style} is \code{'value'} or \code{'value-l1'}.
#' @param warn_if_nothing_done Logical which can be set to FALSE to remove the warning provided if the function is called with parameters that turn off all functions, a.k.a. \code{sumarrize = TRUE} and \code{round = TRUE}.
#' @details This function rounds and/or transforms given \code{pvals} into any of three summary formats.
#'
#' The function DOES NOT perform multiple hypothesis adjustments.  Such adjustment should be performed before running this function.
#'
#' Rounding is performed is \code{round} is \code{TRUE}, to \code{round_digits} number of digits.
#'
#' Pvalues are summarized when \code{summarize} is \code{TRUE}, by first binning pvalues,
#' then labeling in one of three general ways depending on the value of \code{summarize_style}:
#' \itemize{
#' \item \code{*} = 4 bins: 'NS', \code{*}, \code{**}, \code{***}
#' \item \code{value} = 2 bins: \code{NS}, \code{<pval>}
#' \item \code{value-l1} = 4 bins: \code{NS}, \code{<pval>}, \code{pval < 0.001}, \code{pval < 0.0001}
#' }
#'
#' Bin delineations are given by \code{sig_cut_l1}, \code{sig_cut_l2}, and \code{sig_cut_l3}:
#' \itemize{
#' \item Non-significant: \code{pvals > sig_cut_l1} (Default: \code{pvals > 0.05})
#' \item Level1 significance: \code{sig_cut_l1 >= pvals >= sig_cut_l2} (Default: \code{0.05 > pvals >= 0.001})
#' \item Level2 significance: \code{sig_cut_l2 > pvals >= sig_cut_l3} (Default: \code{0.001 > pvals >= 0.0001})
#' \item Level3 significance: \code{sig_cut_l3 > pvals} (Default: \code{pvals < 0.0001})
#' }
#'
#' Note that exact styles of these labels can be adjusted via \code{non_sig_label}, \code{value_style_p_string}, \code{value_style_symbol}, and \code{value_style_compact}.
#'
#' @examples
#' ex_pvals <- c(0.5, 0.055, 0.03, 0.003, 0.0003, 0.00003)
#'
#' p_value_transform(ex_pvals)
#'
#' # To only round pvalues, set 'summarize = FALSE'
#' p_value_transform(ex_pvals, summarize = FALSE)
#' # Note that if 'round_digits' is set too small, you can get labels of '0'
#' p_value_transform(ex_pvals, summarize = FALSE, round_digits = 3)
#'
#' # Change 'summarize_style' to 'value' or 'value-l1' to have a different
#' # summarization scheme.
#' p_value_transform(ex_pvals, summarize_style = 'value')
#' p_value_transform(ex_pvals, summarize_style = 'value-l1')
#'
#' # If your pvalues are not multiple hypothesis adjusted, it is best to change
#' # 'value_style_p_string' to 'p':
#' p_value_transform(ex_pvals, summarize_style = 'value',
#'     value_style_p_string = 'p')
#'
#' # To entirely remove the 'pval = ' portion from labels, use all of
#' # 'value_style_p_string', 'value_style_symbol' and 'value_style_compact':
#' p_value_transform(ex_pvals, summarize_style = 'value',
#'     value_style_p_string = '',
#'     value_style_symbol = '',
#'     value_style_compact = TRUE)
#' p_value_transform(ex_pvals, summarize_style = 'value-l1',
#'     value_style_p_string = '',
#'     value_style_symbol = '',
#'     value_style_compact = TRUE)
#'
#' # Change 'sig_cut_l1', 'sig_cut_l2', or 'sig_cut_l3' to adjust significance and
#' # binning cutoffs. E.g. the below example will highlight and show all p-values
#' # within the murky range of 0.05 to 0.1:
#' p_value_transform(ex_pvals,
#'     summarize_style = 'value-l1',
#'     sig_cut_l1 = 0.1,
#'     sig_cut_l2 = 0.05)
#' p_value_transform(ex_pvals,
#'     summarize_style = 'value-l1',
#'     sig_cut_l1 = 0.1,
#'     sig_cut_l2 = 0.05,
#'     sig_cut_l3 = 0.01,
#'     value_style_p_string = '',
#'     value_style_symbol = '',
#'     value_style_compact = TRUE)
#' @author Dan Bunis
p_value_transform <- function(
    pvals,
    summarize = TRUE,
    round = TRUE,
    round_digits = 5,
    sig_cut_l1 = 0.05,
    sig_cut_l2 = 0.001,
    sig_cut_l3 = 0.0001,
    non_sig_label = "NS",
    summarize_style = c('*', 'value', 'value-l1'),
    value_style_p_string = 'padj',
    value_style_symbol = '=',
    value_style_compact = FALSE,
    warn_if_nothing_done = TRUE
) {

    if (!summarize && !round && warn_if_nothing_done) {
        warning("Nothing to do. 'summarize' and 'round' are both FALSE.")
        return (pvals)
    }
    summarize_style <- match.arg(summarize_style)
    if (sig_cut_l2 >= sig_cut_l1) {
        stop("'sig_cut_l2' should be smaller than 'sig_cut_l1'")
    }

    if (sig_cut_l3 >= sig_cut_l2) {
        stop("'sig_cut_l3' should be smaller than 'sig_cut_l2'")
    }
    if (sig_cut_l3 >= sig_cut_l1) {
        stop("'sig_cut_l3' should be smaller than 'sig_cut_l1'")
    }

    rounded <- round(pvals, round_digits)
    p_use <- out <- if (round) { rounded } else { pvals }

    if (summarize) {
        level <- ifelse(
            pvals > sig_cut_l1,
            0,
            ifelse(
                pvals >= sig_cut_l2,
                1,
                ifelse(
                    pvals >= sig_cut_l3,
                    2,
                    3
                )
            )
        )
        sep <- ifelse(value_style_compact, "", " ")
        if (summarize_style=='value') {
            out[level==0] <- non_sig_label
            out[level>0] <- paste(
                sep = sep, value_style_p_string, value_style_symbol,
                format(p_use[level>0], drop0trailing = TRUE, scientific=FALSE))
        } else if (summarize_style=='value-l1') {
            out[level==0] <- non_sig_label
            out[level==1] <- paste(
                sep = sep, value_style_p_string, value_style_symbol,
                format(p_use[level==1], drop0trailing = TRUE, scientific=FALSE))
            out[level==2] <- paste(
                sep = sep, value_style_p_string, "<",
                format(sig_cut_l2, drop0trailing = TRUE, scientific=FALSE))
            out[level==3] <- paste(
                sep = sep, value_style_p_string, "<",
                format(sig_cut_l3, drop0trailing = TRUE, scientific=FALSE))
        } else {
            out[level==0] <- non_sig_label
            out[level==1] <- "*"
            out[level==2] <- "**"
            out[level==3] <- "***"
        }
    } else {
        out <- format(out, drop0trailing = TRUE, scientific=FALSE)
    }

    out
}
