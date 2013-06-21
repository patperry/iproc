get.factors <- function(terms)
{
    v <- attr(terms, "variables")[-1L]

    factors <- attr(terms, "factors")
    if (length(factors) == 0L) {
        factors <- matrix(0L, attr(terms, "response"), 0L)
        colnames(factors) <- attr(terms, "term.labels")
        if (attr(terms, "response") > 0L) {
            rownames(factors) <- deparse(v[[1L]])
        } else {
            rownames(factors) <- character(0L)
        }
    }

    factors
}
