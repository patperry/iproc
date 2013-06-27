
summary.recv.model <- function (object, correlation = FALSE, symbolic.cor = FALSE, ...) 
{
    dispersion <- 1
    est.disp <- FALSE
    df.r <- object$df.residual
    coef <- object$coefficients
    p <- length(coef)
    if (p > 0) {
        p1 <- 1L:p
        covmat.unscaled <- vcov(object)
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef/s.err
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef), c(dn, 
                "z value", "Pr(>|z|)"))
        }
        else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef), c(dn, 
                "t value", "Pr(>|t|)"))
        }
        else {
            coef.table <- cbind(coef, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef), c(dn, 
                "t value", "Pr(>|t|)"))
        }
        df.f <- object$rank
    }
    else {
        coef.table <- matrix(, 0L, 4L)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error", 
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
        df.f <- object$rank
    }
    keep <- match(c("call", "deviance", "aic", 
        "df.residual", "null.deviance", "df.null", 
        "iter"), names(object), 0L)
    ans <- c(object[keep], list(
             coefficients = coef.table,
        dispersion = dispersion, df = c(object$rank, df.r, df.f), 
        cov.unscaled = covmat.unscaled, cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.recv.model"
    return(ans)
}


print.summary.recv.model <-
    function (x, digits = max(3, getOption("digits") - 3),
	      symbolic.cor = x$symbolic.cor,
	      signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n",
	paste(deparse(x$call), sep="\n", collapse = "\n"), "\n", sep="")

    if(nrow(x$coefficients) == 0L) {
        cat("\nNo Coefficients\n")
    } else {
        ## df component added in 1.8.0
        ## partial matching problem here.
        df <- if ("df" %in% names(x)) x[["df"]] else NULL
        if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
            cat("\nCoefficients: (", nsingular,
                " not defined because of singularities)\n", sep = "")
        else cat("\nCoefficients:\n")
        coefs <- x$coefficients
        if(!is.null(aliased <- x$aliased) && any(aliased)) {
            cn <- names(aliased)
            coefs <- matrix(NA, length(aliased), 4L,
                            dimnames=list(cn, colnames(coefs)))
            coefs[!aliased, ] <- x$coefficients
        }
        printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                     na.print="NA", ...)
    }
    ##
    cat("\n(Dispersion parameter",
	" taken to be ", format(x$dispersion), ")\n\n",
	apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                          "deviance:"),
		    format(unlist(x[c("null.deviance","deviance")]),
			   digits= max(5, digits+1)), " on",
		    format(unlist(x[c("df.null","df.residual")])),
		    " degrees of freedom\n"),
	      1L, paste, collapse=" "), sep="")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (",mess, ")\n", sep="")
    cat("AIC: ", format(x$aic, digits= max(4, digits+1)),"\n\n",
	"Number of Newton iterations: ", x$iter,
	"\n", sep="")

    correl <- x$correlation
    if(!is.null(correl)) {
# looks most sensible not to give NAs for undefined coefficients
#         if(!is.null(aliased) && any(aliased)) {
#             nc <- length(aliased)
#             correl <- matrix(NA, nc, nc, dimnames = list(cn, cn))
#             correl[!aliased, !aliased] <- x$correl
#         }
	p <- NCOL(correl)
	if(p > 1) {
	    cat("\nCorrelation of Coefficients:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
		print(symnum(correl, abbr.colnames = NULL))
	    } else {
		correl <- format(round(correl, 2), nsmall = 2, digits = digits)
		correl[!lower.tri(correl)] <- ""
		print(correl[-1, -p, drop=FALSE], quote = FALSE)
	    }
	}
    }
    cat("\n")
    invisible(x)
}



