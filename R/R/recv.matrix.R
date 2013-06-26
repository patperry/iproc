
recv.matrix.variable <- function(var, contrasts = NULL)
{
    if (inherits(var, "data.frame") && !inherits(var, "response")) {
        mm <- model.matrix(var, contrasts.arg=contrasts)
        ret <- mm

        int <- match("(Intercept)", colnames(ret))
        if (!is.na(int)) {
            contrasts <- attr(ret, "contrasts")
            ret <- ret[,-int,drop=FALSE]
            attr(ret, "contrasts") <- contrasts[[1L]]
        }

        if (inherits(var, "send"))
            colnames(ret) <- paste("s(", colnames(ret), ")", sep="")

        cl <- class(var)
        cl[match("data.frame", cl)] <- "matrix"
        class(ret) <- cl
    } else {
        ret <- var
    }

    ret
}


recv.matrix <- function(frame, contrasts = NULL)
{
    variables <- list()
    for (i in seq_along(frame$variables)) {
        var <- frame$variables[[i]]
        if (i == frame$response) {
            variables[[i]] <- var
        } else {
            variables[[i]] <- recv.matrix.variable(var, contrasts)
        }
    }
    names(variables) <- names(frame$variables)

    mat <- frame
    mat$variables <- variables

    response <- frame$response

    if (response != 0L) {
        mat$variables <- mat$variables[-response]

	if (!identical(mat$factors, integer(0)))
		mat$factors <- mat$factors[-response,,drop=FALSE]
    }

    class(mat) <- "recv.matrix"
    mat
}


model.matrix.recv.frame <- function(object, time, sender, ...)
{
    stop("not implemented")
}

