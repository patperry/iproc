

iproc.fit <- function(model0, messages, penalty) {
    penalty <- as.numeric(penalty)
    model <- .Call("Riproc_fit", model0, messages, penalty)
    model
}
