
loglik.default <- function(model, messages) {
    .Call("Riproc_loglik_new", model, messages)
}

value.loglik <- function(loglik) {
    .Call("Riproc_loglik_value", loglik)
}
