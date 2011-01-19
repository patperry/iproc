
loglik.default <- function(model, messages) {
    .Call("Riproc_loglik_new", model, messages)
}

value.loglik <- function(loglik) {
    .Call("Riproc_loglik_value", loglik)
}

grad.loglik <- function(loglik) {
    .Call("Riproc_loglik_grad", loglik)
}
