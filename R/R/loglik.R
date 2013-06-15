
loglik.default <- function(model, messages = NULL) {
    .Call("Riproc_loglik_new", model, messages)
}

insert.loglik <- function(loglik, cursor) {
    .Call("Riproc_loglik_insert", loglik, cursor);
}

value.loglik <- function(loglik) {
    .Call("Riproc_loglik_value", loglik)
}

grad.loglik <- function(loglik) {
    .Call("Riproc_loglik_grad", loglik)
}
