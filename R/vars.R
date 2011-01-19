
vars.default <- function(senders, receivers = senders) {
    .Call("Riproc_vars_new", senders, receivers)
}

dim.vars <- function(vars) {
    .Call("Riproc_vars_dim", vars)
}

nsender.vars <- function(vars) {
    .Call("Riproc_vars_nsender", vars)
}

nreceiver.vars <- function(vars) {
    .Call("Riproc_vars_nreceiver", vars)
}

senders.vars <- function(vars) {
    .Call("Riproc_vars_senders", vars)
}

receivers.vars <- function(vars) {
    .Call("Riproc_vars_receivers", vars)
}

as.matrix.vars <- function(vars, sender, cursor = NULL) {
    p <- dim(vars)
    mul.vars(vars, diag(p), sender, cursor)
}

mul.vars <- function(vars, x, sender, cursor = NULL) {
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    sender <- as.integer(sender)
    .Call("Riproc_vars_mul", vars, x, sender, cursor)
}

tmul.vars <- function(vars, x, sender, cursor = NULL) {
    x <- as.matrix(x)
    storage.mode(x) <- "numeric"
    sender <- as.integer(sender)
    .Call("Riproc_vars_tmul", vars, x, sender, cursor)
}

