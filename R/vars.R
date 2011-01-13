
vars <- function(senders, receivers = senders) {
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
