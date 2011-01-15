
model.default <- function(vars, coef, has.loops = FALSE) {
    coef <- as.numeric(coef)
    has.loops <- as.logical(has.loops)

    if (length(has.loops) != 1)
        stop("'has.loops' should have length 1")
    
    .Call("Riproc_model_new", vars, coef, has.loops);
}

vars.model <- function(model) {
    .Call("Riproc_model_vars", model)
}

coef.model <- function(model) {
    .Call("Riproc_model_coefs", model)
}

has.loops.model <- function(model) {
    .Call("Riproc_model_has_loops", model)
}

dim.model <- function(model) {
    .Call("Riproc_model_dim", model)
}

nsender.model <- function(model) {
    .Call("Riproc_model_nsender", model)
}

nreceiver.model <- function(model) {
    .Call("Riproc_model_nreceiver", model)
}

probs.model <- function(model, isend, cursor = NULL) {
    lp <- log.probs.model(model, isend, cursor)
    exp(lp)
}

log.probs.model <- function(model, isend, cursor = NULL) {
    isend <- as.integer(isend)
    lpt <- .Call("Riproc_model_log_probs", model, isend, cursor)
    t(lpt)
}
