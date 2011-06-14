
model.default <- function(design, coef) {
    coef <- as.numeric(coef)

    model <- .Call("Riproc_model_new", design, coef)
    attr(model, "attributes.design") <- attributes(design)
    model
}

iproc.design.model <- function(object, ...) {
    model <- object
    design <- .Call("Riproc_model_design", model)
    attributes(design) <- attr(model, "attributes.design")
    design
}

coef.model <- function(model) {
    .Call("Riproc_model_coefs", model)
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

probs.model <- function(model, sender, cursor = NULL) {
    lp <- log.probs.model(model, sender, cursor)
    exp(lp)
}

log.probs.model <- function(model, sender, cursor = NULL) {
    sender <- as.integer(sender)
    lpt <- .Call("Riproc_model_log_probs", model, sender, cursor)
    t(lpt)
}
