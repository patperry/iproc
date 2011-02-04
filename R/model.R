
model.default <- function(frame, coef, has.loops = FALSE) {
    coef <- as.numeric(coef)
    has.loops <- as.logical(has.loops)

    if (length(has.loops) != 1)
        stop("'has.loops' should have length 1")
    
    model <- .Call("Riproc_model_new", frame, coef, has.loops);
    attr(model, "attributes.frame") <- attributes(frame)
    model
}

iproc.frame.model <- function(object, ...) {
    model <- object
    frame <- .Call("Riproc_model_frame", model)
    attributes(frame) <- attr(model, "attributes.frame")
    frame
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

probs.model <- function(model, sender, cursor = NULL) {
    lp <- log.probs.model(model, sender, cursor)
    exp(lp)
}

log.probs.model <- function(model, sender, cursor = NULL) {
    sender <- as.integer(sender)
    lpt <- .Call("Riproc_model_log_probs", model, sender, cursor)
    t(lpt)
}
