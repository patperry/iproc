
loglik.model <- function(model, messages) {
    .Call("Riproc_model_loglik", model, messages)
}
