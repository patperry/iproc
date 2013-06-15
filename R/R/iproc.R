
iproc <- function(messages, senders, receivers = senders, recip.intervals = NULL,
                  has.loops = FALSE, data.senders, data.receivers = data.senders,
                  penalty = 0.0, control = list(...), ...)
{
    cl <- match.call()
    design <- match.call(expand.dots = FALSE)
    dargs <- match(c("senders", "receivers", "recip.intervals",
                     "data.senders", "data.receivers"), names(design), 0L)
    design <- design[c(1L, dargs)]
    design[[1L]] <- as.name("iproc.design")
    design <- eval(design, parent.frame())

    control <- do.call("iproc.control", control)
    
    beta0 <- rep(0, ncol(design))
    model0 <-  model(design, beta0, has.loops)
    model <- iproc.fit(model0, messages, penalty, control)
    browser()
}
