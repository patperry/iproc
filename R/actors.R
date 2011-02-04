
# object should be a formula or a terms object
# no support for subset, since it makes ids confusing
actors.default <- function(object, data, na.action, contrasts = NULL, ...)
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("object", "data", "na.action"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    names(mf)[[2]] <- "formula"

    mf <- eval(mf, parent.frame())

    #mf <- lm(object, quote(data), quote(na.action), method = "model.frame")
    mt <- attr(mf, "terms")
    if (attr(mt, "response") > 0)
        stop("response should be left unspecified, but was '",
             names(mf)[[attr(mt, "response")]], "'")
    if (nrow(mf) == 0)
        stop("number of actors is 0, should be positive")


    x <- model.matrix(mt, mf, contrasts)
    actrs <- .Call("Riproc_actors_new", t(x))
    actrs
}

actors.enron <- function(object, ...) {
    actors( ~ gender*seniority*department, object$employees, ...)
}

ngroup.actors <- function(object, ...) {
    .Call("Riproc_actors_ngroup", object)
}

dim.actors <- function(object, ...) {
    n <- .Call("Riproc_actors_size", object)
    p <- .Call("Riproc_actors_dim", object)
    c(n,p)
}

group.actors <- function(object, ...) {
    ids <- seq_len(nrow(object))
    .Call("Riproc_actors_group", object, ids)
}

group.traits.actors <- function(object, ...) {
    group.ids <- seq_len(ngroup(object))
    xt <- .Call("Riproc_actors_group_traits", object, group.ids)
    t(xt)
}

as.matrix.actors <- function(x, ...) {
    ids <- seq_len(nrow(x))
    mat <- .Call("Riproc_actors_traits", x, ids)
    t(mat)
}
    
mul.actors <- function(x, y, ...) {
    actors <- x
    matrix <- as.matrix(y)
    storage.mode(matrix) <- "numeric"
    .Call("Riproc_actors_mul", actors, matrix)
}

tmul.actors <- function(x, y, ...) {
    actors <- x
    matrix <- as.matrix(y)
    storage.mode(matrix) <- "numeric"
    .Call("Riproc_actors_tmul", actors, matrix)
}
