
actors <- function (classes, class.vectors) {
    classes <- as.integer(classes)
    class.vectors <- as.matrix(class.vectors)
    storage.mode(class.vectors) <- 'numeric'
    
    n <- length(classes)
    k <- nrow(class.vectors)
    p <- ncol(class.vectors)

    if (!(all(1 <= classes & classes <= k)))
        stop("classs label out of bounds")
    
    class.vectors.t <- t(class.vectors)
    .Call("Riproc_actors_new", classes, class.vectors.t)
}

nclass.actors <- function(actors) {
    .Call("Riproc_actors_nclass", actors)
}

size.actors <- function(actors) {
    .Call("Riproc_actors_size", actors)
}

dim.actors <- function(actors) {
    .Call("Riproc_actors_dim", actors)
}

vector.actors <- function(actors, ids) {
    ids <- as.integer(ids)
    xt <- .Call("Riproc_actors_vector", actors, ids)
    t(xt)
}

class.actors <- function(actors, ids) {
    ids <- as.integer(ids)
    .Call("Riproc_actors_class", actors, ids)
}

class.vector.actors <- function(actors, class.ids) {
    class.ids <- as.integer(class.ids)
    xt <- .Call("Riproc_actors_class_vector", actors, class.ids)
    t(xt)
}

vars <- function(senders, receivers) {

}

model <- function(vars, coefs, has.loops=FALSE) {

}
