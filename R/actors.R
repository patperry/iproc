
actors <- function (classes, class.traits) {
    classes <- as.integer(classes)
    class.traits <- as.matrix(class.traits)
    storage.mode(class.traits) <- 'numeric'
    
    n <- length(classes)
    k <- nrow(class.traits)
    p <- ncol(class.traits)

    if (!(all(1 <= classes & classes <= k)))
        stop("classs label out of bounds")
    
    class.traits.t <- t(class.traits)
    .Call("Riproc_actors_new", classes, class.traits.t)
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

traits.actors <- function(actors, ids) {
    ids <- as.integer(ids)
    xt <- .Call("Riproc_actors_traits", actors, ids)
    t(xt)
}

class.actors <- function(actors, ids) {
    ids <- as.integer(ids)
    .Call("Riproc_actors_class", actors, ids)
}

class.traits.actors <- function(actors, class.ids) {
    class.ids <- as.integer(class.ids)
    xt <- .Call("Riproc_actors_class_traits", actors, class.ids)
    t(xt)
}
