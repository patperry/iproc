
actors.default <- function(traits) {
    traits <- as.matrix(traits)
    traits.t <- t(traits)

    .Call("Riproc_actors_new", traits.t)
}

actors.enron <- function(enron) {
    group <- function(emp) {
        (6 * (as.integer(emp$gender) - 1)
         + 2 * (as.integer(emp$department) - 1)
         + as.integer(emp$seniority))
    }

    g <- rep(NA, nrow(enron$employees))
    for (i in seq_along(g)) {
        g[i] <- group(enron$employees[i,])
    }

    gt <- diag(12)

    traits <- gt[g,,drop=FALSE]

    actors.default(traits)
}

ngroup.actors <- function(actors) {
    .Call("Riproc_actors_ngroup", actors)
}

size.actors <- function(actors) {
    .Call("Riproc_actors_size", actors)
}

dim.actors <- function(actors) {
    .Call("Riproc_actors_dim", actors)
}

traits.actors <- function(actors, ids = NULL) {
    if (is.null(ids)) {
        ids <- seq_len(size(actors))
    } else {
        ids <- as.integer(ids)
    }
    xt <- .Call("Riproc_actors_traits", actors, ids)
    t(xt)
}

group.actors <- function(actors, ids = NULL) {
    if (is.null(ids)) {
        ids <- seq_len(size(actors))
    } else {
        ids <- as.integer(ids)
    }
    .Call("Riproc_actors_group", actors, ids)
}

group.traits.actors <- function(actors, group.ids = NULL) {
    if (is.null(group.ids)) {
        group.ids <- seq_len(ngroup(actors))
    } else {
        group.ids <- as.integer(group.ids)
    }
    
    xt <- .Call("Riproc_actors_group_traits", actors, group.ids)
    t(xt)
}

as.matrix.actors <- function(actors) {
    traits(actors)
}
    
mul.actors <- function(actors, matrix) {
    matrix <- as.matrix(matrix)
    storage.mode(matrix) <- "numeric"
    .Call("Riproc_actors_mul", actors, matrix)
}

tmul.actors <- function(actors, matrix) {
    matrix <- as.matrix(matrix)
    storage.mode(matrix) <- "numeric"
    .Call("Riproc_actors_tmul", actors, matrix)
}
