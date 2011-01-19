
actors.default <- function(groups, group.traits) {
    groups <- as.integer(groups)
    group.traits <- as.matrix(group.traits)
    storage.mode(group.traits) <- 'numeric'
    
    n <- length(groups)
    k <- nrow(group.traits)
    p <- ncol(group.traits)

    if (!(all(1 <= groups & groups <= k)))
        stop("groups label out of bounds")
    
    group.traits.t <- t(group.traits)
    .Call("Riproc_actors_new", groups, group.traits.t)
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

    actors.default(g, gt)
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
