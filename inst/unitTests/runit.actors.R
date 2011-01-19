
require(RUnit)
require(iproc)
data(enron)

groups <- group.traits <- a <- NULL

.setUp <- function() {
    e <- enron$employees
    groups <<- (6 * (as.integer(e$gender) - 1)
               + 2 * (as.integer(e$department) - 1)
               + as.integer(e$seniority))
    group.traits <<- diag(12)
    a <<- actors(enron)
    set.seed(0)
}

.tearDown <- function() {
    groups <<- group.traits <<- a <<- NULL
    gc()
}


test.dimensions <- function() {
    checkEquals(size(a), length(groups))
    checkEquals(ngroup(a), nrow(group.traits))
    checkEquals(dim(a), ncol(group.traits))
}

test.traits <- function() {
    for (i in seq_along(groups)) {
        checkEquals(traits(a, i), group.traits[groups[i],,drop=FALSE])
    }

    ids <- c(3,2)
    checkEquals(traits(a, ids),
                group.traits[groups[ids],,drop=FALSE])
}

test.group <- function() {
    for (i in seq_along(groups)) {
        checkEquals(group(a, i), groups[i])
    }

    ids <- c(4, 1, 2, 1, 1)
    checkEquals(group(a, ids), groups[ids])
}

test.group.traits <- function() {
    for (i in seq_len(groups)) {
        checkEquals(group.traits(a, i), group.traits[i,,drop=FALSE])
    }

    ids <- c(1, 2, 1)
    checkEquals(group.traits(a, ids),
                group.traits[ids,,drop=FALSE])
}

test.mul <- function() {
    x <- sample(-2:2, dim(a), replace = TRUE)
    checkEquals(mul(a, x), traits(a) %*% x)
}

test.tmul <- function() {
    x <- sample(-2:2, size(a), replace = TRUE)
    checkEquals(tmul(a, x), t(traits(a)) %*% x)
}
