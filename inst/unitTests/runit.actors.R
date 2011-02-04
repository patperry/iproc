
require(RUnit)
require(iproc)
data(enron)

a <- traits <- ngroup <- NULL

.setUp <- function() {
    e <- enron$employees
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
    traits <<- gt[g,,drop=FALSE]
    ngroup <<- 12
    a <<- actors(enron)
    set.seed(0)
}

.tearDown <- function() {
    a <<- NULL
    gc()
}


test.dimensions <- function() {
    checkEquals(size(a), nrow(traits))
    checkEquals(ngroup(a), ngroup)
    checkEquals(dim(a), ncol(traits))
}

test.traits <- function() {
    for (i in seq_len(size(a))) {
        checkEquals(traits(a, i), traits[i,,drop=FALSE])
    }

    checkEquals(traits(a), traits)
    
    ids <- c(3,2)
    checkEquals(traits(a, ids), traits[ids,,drop=FALSE])
}

test.group.traits <- function() {
    for (i in seq_len(size(a))) {
        g <- group(a, i)
        checkEquals(group.traits(a, g), traits(a, i))
    }

    ids <- 5:10
    gs <- group(a, ids)
    checkEquals(group.traits(a, gs), traits(a, ids))
}

test.as.matrix <- function() {
    checkEquals(as.matrix(a), traits(a))
}

test.mul <- function() {
    x <- sample(-2:2, dim(a), replace = TRUE)
    checkEquals(mul(a, x), as.matrix(a) %*% x)
}

test.tmul <- function() {
    x <- sample(-2:2, size(a), replace = TRUE)
    checkEquals(tmul(a, x), t(as.matrix(a)) %*% x)
}
