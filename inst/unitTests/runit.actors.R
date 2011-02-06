
require(RUnit)
require(iproc)
data(enron)

a <- traits <- ngroup <- NULL

.setUp <- function() {
    e <- enron$employees
    mm <- model.matrix(~ gender*seniority*department, enron$employees)
    attr(mm, "dimnames") <- NULL
    attr(mm, "assign") <- NULL
    attr(mm, "contrasts") <- NULL

    traits <<- mm
    ngroup <<- 12
    a <<- actors(~ gender *seniority * department, enron$employees)
    set.seed(0)
}

.tearDown <- function() {
    a <<- NULL
    gc()
}

test.dimensions <- function() {
    checkEquals(nrow(a), nrow(traits))
    checkEquals(ncol(a), ncol(traits))
    checkEquals(ngroup(a), ngroup)
}

test.group.traits <- function() {
    g <- group(a)
    gt <- group.traits(a)
    checkEquals(gt[g,,drop = FALSE], traits)
}

test.as.matrix <- function() {
    checkEquals(as.matrix(a), traits)
}

test.mul <- function() {
    x <- sample(-2:2, ncol(a), replace = TRUE)
    checkEquals(mul(a, x), as.matrix(a) %*% x)
}

test.tmul <- function() {
    x <- sample(-2:2, nrow(a), replace = TRUE)
    checkEquals(tmul(a, x), t(as.matrix(a)) %*% x)
}
