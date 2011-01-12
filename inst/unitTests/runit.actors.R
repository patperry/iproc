
require(RUnit)
require(iproc)

basic <- function(test) {
    setup <- function () {
        v1 <<- c( 1, 100, 2)
        v2 <<- c(-3,   8, 7)
        classes <<- c(1, 1, 2, 1)
        class.vectors <<- rbind(v1, v2, deparse.level = 0)
        a <<- actors(classes, class.vectors)
    }

    teardown <- function() {
        rm(a, class.vectors, classes, v1, v2, inherits = TRUE)
        gc()
    }

    function() { setup(); test(); teardown() }
}

test.dimensions <- basic(function() {
    checkEquals(size.actors(a), length(classes))
    checkEquals(nclass.actors(a), nrow(class.vectors))
    checkEquals(dim.actors(a), ncol(class.vectors))
})

test.vector <- basic(function() {
    for (i in seq_along(classes)) {
        checkEquals(vector.actors(a, i), class.vectors[classes[i],,drop=FALSE])
    }

    ids <- c(3,2)
    checkEquals(vector.actors(a, ids),
                class.vectors[classes[ids],,drop=FALSE])
})

test.class <- basic(function() {
    for (i in seq_along(classes)) {
        checkEquals(class.actors(a, i), classes[i])
    }

    ids <- c(4, 1, 2, 1, 1)
    checkEquals(class.actors(a, ids), classes[ids])
})

test.class.vector <- basic(function() {
    for (i in seq_len(classes)) {
        checkEquals(class.vector.actors(a, i), class.vectors[i,,drop=FALSE])
    }

    ids <- c(1, 2, 1)
    checkEquals(class.vector.actors(a, ids),
                class.vectors[ids,,drop=FALSE])
})
