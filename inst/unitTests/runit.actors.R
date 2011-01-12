
require(RUnit)
require(iproc)

basic <- function(test) {
    setup <- function () {
        v1 <<- c( 1, 100, 2)
        v2 <<- c(-3,   8, 7)
        classes <<- c(1, 1, 2, 1)
        class.traits <<- rbind(v1, v2, deparse.level = 0)
        a <<- actors(classes, class.traits)
    }

    teardown <- function() {
        rm(a, class.traits, classes, v1, v2, inherits = TRUE)
        gc()
    }

    function() { setup(); test(); teardown() }
}

test.dimensions <- basic(function() {
    checkEquals(size.actors(a), length(classes))
    checkEquals(nclass.actors(a), nrow(class.traits))
    checkEquals(dim.actors(a), ncol(class.traits))
})

test.traits <- basic(function() {
    for (i in seq_along(classes)) {
        checkEquals(traits.actors(a, i), class.traits[classes[i],,drop=FALSE])
    }

    ids <- c(3,2)
    checkEquals(traits.actors(a, ids),
                class.traits[classes[ids],,drop=FALSE])
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
        checkEquals(class.traits.actors(a, i), class.traits[i,,drop=FALSE])
    }

    ids <- c(1, 2, 1)
    checkEquals(class.traits.actors(a, ids),
                class.traits[ids,,drop=FALSE])
})
