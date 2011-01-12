
require(RUnit)
require(iproc)

basic <- function(test) {
    setup <- function () {
        within(list(), {
            v1 <- c( 1, 100, 2)
            v2 <- c(-3,   8, 7)
            groups <- c(1, 1, 2, 1)
            group.traits <- rbind(v1, v2, deparse.level = 0)
            a <- actors(groups, group.traits)
        })
    }

    teardown <- function() {
        gc()
    }

    function() {
        env <- setup()
        attach(env)
        test()
        teardown()
        detach(env)
    }
}

test.dimensions <- basic(function() {
    checkEquals(size(a), length(groups))
    checkEquals(ngroup(a), nrow(group.traits))
    checkEquals(dim(a), ncol(group.traits))
})

test.traits <- basic(function() {
    for (i in seq_along(groups)) {
        checkEquals(traits(a, i), group.traits[groups[i],,drop=FALSE])
    }

    ids <- c(3,2)
    checkEquals(traits(a, ids),
                group.traits[groups[ids],,drop=FALSE])
})

test.group <- basic(function() {
    for (i in seq_along(groups)) {
        checkEquals(group(a, i), groups[i])
    }

    ids <- c(4, 1, 2, 1, 1)
    checkEquals(group(a, ids), groups[ids])
})

test.group.traits <- basic(function() {
    for (i in seq_len(groups)) {
        checkEquals(group.traits(a, i), group.traits[i,,drop=FALSE])
    }

    ids <- c(1, 2, 1)
    checkEquals(group.traits(a, ids),
                group.traits[ids,,drop=FALSE])
})
