

test.basic <- function() {
    a <- actors(c(1,1,2,1), matrix(1:6, 2, 3))

    checkEquals(size.actors(a), 4)
    checkEquals(nclass.actors(a), 2)
    checkEquals(dim.actors(a), 3)
}
