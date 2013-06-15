
hash.numeric <- function(x) {
    x <- as.numeric(x)
    .Call("Riproc_hash_numeric", x)
}
