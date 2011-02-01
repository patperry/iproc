
cursor.messages <- function(messages) {
    .Call("Riproc_cursor_new", messages)
}

reset.cursor <- function(cursor) {
    .Call("Riproc_cursor_reset", cursor)
}

advance.cursor <- function(cursor) {
    .Call("Riproc_cursor_advance", cursor)
}

started.cursor <- function(cursor) {
    .Call("Riproc_cursor_started", cursor)
}

finished.cursor <- function(cursor) {
    .Call("Riproc_cursor_finished", cursor)
}

time.cursor <- function(cursor) {
    .Call("Riproc_cursor_time", cursor)
}

nties.cursor <- function(cursor) {
    .Call("Riproc_cursor_nties", cursor)
}

from.cursor <- function(cursor) {
    .Call("Riproc_cursor_from", cursor)
}

to.cursor <- function(cursor) {
    .Call("Riproc_cursor_to", cursor)
}

print.cursor <- function(x, ...) {
    if (!started(x)) {
        cat("(unstarted cursor)\n")
    } else if (finished(x)) {
        cat("(finished cursor)\n")
    } else {
        cat("cursor at time ", time(x), ":\n", sep='')
        from <- from(x)
        to <- to(x)

        for (i in seq_along(from)) {
            cat("    ", from[[i]], " -> {", sep='')
            nto <- length(to[[i]])
            if (nto > 0) {
                cat(" ", to[[i]][1], sep='')
            }
            for (j in (1 + seq_len(nto - 1))) {
                cat(", ", to[[i]][j], sep='')
            }
            cat(" }\n")
        }
    }
    invisible(x)
}
