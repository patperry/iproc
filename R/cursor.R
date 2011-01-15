
cursor.messages <- function(messages) {
    .Call("Riproc_cursor_new", messages)
}

reset.cursor <- function(cursor) {
    .Call("Riproc_cursor_reset", cursor)
}

advance.cursor <- function(cursor) {
    .Call("Riproc_cursor_advance", cursor)
}
