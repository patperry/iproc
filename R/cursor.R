
cursor.messages <- function(messages) {
    .Call("Riproc_cursor_new", messages)
}

reset.cursor <- function(cursor) {
    .Call("Riproc_cursor_reset", cursor)
}

advance.cursor <- function(cursor) {
    .Call("Riproc_cursor_advance", cursor)
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
