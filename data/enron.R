
enron <- local({
    pkg <- "iproc"
    if (Sys.getenv("RCMDCHECK") == "FALSE") {
        ## PKG/tests/../inst/extdata/enron
        path <- file.path(getwd(), "..", "inst", "extdata", "enron")
    } else {
        path <- system.file(package = pkg, "extdata", "enron")
    }

    employees <- read.table(file.path(path, "employees.tsv"),
                            sep = "\t",
                            quote = "",
                            col.names = c("id",
                                          "name",
                                          "department",
                                          "long.department",
                                          "title",
                                          "gender",
                                          "seniority"),
                            colClasses = c("integer",
                                           "character",
                                           "factor",
                                           "character",
                                           "character",
                                           "factor",
                                           "factor"))

    messages <- read.table(file.path(path, "messages.tsv"),
                           sep = "\t",
                           quote = "",
                           col.names = c("id",
                                         "filename",
                                         "unix.time",
                                         "subject",
                                         "sender.id"),
                           colClasses = c("integer",
                                          "character",
                                          "integer",
                                          "character",
                                          "integer"))
    time <- as.POSIXct(messages$unix.time, origin = "1970-01-01", tz = "UTC")
    messages$time <- time
    messages <- messages[,c(1,2,6,3,4,5)]

    recipients <- read.table(file.path(path, "recipients.tsv"),
                             sep = "\t",
                             quote = "",
                             col.names = c("message.id",
                                           "position",
                                           "receiver.id"),
                             colClasses = c("integer",
                                            "integer",
                                            "integer"))                            

    enron <- list(employees = employees,
                  messages = messages,
                  recipients = recipients)
    class(enron) <- "enron"
    enron
})
