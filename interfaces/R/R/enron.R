
read.enron.nocache <- function() {
    pkg <- "iproc"
    path <- system.file(package = pkg, "extdata", "enron")
    
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
                                           "character",
                                           "character",
                                           "character",
                                           "factor",
                                           "factor"))
    employees$department <- factor(employees$department,
                                   levels = c("Legal", "Trading", "Other"))

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
    n <- nrow(messages)
    receiver.id <- vector("list", n)
    offset <- 1
    last <- nrow(recipients)
    for (i in seq_len(n)) {
        if (offset > last || recipients$message.id[offset] != i) {
            receiver.id[[i]] <- integer(0)
        } else {
            end <- offset + 1
            while(end <= last && recipients$message.id[end] == i)
                end <- end + 1
            receiver.id[[i]] <- recipients$receiver.id[offset:(end-1)]
            offset <- end
        }
    }
    messages$receiver.id <- receiver.id
    messages$nreceiver <- sapply(receiver.id, length)
    
    enron <- list(employees = employees,
                  messages = messages)
    class(enron) <- "enron"

    enron
}

