
require(RUnit)
require(iproc)

data(enron)

test.employees.nrow <- function() {
    checkEquals(nrow(enron$employees), 156)
}

test.employees.data <- function() {
    e <- enron$employees[c(1,34,156),]
    name <- c("John Arnold", "Drew Fossum", "Phillip K. Allen")
    department <- factor(c("Trading", "Legal", "Trading"),
                         levels = c("Legal", "Trading", "Other"))
    long.department <- c("ENA Gas Financial", "ETS", "ENA Gas West")
    title <- c("VP Trading", "VP & Gen Cnsl", "Mng Dir Trading")
    gender <- factor(c("Male", "Male", "Male"),
                     levels = c("Female", "Male"))
    seniority <- factor(c("Senior", "Senior", "Senior"),
                        levels = c("Junior", "Senior"))
    
    checkEquals(e$name, name)
    checkEquals(e$department, department)
    checkEquals(e$long.department, long.department)
    checkEquals(e$title, title)
    checkEquals(e$gender, gender)
    checkEquals(e$seniority, seniority)
}

test.messages.nrow <- function() {
    checkEquals(nrow(enron$messages), 21635)
}

test.messages.unix.time <- function() {
    checkEquals(enron$messages$unix.time,
                as.integer(enron$messages$time))
}

test.messages.nreceiver <- function() {
    checkEquals(enron$messages$nreceiver,
                sapply(enron$messages$receiver.id, length))
}

test.messages.data <- function() {
    m <- enron$messages[c(1, 1971, 21635),]
    filename <- c("taylor-m/sent/11",
                  "shackleton-s/sent/1304",
                  "germany-c/inbox/23")
    time <- as.POSIXct(strptime(c("13 Nov 1998 04:07:00",
                                  "2 May 2000 04:00:00",
                                  "21 Jun 2002 13:37:34"),
                                format = '%d %b %Y %H:%M:%S', tz = 'UTC'))
    subject <- c("Cd$ CME letter",
                 "Doctor's appointment on Wednesday, May 3",
                 "Master Termination Log")
    sender.id <- c(138,
                   120,
                   92)
    receiver.id <- list(c(59),
                        c(138, 59, 130, 4),
                        c(39, 120, 4))

    checkEquals(m$filename, filename)
    checkEquals(m$time, time)
    checkEquals(m$subject, subject)
    checkEquals(m$sender.id, sender.id)
    checkEquals(m$receiver.id, receiver.id)
}

