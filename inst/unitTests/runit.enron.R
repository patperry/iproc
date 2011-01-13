
require(RUnit)
require(iproc)

.setUp <- function() {
    data(enron)
}

test.employees.nrow <- function() {
    checkEquals(nrow(enron$employees), 156)
}

test.employees.data <- function() {
    e <- enron$employees[c(1,34,156),]
    name <- c("John Arnold", "Drew Fossum", "Phillip K. Allen")
    department <- factor(c("Trading", "Legal", "Trading"),
                         levels = c("Legal", "Other", "Trading"))
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

test.recipients.nrow <- function() {
    checkEquals(nrow(enron$recipients), 38388)
}
