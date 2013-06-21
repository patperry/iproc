
context("Mesg")

test_that("data argument works", {
    data <- data.frame(date = enron.messages$time)
    data$from <- enron.messages$sender
    data$to <- enron.messages$receiver
    data <- data[1:100,]

    m <- Mesg(date, from, to, data)

    expect_that(nrow(m), is_identical_to(nrow(data)))
    expect_that(m$time, is_identical_to(data$date))
    expect_that(m$sender, is_identical_to(data$from))
    expect_that(m$receiver, is_identical_to(data$to))
})


test_that("arguments can be read from parent environment", {
    date <- enron.messages$time[1:50]
    from <- enron.messages$sender[1:50]
    to <- enron.messages$receiver[1:50]

    m <- Mesg(date, from, to)

    expect_that(nrow(m), is_identical_to(length(date)))
    expect_that(m$time, is_identical_to(date))
    expect_that(m$sender, is_identical_to(from))
    expect_that(m$receiver, is_identical_to(to))
})


test_that("attributes work", {
    m <- Mesg(time, sender, receiver, enron.messages[1:60,],
              message.attr=subject)

    expect_that(m$message.attr, is_identical_to(enron.messages$subject[1:60]))
})
