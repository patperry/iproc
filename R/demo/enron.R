

# Restrict attention to messages with 5 or fewer receivers
message.data <- subset(enron.messages, sapply(receiver, length) <= 5)

# Restrict attention to 3 covariates
actor.data <- enron.actors[, c("department", "seniority", "gender")]

# fit a model with only static covariates
model0 <- recv.model(Mesg(time, sender, receiver) ~ s(.) * .,
                     message.data, actor.data)

# Show the fitted coefficients
summary(model0)


# Fit a model with dynamic covariates
intervals <- c(7.5 * 60 * 4^(1:7))
horizon <- min(message.data$time) + max(intervals)

model1 <- recv.model(Mesg(time, sender, receiver)
                     ~ s(.) * .
                     + isend(w) + irecv(w)
                     + nsend(ws) + nrecv(ws),
                     message.data, actor.data,
                     response.data = subset(message.data, time >= horizon),
                     w = max(intervals), ws = intervals)


# Extract the coefficients
beta1 <- coef(model1)
names1 <- names(beta1)

# Compute their standard errors
cov1 <- vcov(model1)
se1 <- sqrt(diag(cov1))

# Compute a table of the static coefficients
tab <- matrix(NA, 5, 4)
colnames(tab) <- c("departmentLegal", "departmentTrading", "seniorityJunior", "genderFemale")
rownames(tab) <- c("(Intercept)", paste("s(", colnames(tab), ")", sep=""))
tab[1,] <- beta1[colnames(tab)]
for (i in 2:5) {
    for (j in 1:4) {
        name <- paste(rownames(tab)[i], ":", colnames(tab)[j], sep="")
        tab[i,j] <- beta1[name]
    }
}

# Abbreviate the row and column names
colnames(tab) <- c("L", "T", "J", "F")
rownames(tab) <- c("1", "L", "T", "J", "F")
round(tab, 2)


# Plot the decays of the "nsend" and "nrecv" effects
beta1.nsend <- beta1[grep("nsend", names1)]
se1.nsend <- se1[grep("nsend", names1)]
beta1.nrecv <- beta1[grep("nrecv", names1)]
se1.nrecv <- se1[grep("nrecv", names1)]


par(mfrow=c(1, 2))

plot(log10(intervals), beta1.nsend)
segments(log10(intervals), beta1.nsend - se1.nsend,
         log10(intervals), beta1.nsend + se1.nsend)

plot(log10(intervals), beta1.nrecv)
segments(log10(intervals), beta1.nrecv - se1.nrecv,
         log10(intervals), beta1.nrecv + se1.nrecv)



# Fit a model with triadic covariates, but only for messages from
# sender 138 (all messages are used to compute the covariates)
model2 <- recv.model(Mesg(time, sender, receiver)
                     ~ s(.) * .
                     + isend(w) + irecv(w)
                     + nsend(ws) + nrecv(ws)
                     + nsend2(ws,ws) + nrecv2(ws,ws)
                     + nsib(ws,ws) + ncosib(ws,ws),
                     message.data, actor.data,
                     response.data = subset(message.data, time >= horizon & sender == 138),
                     w = max(intervals), ws = intervals, method="recv.frame")


# Note: the estimated parameters are not bias-corrected for multiple
# recipients, and the models are slightly different than those
# appearing in the JRSS paper.
