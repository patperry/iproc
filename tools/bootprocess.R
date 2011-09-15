# tools/bootprocess.R
# -------------------

source("tools/process.R")

fit0 <- fromJSON("output/fit-dynamic.json")
coefs0 <- get.coefs(fit0)

nreps <- length(list.files("output/boot-dynamic", "*.json"))
coefs <- array(NA, c(dim(coefs0), nreps))
dimnames(coefs) <- c(dimnames(coefs0), NULL)

for (r in seq_len(nreps)) {
    fn <- paste("output/boot-dynamic/fit", r, ".json", sep='')
    fit <- fromJSON(fn)
    coefs[,,r] <- get.coefs(fit)
    cat('.')
}

tx <- get.transform()
tcoefs0 <- coefs0 %*% tx

tcoefs <- array(NA, dim(coefs))
for (r in seq_len(nreps)) {
    tcoefs[,,r] <- coefs[,,r] %*% tx
}

tix0 <- rep(1:9, 10) + rep((1:10 - 1) * nrow(tcoefs0), each=9)
tix1 <- 12:nrow(tcoefs0)
tix <- c(tix0, tix1)

bias <- array(NA, dim(tcoefs))
for (r in seq_len(nreps)) {
    bias[,,r] <- tcoefs[,,r] - tcoefs0
}

cov <- get.cov(fit0)
vtx <- get.transform(TRUE)

ix0 <- rep(1:11, 12) + rep((1:12 - 1) * nrow(coefs0), each=11)

cov0 <- cov[ix0, ix0]
tcov0 <- vtx %*% cov0 %*% t(vtx)
cov1 <- cov[12:(nrow(cov)/12), 12:(nrow(cov)/12)]
se0 <- as.vector(matrix(sqrt(pmax(0, diag(tcov0))), 11, 12)[1:9, 1:10])
se1 <- sqrt(pmax(0, diag(cov1)))
se <- c(se0, se1)

bias.mean <- apply(bias, c(1,2), mean)[tix]
bias.sd <- apply(bias, c(1,2), sd)[tix]

save(bias.mean, bias.sd, se, bias, file="boot.rda")


par(las = 1)

ix <- seq_along(bias.mean)

color.segments <- 2
color.points <- 1

plot(range(ix),
     c(min((bias.mean - bias.sd) / se),
       max((bias.mean + bias.sd) / se)), 
     t = "n", axes = FALSE,
     xlab = "Coefficient Index",
     ylab = "Normalized Residual") 
#abline(v = nstatic + 1, lty = 2, col = "gray")
#xaxis.at <- c(1, 73, 145)
#axis(1, at = xaxis.at, labels = TRUE)
#axis(3, at = xaxis.at, labels = FALSE)
axis(2, labels = TRUE)
axis(4, labels = FALSE)
box()

segments(ix, ((bias.mean - bias.sd) / se)[ix],
         ix, ((bias.mean + bias.sd) / se)[ix],
         col = color.segments)
segments(ix - 0.45, ((bias.mean - bias.sd) / se)[ix],
         ix + 0.45, ((bias.mean - bias.sd) / se)[ix],
         col = color.segments)
segments(ix - 0.45, ((bias.mean + bias.sd) / se)[ix],
         ix + 0.45, ((bias.mean + bias.sd) / se)[ix],
         col = color.segments)

points(ix, (bias.mean / se)[ix],
       pch = 16, cex = 0.5, col = color.points)

