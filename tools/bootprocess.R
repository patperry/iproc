# tools/bootprocess.R
# -------------------

source("tools/process.R")

pdffile <- "figures/boot-resid.pdf"
pdfout <- TRUE

fit0 <- fromJSON("output/fit-dynamic.json")
coefs0 <- get.coefs(fit0)
nobs0 <- jmat(fit0$observed)
nexp0 <- jmat(fit0$expected)

nreps <- length(list.files("output/boot-dynamic", "*.json"))
coefs <- array(NA, c(dim(coefs0), nreps))
nobs <- array(NA, c(dim(nobs0), nreps))
nexp <- array(NA, c(dim(nexp0), nreps))

dimnames(coefs) <- c(dimnames(coefs0), NULL)

for (r in seq_len(nreps)) {
    fn <- paste("output/boot-dynamic/fit", r, ".json", sep='')
    fit <- fromJSON(fn)
    coefs[,,r] <- get.coefs(fit)
    nobs[,,r] <- jmat(fit$observed)
    nexp[,,r] <- jmat(fit$expected)
    cat('.')
}

resid <- (nobs - nexp) / sqrt(nexp)
resid[nobs == 0 & nexp == 0] <- 0
x2 <- apply(resid^2, 3, sum)
resid0 <- (nobs0 - nexp0) / sqrt(nexp0)
resid0[nobs == 0 & nexp0 == 0] <- 0
x20 <- sum(resid0^2)


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
ix.stat <- 1:90
ix.recv <- 79 + match("IRecv", rownames(coefs0)) + 0:7
ix.send <- 79 + match("ISend", rownames(coefs0)) + 0:7
ix.recv2 <- 79 + match("IRecv2", rownames(coefs0)) + 0:49
ix.send2 <- 79 + match("ISend2", rownames(coefs0)) + 0:49
ix.sib <- 79 + match("ISib", rownames(coefs0)) + 0:49
ix.cosib <- 79 + match("ICosib", rownames(coefs0)) + 0:49
ix <- c(ix.stat, ix.send, ix.recv, ix.send2, ix.recv2, ix.sib, ix.cosib)
sz <- c(length(ix.stat),
        length(ix.send), length(ix.recv),
        length(ix.send2), length(ix.recv2), length(ix.sib), length(ix.cosib))


purple <- rgb(158, 154, 200, max = 255)
orange <- rgb(253, 141,  60, max = 255)
color.points <- orange
color.segments <- purple


margin <- 0.82
padding <- margin/2

if (pdfout) {
    width.fig <- 8.3
    height.fig <- 0.5 * width.fig - margin + padding
    pdf(pdffile, width = width.fig, height = height.fig)
}

fin <- par()$fin
width = fin[1] - 2 * margin
height = fin[2] - margin - padding

par(las = 1)
par(plt = c(c(margin, margin + width)/fin[1],
            c(margin, margin + height)/fin[2]), new = FALSE)

plot(range(ix),
     c(min((bias.mean - bias.sd) / se),
       max((bias.mean + bias.sd) / se)), 
     xlim=c(0.5, length(ix) + 0.5),
     t = "n", axes = FALSE,
     xlab = "Term",
     ylab = "Normalized Residual") 
xtick <- c(0, cumsum(sz)) + 0.5
abline(v = xtick, lty = 1, col = "gray")
axis(1, at = xtick, labels = FALSE)
#axis(3, at = xaxis.at, labels = FALSE)
mtext("Static", 1, line=0.25, at=mean(xtick[1 + 0:1]), adj=1, las=2, cex=0.75)
mtext("Send", 1, line=0.25, at=mean(xtick[2 + 0:1]), adj=1, las=2, cex=0.75)
mtext("Receive", 1, line=0.25, at=mean(xtick[3 + 0:1]), adj=1, las=2, cex=0.75)
mtext("2-Send", 1, line=0.25, at=mean(xtick[4 + 0:1]), adj=1, las=2, cex=0.75)
mtext("2-Receive", 1, line=0.25, at=mean(xtick[5 + 0:1]), adj=1, las=2, cex=0.75)
mtext("Sibling", 1, line=0.25, at=mean(xtick[6 + 0:1]), adj=1, las=2, cex=0.75)
mtext("Cosibling", 1, line=0.25, at=mean(xtick[7 + 0:1]), adj=1, las=2, cex=0.75)

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


if (pdfout) dev.off()
