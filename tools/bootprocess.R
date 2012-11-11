# tools/bootprocess.R
# -------------------

source("tools/utils.R")

pdffile <- "boot-resid.pdf"
pdfout <- TRUE

fit0 <- read.h5out("output/full/full.h5")
coefs0 <- fit0$coefficients
nobs0 <- fit0$observed.counts
nexp0 <- fit0$fitted.counts

nreps <- length(list.files("output/full", "fit.+[.]h5"))
coefs <- array(NA, c(dim(coefs0), nreps))
nobs <- array(NA, c(dim(nobs0), nreps))
nexp <- array(NA, c(dim(nexp0), nreps))


for (r in seq_len(nreps)) {
    fn <- paste("output/full/fit", r, ".h5", sep='')
    fit <- read.h5out(fn)
    coefs[,r] <- fit$coefficients
    #nobs[,,r] <- fit$observed.counts
    #nexp[,,r] <- fit$fitted.counts
    cat('.')
}

resid <- (nobs - nexp) / sqrt(nexp)
resid[nobs == 0 & nexp == 0] <- 0
x2 <- apply(resid^2, 3, sum)
resid0 <- (nobs0 - nexp0) / sqrt(nexp0)
resid0[nobs0 == 0 & nexp0 == 0] <- 0
x20 <- sum(resid0^2)



bias <- array(NA, dim(coefs))
for (r in seq_len(nreps)) {
    bias[,r] <- coefs[,r] - coefs0
}

cov <- get.cov(fit0)
se <- sqrt(pmax(0, diag(cov)))

bias.mean <- apply(bias, 1, mean)
bias.sd <- apply(bias, 1, sd)

save(bias.mean, bias.sd, se, bias, file="boot.rda")


par(las = 1)

ix <- seq_along(bias.mean)

nm <- names(coefs0)
ni <- length(fit0$intervals)
i1 <- seq(from=0, by = 1, length.out = ni)
i2 <- seq(from=0, by = 1, length.out = ni^2)

ix.recv <- c(match("IRecv", nm), match("NRecv[1]", nm) + i1)
ix.send <- c(match("ISend", nm), match("NSend[1]", nm) + i1)
ix.recv2 <- c(match("IRecv2", nm), match("NRecv2[1,1]", nm) + i2)
ix.send2 <- c(match("ISend2", nm), match("NSend2[1,1]", nm) + i2)
ix.sib <- c(match("ISib", nm), match("NSib[1,1]", nm) + i2)
ix.cosib <- c(match("ICosib", nm), match("NCosib[1,1]", nm) + i2)
ix <- c(ix.send, ix.recv, ix.send2, ix.recv2, ix.sib, ix.cosib)
ix.stat <- setdiff(seq_along(nm), ix)
ix <- c(ix.stat, ix)

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
