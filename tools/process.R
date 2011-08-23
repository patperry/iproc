require(RJSONIO)
jmat <- function(j) { matrix(j$data, j$nrow, j$ncol) }

#fit <- fromJSON("~/Projects/iproc/fit.json")

get.coefs <- function(fit) {
    coefs <- jmat(fit$coefficients)
    rownames(coefs) <- fit$variate_names
    colnames(coefs) <- fit$cohort_names
    coefs
}

get.duals <- function(fit) {
    duals <- fit$duals
    names(duals) <- fit$constraint_names
    duals
}

get.imat <- function(fit, duals = FALSE) {
        n <- fit$count
        ntot <- sum(n)
        imatc <- lapply(fit$information, jmat)
        ce <- jmat(fit$constraints)
        dim <- nrow(imatc[[1]])
        nc <- length(imatc)
        ne <- ncol(ce)

        m <- dim*nc
        if (duals)
            m <- m + ne

        imat <- matrix(0, m, m)
        for (ic in 1:nc) {
            off <- dim * (ic - 1)
            imat[off + 1:dim, off + 1:dim] <- (n[[ic]]/ntot) * imatc[[ic]]
        }
        imat
}

get.kkt <- function(fit) {
        n <- fit$count
        ntot <- sum(n)
        ce <- jmat(fit$constraints)

        imatc <- lapply(fit$information, jmat)
        dim <- nrow(imatc[[1]])
        nc <- length(imatc)
        ne <- ncol(ce)

        m <- dim * nc + ne
        kkt <- matrix(0, m, m)
        kkt[1:(dim*nc), 1:(dim*nc)] <- get.imat(fit, duals=FALSE)
        if (ne > 0) {
                kkt[1:(dim*nc), dim * nc + 1:ne] <- ce
                kkt[dim * nc + 1:ne, 1:(dim*nc)] <- t(ce)
        }
        kkt
}

get.cov <- function(fit) {
        kkt <- get.kkt(fit)
        imat <- get.imat(fit, duals=TRUE)
        n <- sum(fit$count)

	s <- diag(kkt)
	s[s < 1e-8] <- 1.0
	s <- sqrt(s)

	skkt <- diag(1/s) %*% kkt %*% diag(1/s)
	simat <- diag(1/s) %*% imat %*% diag(1/s)

        scov <- solve(skkt, t(solve(skkt, simat / n)))
        cov <- diag(1/s) %*% scov %*% diag(1/s)
	cov
}

get.se <- function(fit) {
       imatc <- lapply(fit$information, jmat)
       ce <- jmat(fit$constraints)
       dim <- nrow(imatc[[1]])
       nc <- length(imatc)
       ne <- ncol(ce)

       cov <- get.cov(fit)
       se2 <- diag(cov)
       se2[se2 < 0] <- 0
       se <- sqrt(se2)

       se.coefs <- matrix(se[1:(dim*nc)], dim, nc)
       se.duals <- se[dim*nc + 1:ne]
       list(coefs = se.coefs, duals = se.duals)
}

