require(RJSONIO)
jmat <- function(j) { matrix(j$data, j$nrow, j$ncol) }

#fit <- fromJSON("~/Projects/iproc/fit.json")

get.traits <- function() {
    true <- 1
    false <- 0
    L <- c(true, true, true, true, false, false,
      	   false, false, false, false, false, false)
    T <- c(false, false, false, false, true, true,
      	   true, true, false, false, false, false)
    J <- c(true, true, false, false, true, true,
      	   false, false, true, true, false, false)
    F <- c(true, false, true, false, true, false,
      	   true, false, true, false, true, false)
    x <- cbind(1,
	       L, T, J, F,
	       L * J, T * J, L * F, T * F, J * F,
	       L * J * F, T * J * F)
    colnames(x) <- c("(intercept)",
    		     "L", "T", "J", "F",
    		     "L*J", "T*J", "L*F", "T*F", "J*F",
		     "L*J*F", "T*J*F")
    rownames(x) <- c("LJF", "LJM", "LSF", "LSM", "TJF",
    		     "TJM", "TSF", "TSM", "OJF", "OJM", "OSF", "OSM")
    x
}

get.transform <- function(vec = FALSE) {
    x <- get.traits()
    a <- solve(t(x))
    if (vec) {
        a <- kronecker(t(a), diag(ncol(x) - 1))
    }
    a
}


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

get.cov <- function(fit, duals=FALSE) {
        kkt <- get.kkt(fit)
        imat <- get.imat(fit, duals=TRUE)
        n <- sum(fit$count)

	s <- diag(kkt)
	s[s < 1e-8] <- 1.0
	s <- sqrt(s)

        skkt <- sweep(kkt / s, 2, s, "/")
	simat <- sweep(imat / s, 2, s, "/")

        qr.skkt <- qr(skkt)
        scov <- qr.solve(qr.skkt, t(qr.solve(qr.skkt, simat / n)))
        #scov <- solve(skkt, t(solve(skkt, simat / n)))
        cov <- sweep(scov / s, 2, s, "/")
	cov
}

get.se <- function(fit) {
       imatc <- lapply(fit$information, jmat)
       ce <- jmat(fit$constraints)
       dim <- nrow(imatc[[1]])
       nc <- length(imatc)
       ne <- ncol(ce)

       cov <- get.cov(fit, duals = TRUE)
       se2 <- diag(cov)
       se2[se2 < 0] <- 0
       se <- sqrt(se2)

       se.coefs <- matrix(se[1:(dim*nc)], dim, nc)
       se.duals <- se[dim*nc + 1:ne]
       list(coefs = se.coefs, duals = se.duals)
}

