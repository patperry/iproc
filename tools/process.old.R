require(RJSONIO)
require(Matrix)

jmat <- function(j) { matrix(j$data, j$nrow, j$ncol) }
jsvec <- function(j) { spMatrix(j$dim, 1, j$pattern, rep(1L, j$count), j$data) }

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

get.constraints <- function(fit) {
    nce <- length(fit$constraints)
    dim <- fit$coefficients$nrow * fit$coefficients$ncol
    nnz <- lapply(fit$constraints, function(x) x$count)
    i <- rep(1:nce, nnz)
    j <- c(lapply(fit$constraints, function(x) x$pattern), recursive=TRUE)
    x <- c(lapply(fit$constraints, function(x) x$data), recursive=TRUE)
    a <- sparseMatrix(i, j, x = x, dims = c(nce, dim))
    b <- fit$constraint_values
    list(matrix = a, values = b)
}

get.imat <- function(fit, duals = FALSE) {
    dim <- fit$coefficients$nrow
    nc <- fit$coefficients$ncol
    nce <- length(fit$constraints)

    n <- fit$count
    ntot <- sum(n)
    imatc <- lapply(fit$information, jmat)


    m <- dim*nc
    if (duals)
        m <- m + nce

    imat <- sparseMatrix(integer(), integer(), x = numeric(), dims = c(m, m))
    for (ic in 1:nc) {
        off <- dim * (ic - 1)
        imat[off + 1:dim, off + 1:dim] <- (n[[ic]]/ntot) * imatc[[ic]]
    }
    imat
}

subspaces <- function(a, etol = 1e-8, vartol = 1e-8) {
    n <- nrow(a)
    s2 <- diag(a)
    s2[s2 < vartol] <- 1
    s <- sqrt(s2)
    as <- diag(1/s) %*% a %*% diag(1/s)
    e <- eigen(as, symmetric = TRUE)
    rank <- sum(e$values > etol)
    range <- e$vectors[,seq_len(rank), drop = FALSE]
    null <- e$vectors[,rank + seq_len(n - rank), drop = FALSE]
    list(scale = s, rank = rank, range = range, null = null)
}

get.kkt.inv <- function(fit, z) {
    dim <- fit$coefficients$nrow
    nc <- fit$coefficients$ncol
    nce <- length(fit$constraints)

    imatc <- lapply(fit$information, jmat)
    subsp <- lapply(imatc, subspaces)
    ce <- get.constraints(fit)$matrix

    xc <- lapply(seq_len(nc), function(i) z[(i-1)*dim + seq_len(dim)])
    y <- z[dim*nc + seq_len(nce)]


    imatc1 <- lapply(seq_len(nc), function(i) {
        s <- subsp[[i]]$scale
        u <- subsp[[i]]$range
        t(u) %*% diag(1/s) %*% imatc[[i]] %*% diag(1/s) %*% u
    })

    ce1 <- lapply(seq_len(nc), function(i) {
        s <- subsp[[i]]$scale
        u <- subsp[[i]]$range
        ce[,(i-1)*dim + seq_len(dim)] %*% diag(1/s) %*% u
    })
    ce1 <- do.call(cbind, lapply(ce1, as.matrix))

    ce2 <- lapply(seq_len(nc), function(i) {
        s <- subsp[[i]]$scale
        u <- subsp[[i]]$null
        ce[,(i-1)*dim + seq_len(dim)] %*% diag(1/s) %*% u
    })
    ce2 <- do.call(cbind, lapply(ce2, as.matrix))

    xc1 <- lapply(seq_len(nc), function(i) {
        s <- subsp[[i]]$scale
        u <- subsp[[i]]$range
        t(u) %*% diag(1/s) %*% xc[[i]]
    })

    xc2 <- lapply(seq_len(nc), function(i) {
        s <- subsp[[i]]$scale
        u <- subsp[[i]]$null
        t(u) %*% diag(1/s) %*% xc[[i]]
    })

    rank <- sapply(subsp, function(x) x$rank)
    cumrank <- c(0, cumsum(rank))

    ce1_vi <- lapply(seq_len(nc), function(i) {
        a <- ce1[,cumrank[i] + seq_len(rank[i]), drop = FALSE]
        t(solve(imatc1[[i]], t(a)))
    })
    ce1_vi <- do.call(cbind, lapply(ce1_vi, as.matrix))
    ce1_vi_ce1 <- ce1_vi %*% t(ce1)    

    v2 <- ce1_vi_ce1 + ce2 %*% t(ce2)
    s2 <- diag(v2)
    s2[s2 < 1e-8] <- 1
    s2 <- sqrt(s2)

}

get.kkt <- function(fit) {
    dim <- fit$coefficients$nrow
    nc <- fit$coefficients$ncol
    nce <- length(fit$constraints)

    n <- fit$count
    ntot <- sum(n)
    ce <- get.constraints(fit)$matrix

    m <- dim * nc + nce
    kkt <- get.imat(fit, duals=TRUE)
    if (nce > 0) {
        kkt[1:(dim*nc), dim * nc + 1:nce] <- t(ce)
        kkt[dim * nc + 1:nce, 1:(dim*nc)] <- ce
    }
    kkt
}

get.cov <- function(fit, duals=FALSE) {
    dim <- fit$coefficients$nrow
    nc <- fit$coefficients$ncol

    kkt <- get.kkt(fit)
    imat <- get.imat(fit, duals=TRUE)
    n <- sum(fit$count)

    s <- diag(kkt)
    s[s < 1e-8] <- 1.0
    s <- sqrt(s)

    skkt <- t(t(kkt) / s) / s
    simat <- t(t(imat) / s) / s
    skkti <- solve(skkt)
    scov <- skkti %*% (simat / n) %*% skkti
    cov <- t(t(scov) / s) / s

    if (!duals) {
        m <- dim * nc
        cov <- cov[seq_len(m), seq_len(m), drop = FALSE]
    }
    cov
}

get.se <- function(fit) {
    dim <- fit$coefficients$nrow
    nc <- fit$coefficients$ncol
    nce <- length(fit$constraints)

    cov <- get.cov(fit, duals = TRUE)
    se2 <- diag(cov)
    se2[se2 < 0] <- 0
    se <- sqrt(se2)

    se.coefs <- matrix(se[seq_len(dim*nc)], dim, nc)
    se.duals <- se[dim*nc + seq_len(nce)]
    list(coefs = se.coefs, duals = se.duals)
}

