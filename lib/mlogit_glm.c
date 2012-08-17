


#if 0


void mlogit_mean_init(struct mlogit_mean *m, size_t dim, const double *mean0)
{
	m->dim = dim;
	m->mean = xmalloc(dim * sizeof(*m->mean));
	m->xbuf = xmalloc(dim * sizeof(*m->xbuf));
	
	memcpy(m->mean, mean0, dim * sizeof(*m->mean));
}


void mlogit_mean_deinit(struct mlogit_mean *m)
{
	free(m->xbuf);
	free(m->mean);
}


void mlogit_mean_update(struct mlogit_mean *m, const struct mlogit *mlogit,
			const double *x1, const double *dx,
			const struct vpattern *ix)
{
	size_t n = m->dim;
	double *buf = m->xbuf;
	double eta0 = mlogit->eta0;
	double eta_max = mlogit->eta_max;
	double phi = mlogit->phi;
	double expm1_deta = mlogit->expm1_deta;
	
	// buf := expm1(deta) * (x1 - mean0)
	memcpy(buf, x1, n * sizeof(*buf));
	blas_daxpy(n, -1.0, m->mean, 1, buf, 1); 
	blas_dscal(n, expm1_deta, buf, 1);
	
	// buf += dx
	if (ix) {
		sblas_daxpyi(1.0, dx, ix, buf);
	} else {
		blas_daxpy(n, 1.0, dx, 1, buf, 1);
	}
	
	blas_daxpy(n, exp(eta0 - eta_max - phi), buf, 1, m->mean, 1);
}

#endif


