#ifndef FIXTURES_RECV_MODEL_H
#define FIXTURES_RECV_MODEL_H

struct recv_model_fixture {
	struct recv_params params;
	size_t recv_trait_dim;
	size_t recv_tvar_dim;
	size_t dyad_trait_dim;
	size_t dyad_tvar_dim;
	int exclude_loops;
};

void recv_model_fixture_setup(struct recv_model_fixture	*f,
			      const struct design_fixture *r,
			      const struct design2_fixture *d);
void recv_model_fixture_set_exlude_loops(struct recv_model_fixture *f, int exclude_loops);
void recv_model_fixture_set_rand(struct recv_model_fixture *f);
void recv_model_fixture_teardown(struct recv_model_fixture *f);

void recv_model_test_setup(struct recv_model *m, struct design *r,
			   struct design2 *d,
			   const struct recv_model_fixture *f);
void recv_model_test_teardown(struct recv_model *m);


#endif /* FIXTURES_RECV_MODEL_H */
