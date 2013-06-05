#ifndef FIXTURES_SEND_MODEL_H
#define FIXTURES_SEND_MODEL_H

struct send_model_fixture {
	struct send_params params;
};

void send_model_fixture_setup(struct send_model_fixture	*f,
			      const struct design_fixture *d);
void send_model_fixture_teardown(struct send_model_fixture *f);

void send_model_test_setup(struct send_model *m, struct design *d,
			   const struct send_model_fixture *f);
void send_model_test_teardown(struct send_model *m);


#endif /* FIXTURES_SEND_MODEL_H */