#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <setjmp.h>
#include <cmockery.h>
#include <yajl/yajl_parse.h> 

#include "actors.h"

#define EMPLOYEES_FILE "data/enron/employees.json"
#define EMPLOYEES_SIZE  156
#define EMPLOYEES_DIM   5

static struct matrix employees;
static struct matrix matrix;
static struct actors actors;

enum gender_code { GENDER_NA = -1, GENDER_MALE, GENDER_FEMALE };
enum seniority_code { SENIORITY_NA = -1, SENIORITY_JUNIOR, SENIORITY_SENIOR };
enum department_code { DEPARTMENT_NA = -1, DEPARTMENT_LEGAL, DEPARTMENT_TRADING, DEPARTMENT_OTHER };
enum employee_map_key {
	MAP_KEY_NONE = -1,
	MAP_KEY_ID,
	MAP_KEY_DEPARTMENT,
	MAP_KEY_GENDER,
	MAP_KEY_SENIORITY,
	MAP_KEY_OTHER
};

struct employee_parse {
	struct matrix *matrix;
	ssize_t id;
	enum gender_code gender;
	enum seniority_code seniority;
	enum department_code department;
	enum employee_map_key map_key;
};

static int parse_integer(void *ctx, long long integerVal)
{
	struct employee_parse *parse = ctx;

	switch (parse->map_key) {
	case MAP_KEY_ID:
		parse->id = (ssize_t)integerVal;
		break;
	case MAP_KEY_OTHER:
		break;
	default:
		assert(0);
		return 0;
	}
	parse->map_key = MAP_KEY_NONE;
	return 1;
}

static int parse_string(void *ctx, const unsigned char *stringVal, size_t stringLen)
{
	struct employee_parse *parse = ctx;
	const char *sstringVal = (const char *)stringVal;

	switch (parse->map_key) {
	case MAP_KEY_DEPARTMENT:
		if (strncmp("Legal", sstringVal, stringLen) == 0) {
			parse->department = DEPARTMENT_LEGAL;
		} else if (strncmp("Trading", sstringVal, stringLen) == 0) {
			parse->department = DEPARTMENT_TRADING;
		} else if (strncmp("Other", sstringVal, stringLen) == 0) {
			parse->department = DEPARTMENT_OTHER;
		} else {
			parse->department = DEPARTMENT_NA;
		}
		break;
	case MAP_KEY_GENDER:
		if (strncmp("Male", sstringVal, stringLen) == 0) {
			parse->gender = GENDER_MALE;
		} else if (strncmp("Female", sstringVal, stringLen) == 0) {
			parse->gender = GENDER_FEMALE;
		} else {
			parse->gender = GENDER_NA;
		}
		break;
	case MAP_KEY_SENIORITY:
		if (strncmp("Junior", sstringVal, stringLen) == 0) {
			parse->seniority = SENIORITY_JUNIOR;
		} else if (strncmp("Senior", sstringVal, stringLen) == 0) {
			parse->seniority = SENIORITY_SENIOR;
		} else {
			parse->seniority = SENIORITY_NA;
		}
		break;
	case MAP_KEY_OTHER:
		break;
	default:
		assert(0);
		return 0;
	}
	parse->map_key = MAP_KEY_NONE;
	return 1;
}

static int parse_start_map(void *ctx)
{
	struct employee_parse *parse = ctx;

	parse->id = -1;
	parse->gender = GENDER_NA;
	parse->seniority = SENIORITY_NA;
	parse->department = DEPARTMENT_NA;
	parse->map_key = MAP_KEY_NONE;
	return 1;
}

static int parse_map_key(void *ctx, const unsigned char *stringVal, size_t stringLen)  
{  
    struct employee_parse *parse = ctx;
    const char *sstringVal = (const char *)stringVal;

    if (strncmp("id", sstringVal, stringLen) == 0) {
	    parse->map_key = MAP_KEY_ID;
    } else if (strncmp("department", sstringVal, stringLen) == 0) {
	    parse->map_key = MAP_KEY_DEPARTMENT;
    } else if (strncmp("gender", sstringVal, stringLen) == 0) {
	    parse->map_key = MAP_KEY_GENDER;
    } else if (strncmp("seniority", sstringVal, stringLen) == 0) {
	    parse->map_key = MAP_KEY_SENIORITY;
    } else {
	    parse->map_key = MAP_KEY_OTHER;
    }

    return 1;
}

static int parse_end_map(void *ctx)
{
	struct employee_parse *parse = ctx;
	struct matrix *matrix = parse->matrix;
	ssize_t i = parse->id - 1;

	matrix_fill_row(matrix, i, 0.0);

	matrix_set(matrix, i, 0, 1.0); // intercept

	switch (parse->department) {
	case DEPARTMENT_LEGAL:
		matrix_set(matrix, i, 1, 1.0);
		break;
	case DEPARTMENT_TRADING:
		matrix_set(matrix, i, 2, 1.0);
		break;
	case DEPARTMENT_OTHER:
		break;
	default:
		assert(0);
		break;
	}

	switch (parse->gender) {
	case GENDER_FEMALE:
		matrix_set(matrix, i, 3, 1.0);
		break;
	case GENDER_MALE:
		break;
	default:
		assert(0);
		break;
	}

	switch (parse->seniority) {
	case SENIORITY_JUNIOR:
		matrix_set(matrix, i, 4, 1.0);
		break;
	case SENIORITY_SENIOR:
		break;
	default:
		assert(0);
		break;
	}

	return 1;
}


static yajl_callbacks parse_callbacks = {
	NULL,
	NULL,
	parse_integer,
	NULL,
	NULL,
	parse_string,
	parse_start_map,
	parse_map_key,
	parse_end_map,
	NULL,
	NULL
};

static bool employee_matrix_init_fread(struct matrix *x, FILE *stream)
{
	unsigned char fileData[65536];
	size_t rd;
	yajl_status stat;
	bool parse_ok = true;

	struct employee_parse parse;
	ssize_t n = EMPLOYEES_SIZE;
	ssize_t p = EMPLOYEES_DIM;
	
	parse.matrix = x;
	if (!matrix_init(parse.matrix, n, p))
		return false;

	yajl_handle hand = yajl_alloc(&parse_callbacks, NULL, (void *) &parse);
	yajl_config(hand, yajl_allow_comments, 1);
	yajl_config(hand, yajl_dont_validate_strings, 1);

	for (;;) {
		rd = fread((void *)fileData, 1, sizeof(fileData) - 1, stream);

		if (rd == 0) {
			if (!feof(stream)) {
				fprintf(stderr, "error on file read.\n");
				parse_ok = false;
			}
			break;
		}
		fileData[rd] = 0;

		stat = yajl_parse(hand, fileData, rd);
		if (stat != yajl_status_ok) break;
	}

	stat = yajl_complete_parse(hand);
	
	if (stat != yajl_status_ok) {
		unsigned char * str = yajl_get_error(hand, 1, fileData, rd);
		fprintf(stderr, "%s", (const char *) str);
		yajl_free_error(hand, str);
		parse_ok = false;
	}

	yajl_free(hand);

	if (!parse_ok) {
		matrix_deinit(parse.matrix);
	}

	return parse_ok;
}


static void enron_setup_fixture(void **state)
{
	print_message("Enron employees\n");
	print_message("---------------\n");

	FILE *f = fopen(EMPLOYEES_FILE, "r");
	if (!f) {
		print_error("Couldn't open employees file '%s'\n", EMPLOYEES_FILE);
		exit(-1);
	}
	if (!employee_matrix_init_fread(&employees, f)) {
		exit(-1);
	}
	fclose(f);
}

static void enron_teardown_fixture(void **state)
{
	matrix_deinit(&employees);
	print_message("\n\n");
}

static void enron_setup(void **state)
{
	matrix_init_copy(&matrix, &employees);
	actors_init(&actors, matrix_ncol(&matrix));

	ssize_t i, m = matrix_nrow(&matrix);
	struct vector row;

	vector_init(&row, matrix_ncol(&matrix));

	for (i = 0; i < m; i++) {
		matrix_get_row(&matrix, i, vector_front(&row));
		actors_add(&actors, &row);
	}

	vector_deinit(&row);
}

static void enron_teardown(void **state)
{
	matrix_deinit(&matrix);
	actors_deinit(&actors);
}

static void enron_test_size(void **state)
{
	assert_int_equal(actors_size(&actors), matrix_nrow(&matrix));
	assert_int_equal(actors_dim(&actors), matrix_ncol(&matrix));
	assert_int_equal(actors_cohorts_size(&actors), 12);
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(enron_test_size, enron_setup, enron_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}


