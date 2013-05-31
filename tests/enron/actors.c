#include "port.h"
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <yajl/yajl_parse.h>
#include "coreutil.h"
#include "matrixutil.h"
#include "xalloc.h"
#include "enron/paths.h"
#include "enron/actors.h"


static const size_t ENRON_NTRAIT[ENRON_TERMS_MAX + 1] = { 0, 4, 9, 11 };

enum department_code { DEPARTMENT_NA =
	    -1, DEPARTMENT_LEGAL, DEPARTMENT_TRADING, DEPARTMENT_OTHER };
enum seniority_code { SENIORITY_NA = -1, SENIORITY_JUNIOR, SENIORITY_SENIOR };
enum gender_code { GENDER_NA = -1, GENDER_MALE, GENDER_FEMALE };
enum employee_map_key {
	MAP_KEY_NONE = -1,
	MAP_KEY_ID,
	MAP_KEY_DEPARTMENT,
	MAP_KEY_GENDER,
	MAP_KEY_SENIORITY,
	MAP_KEY_OTHER
};


static const char *ENRON_TRAIT_NAMES[] = {
	"Leg", "Trad", "Jun", "Fem",
	"Leg:Jun", "Trad:Jun", "Leg:Fem", "Trad:Fem", "Jun:Fem",
	"Leg:Jun:Fem", "Trad:Jun:Fem"
};

struct employee_parse {
	size_t nactor;
	size_t dim;
	double *traits;
	size_t nactor_max;
	size_t terms;

	ptrdiff_t id;
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
		parse->id = (ptrdiff_t)integerVal;
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

static int parse_string(void *ctx, const unsigned char *stringVal,
			size_t stringLen)
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

static int parse_map_key(void *ctx, const unsigned char *stringVal,
			 size_t stringLen)
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

	if (needs_grow(parse->nactor + 1, &parse->nactor_max)) {
		parse->traits = xrealloc(parse->traits,
					  parse->nactor_max
					  * parse->dim
					  * sizeof(parse->traits[0]));
	}

	size_t id = parse->nactor;
	double *x = parse->traits + id * parse->dim;

	assert(parse->id == (ptrdiff_t)parse->nactor + 1);
	assert(parse->department != DEPARTMENT_NA);
	assert(parse->gender != GENDER_NA);
	assert(parse->seniority != SENIORITY_NA);

	parse->nactor++;

	memset(x, 0, parse->dim * sizeof(x[0]));

	if (parse->terms < 1)
		goto out;

	/* department */
	if (parse->department == DEPARTMENT_LEGAL)
		x[ENRON_TRAIT_L] = 1.0;
	if (parse->department == DEPARTMENT_TRADING)
		x[ENRON_TRAIT_T] = 1.0;

	/* seniority */
	if (parse->seniority == SENIORITY_JUNIOR)
		x[ENRON_TRAIT_J] = 1.0;

	/* gender */
	if (parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_F] = 1.0;


	if (parse->terms < 2)
		goto out;

	/* department * seniority */
	if (parse->department == DEPARTMENT_LEGAL
	    && parse->seniority == SENIORITY_JUNIOR)
		x[ENRON_TRAIT_LJ] = 1.0;
	if (parse->department == DEPARTMENT_TRADING
	    && parse->seniority == SENIORITY_JUNIOR)
		x[ENRON_TRAIT_TJ] = 1.0;

	/* department * gender */
	if (parse->department == DEPARTMENT_LEGAL
	    && parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_LF] = 1.0;
	if (parse->department == DEPARTMENT_TRADING
	    && parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_TF] = 1.0;

	/* seniority * gender */
	if (parse->seniority == SENIORITY_JUNIOR
	    && parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_JF] = 1.0;

	if (parse->terms < 3)
		goto out;

	/* department * seniority * gender */
	if (parse->department == DEPARTMENT_LEGAL
	    && parse->seniority == SENIORITY_JUNIOR
	    && parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_LJF] = 1.0;
	if (parse->department == DEPARTMENT_TRADING
	    && parse->seniority == SENIORITY_JUNIOR
	    && parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_TJF] = 1.0;

out:
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


int enron_actors_init_fread(FILE *stream, size_t terms,
			    size_t *count, double **trait_x,
			    const char * const **trait_names, size_t *trait_dim)
{
	unsigned char fileData[65536];
	size_t rd;
	yajl_status stat;
	int parse_ok = 1;

	struct employee_parse parse;

	terms = MIN(terms, ENRON_TERMS_MAX);

	parse.dim = ENRON_NTRAIT[terms];
	parse.nactor = 0;
	parse.traits = NULL;
	parse.nactor_max = 0;
	parse.terms = terms;

	yajl_handle hand = yajl_alloc(&parse_callbacks, NULL, (void *)&parse);
	yajl_config(hand, yajl_allow_comments, 1);
	yajl_config(hand, yajl_dont_validate_strings, 1);

	for (;;) {
		rd = fread((void *)fileData, 1, sizeof(fileData) - 1, stream);

		if (rd == 0) {
			if (!feof(stream)) {
				fprintf(stderr, "error on file read.\n");
				parse_ok = 0;
			}
			break;
		}
		fileData[rd] = 0;

		stat = yajl_parse(hand, fileData, rd);
		if (stat != yajl_status_ok)
			break;
	}

	stat = yajl_complete_parse(hand);

	if (stat != yajl_status_ok) {
		unsigned char *str = yajl_get_error(hand, 1, fileData, rd);
		fprintf(stderr, "%s", (const char *)str);
		yajl_free_error(hand, str);
		parse_ok = 0;
	}

	yajl_free(hand);

	if (!parse_ok) {
		*count = 0;
		*trait_x = NULL;
		free(parse.traits);
		*trait_dim = 0;
		*trait_names = NULL;
	} else {
		*count = parse.nactor;
		*trait_x = parse.traits;
		*trait_dim = parse.dim;
		*trait_names = ENRON_TRAIT_NAMES;
	}

	return parse_ok;
}


int enron_actors_init(size_t terms, size_t *count, double **trait_x,
		      const char * const **trait_names, size_t *trait_dim)
{
	FILE *f = fopen(ENRON_ACTORS_FILE, "r");

	if (!f) {
		fprintf(stderr, "Couldn't open employees file '%s'\n",
			ENRON_ACTORS_FILE);
		return errno;
	}

	if (!enron_actors_init_fread(f, terms, count, trait_x, trait_names, trait_dim)) {
		fprintf(stderr, "Couldn't parse employees file '%s'\n",
			ENRON_ACTORS_FILE);
		fclose(f);
		return EINVAL;
	}

	fclose(f);
	return 0;
}

