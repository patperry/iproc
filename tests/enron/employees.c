#include "port.h"
#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <yajl/yajl_parse.h>
#include "coreutil.h"
#include "matrixutil.h"
#include "xalloc.h"
#include "enron.h"

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

static const char *ENRON_COHORT_NAMES[] = {
	"LJF",
	"LJM",
	"LSF",
	"LSM",
	"TJF",
	"TJM",
	"TSF",
	"TSM",
	"OJF",
	"OJM",
	"OSF",
	"OSM"
};

static ptrdiff_t enron_cohort(enum department_code d,
			      enum seniority_code s,
			      enum gender_code g)
{
	ptrdiff_t c = 0;

	switch(d) {
	case DEPARTMENT_OTHER:
		c += 4;
	case DEPARTMENT_TRADING:
		c += 4;
	case DEPARTMENT_LEGAL:
		break;
	case DEPARTMENT_NA:
		return -1;
	}

	switch(s) {
	case SENIORITY_SENIOR:
		c += 2;
	case SENIORITY_JUNIOR:
		break;
	case SENIORITY_NA:
		return -1;
	}

	switch(g) {
	case GENDER_MALE:
		c += 1;
	case GENDER_FEMALE:
		break;
	case GENDER_NA:
		return -1;
	}

	return c;
}

struct employee_parse {
	size_t nactor;
	size_t ncohort;
	size_t dim;
	double *traits_t;
	size_t *cohorts;
	size_t nactor_max;

	ptrdiff_t id;
	enum gender_code gender;
	enum seniority_code seniority;
	enum department_code department;
	enum employee_map_key map_key;
};

bool enron_legal[ENRON_NCOHORT] =
    { true, true, true, true, false, false, false, false, false, false, false,
false };
bool enron_trading[ENRON_NCOHORT] =
    { false, false, false, false, true, true, true, true, false, false, false,
false };
bool enron_other[ENRON_NCOHORT] =
    { false, false, false, false, false, false, false, false, true, true, true,
true };

bool enron_junior[ENRON_NCOHORT] =
    { true, true, false, false, true, true, false, false, true, true, false,
false };
bool enron_senior[ENRON_NCOHORT] =
    { false, false, true, true, false, false, true, true, false, false, true,
true };

bool enron_female[ENRON_NCOHORT] =
    { true, false, true, false, true, false, true, false, true, false, true,
false };
bool enron_male[ENRON_NCOHORT] =
    { false, true, false, true, false, true, false, true, false, true, false,
true };

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

	if (parse->nactor == parse->nactor_max) {
		size_t nactor_max = ARRAY_GROW(parse->nactor_max, SIZE_MAX); 
		parse->traits_t = xrealloc(parse->traits_t,
					  nactor_max
					  * parse->dim
					  * sizeof(parse->traits_t[0]));
		parse->cohorts = xrealloc(parse->cohorts,
					  nactor_max *
					  sizeof(parse->cohorts[0]));
		parse->nactor_max = nactor_max;
	}

	size_t id = parse->nactor;
	double *x = parse->traits_t + id * parse->dim;

	assert(parse->id == (ptrdiff_t)parse->nactor + 1);
	assert(parse->department != DEPARTMENT_NA);
	assert(parse->gender != GENDER_NA);
	assert(parse->seniority != SENIORITY_NA);

	parse->nactor++;
	parse->cohorts[id] = enron_cohort(parse->department, parse->seniority,
					  parse->gender);


	memset(x, 0, parse->dim * sizeof(x[0]));
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

	/* department * seniority * gender */
	if (parse->department == DEPARTMENT_LEGAL
	    && parse->seniority == SENIORITY_JUNIOR
	    && parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_LJF] = 1.0;
	if (parse->department == DEPARTMENT_TRADING
	    && parse->seniority == SENIORITY_JUNIOR
	    && parse->gender == GENDER_FEMALE)
		x[ENRON_TRAIT_TJF] = 1.0;


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

int enron_employees_init_fread(size_t *nactorp,
			       size_t **cohortsp, size_t *ncohortp, 
			       const char * const **cohort_namesp,
			       double **traitsp, size_t *ntraitp,
			       const char * const **trait_namesp,
			       FILE *stream)
{
	unsigned char fileData[65536];
	size_t rd;
	yajl_status stat;
	bool parse_ok = true;

	struct employee_parse parse;

	parse.dim = ENRON_NTRAIT;
	parse.ncohort = ENRON_NCOHORT;
	parse.nactor = 0;
	parse.traits_t = NULL;
	parse.cohorts = NULL;
	parse.nactor_max = 0;

	yajl_handle hand = yajl_alloc(&parse_callbacks, NULL, (void *)&parse);
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
		if (stat != yajl_status_ok)
			break;
	}

	stat = yajl_complete_parse(hand);

	if (stat != yajl_status_ok) {
		unsigned char *str = yajl_get_error(hand, 1, fileData, rd);
		fprintf(stderr, "%s", (const char *)str);
		yajl_free_error(hand, str);
		parse_ok = false;
	}

	yajl_free(hand);

	if (!parse_ok) {
		*nactorp = 0;
		*cohortsp = NULL;
		*ncohortp = 0;
		free(parse.cohorts);
		*cohort_namesp = NULL;
		*traitsp = NULL;
		free(parse.traits_t);
		*ntraitp = 0;
		*trait_namesp = NULL;
	} else {
		*nactorp = parse.nactor;
		*ncohortp = parse.ncohort;
		*cohortsp = xrealloc(parse.cohorts,
				     parse.nactor * sizeof(parse.cohorts[0]));
		*cohort_namesp = ENRON_COHORT_NAMES;

		size_t n = parse.nactor;
		size_t p = parse.dim;
		struct dmatrix traits_t = { parse.traits_t, MAX(1, p) };
		struct dmatrix traits = { xmalloc(n * p * sizeof(double)), MAX(1, n) };

		matrix_dtrans(p, n, &traits_t, &traits);
		free(traits_t.data);
		*traitsp = traits.data;
		*ntraitp = p;
		*trait_namesp = ENRON_TRAIT_NAMES;
	}

	return parse_ok;
}

int enron_employees_init(size_t *nactorp,
		         size_t **cohortsp, size_t *ncohortp, 
			 const char * const **cohort_namesp,
			 double **traitsp, size_t *ntraitsp,
			 const char * const **trait_namesp)
{
	FILE *f = fopen(ENRON_EMPLOYEES_FILE, "r");

	if (!f) {
		fprintf(stderr, "Couldn't open employees file '%s'\n",
			ENRON_EMPLOYEES_FILE);
		return errno;
	}

	if (!enron_employees_init_fread(nactorp, cohortsp, ncohortp,
				cohort_namesp, traitsp, ntraitsp,
				trait_namesp, f)) {
		fprintf(stderr, "Couldn't parse employees file '%s'\n",
			ENRON_EMPLOYEES_FILE);
		fclose(f);
		return EINVAL;
	}

	fclose(f);
	return 0;
}
