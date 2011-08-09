#include "port.h"
#include <assert.h>
#include <yajl/yajl_parse.h>
#include "strata.h"
#include "enron.h"

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
	struct strata strata;
	struct actors *actors;
	struct vector traits;
	ssize_t id;
	enum gender_code gender;
	enum seniority_code seniority;
	enum department_code department;
	enum employee_map_key map_key;
};


bool enron_legal[ENRON_NCOHORT]   = { true, true, true, true, false, false, false, false, false, false, false, false };
bool enron_trading[ENRON_NCOHORT] = { false, false, false, false, true, true, true, true, false, false, false, false };
bool enron_other[ENRON_NCOHORT]   = { false, false, false, false, false, false, false, false, true, true, true, true };

bool enron_junior[ENRON_NCOHORT] = { true, true, false, false, true, true, false, false, true, true, false, false };
bool enron_senior[ENRON_NCOHORT] = { false, false, true, true, false, false, true, true, false, false, true, true };

bool enron_female[ENRON_NCOHORT] = { true, false, true, false, true, false, true, false, true, false, true, false };
bool enron_male[ENRON_NCOHORT]   = { false, true, false, true, false, true, false, true, false, true, false, true };



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
	struct vector *traits = &parse->traits;
	
	assert(parse->id - 1 == actors_count(parse->actors));
	assert(parse->department != DEPARTMENT_NA);
	assert(parse->gender != GENDER_NA);
	assert(parse->seniority != SENIORITY_NA);
	
	vector_fill(traits, 0.0);
	
	/* department */
	if (parse->department == DEPARTMENT_LEGAL)
		vector_set_item(traits, ENRON_TRAIT_L, 1.0);
	if (parse->department == DEPARTMENT_TRADING)
		vector_set_item(traits, ENRON_TRAIT_T, 1.0);
	
	/* seniority */
	if (parse->seniority == SENIORITY_JUNIOR)
		vector_set_item(traits, ENRON_TRAIT_J, 1.0);

	/* gender */
	if (parse->gender == GENDER_FEMALE)
		vector_set_item(traits, ENRON_TRAIT_F, 1.0);
	
	/* department * seniority */
	if (parse->department == DEPARTMENT_LEGAL
	    && parse->seniority == SENIORITY_JUNIOR)
		vector_set_item(traits, ENRON_TRAIT_LJ, 1.0);
	if (parse->department == DEPARTMENT_TRADING
	    && parse->seniority == SENIORITY_JUNIOR)
		vector_set_item(traits, ENRON_TRAIT_TJ, 1.0);
	
	/* department * gender */
	if (parse->department == DEPARTMENT_LEGAL
	    && parse->gender == GENDER_FEMALE)
		vector_set_item(traits, ENRON_TRAIT_LF, 1.0);
	if (parse->department == DEPARTMENT_TRADING
	    && parse->gender == GENDER_FEMALE)
		vector_set_item(traits, ENRON_TRAIT_TF, 1.0);

	/* seniority * gender */
	if (parse->seniority == SENIORITY_JUNIOR
	    && parse->gender == GENDER_FEMALE)
		vector_set_item(traits, ENRON_TRAIT_JF, 1.0);

	/* department * seniority * gender */
	if (parse->department == DEPARTMENT_LEGAL
	    && parse->seniority == SENIORITY_JUNIOR	    
	    && parse->gender == GENDER_FEMALE)
		vector_set_item(traits, ENRON_TRAIT_LJF, 1.0);
	if (parse->department == DEPARTMENT_TRADING
	    && parse->seniority == SENIORITY_JUNIOR	    
	    && parse->gender == GENDER_FEMALE)
		vector_set_item(traits, ENRON_TRAIT_TJF, 1.0);

	ssize_t cohort = strata_add(&parse->strata, traits);
	assert(cohort < ENRON_NCOHORT);
	actors_add(parse->actors, cohort);
	
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

static void traits_init(struct matrix *traits)
{
	matrix_init(traits, ENRON_NCOHORT, ENRON_NTRAIT);
	ssize_t ic, nc = ENRON_NCOHORT;

	for (ic = 0; ic < nc; ic++) {
		if (enron_legal[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_L, 1.0);
		if (enron_trading[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_T, 1.0);
		if (enron_junior[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_J, 1.0);
		if (enron_female[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_F, 1.0);
		
		if (enron_legal[ic] && enron_junior[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_LJ, 1.0);
		if (enron_trading[ic] && enron_junior[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_TJ, 1.0);
		if (enron_legal[ic] && enron_female[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_LF, 1.0);
		if (enron_trading[ic] && enron_female[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_TF, 1.0);
		if (enron_junior[ic] && enron_female[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_JF, 1.0);

		if (enron_legal[ic] && enron_junior[ic] & enron_female[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_LJF, 1.0);
		if (enron_trading[ic] && enron_junior[ic] & enron_female[ic])
			matrix_set_item(traits, ic, ENRON_TRAIT_TJF, 1.0);
	}
}

bool enron_employees_init_fread(struct actors *employees, struct matrix *traits, FILE *stream)
{
	unsigned char fileData[65536];
	size_t rd;
	yajl_status stat;
	bool parse_ok = true;
	
	struct employee_parse parse;
	ssize_t p = ENRON_NTRAIT;
	
	strata_init(&parse.strata, p);
	vector_init(&parse.traits, p);
	actors_init(employees);
	parse.actors = employees;
	
	yajl_handle hand = yajl_alloc(&parse_callbacks, NULL, (void *) &parse);
	yajl_config(hand, yajl_allow_comments, 1);
	yajl_config(hand, yajl_dont_validate_strings, 1);
	
	
	traits_init(traits);
	ssize_t ic, nc = ENRON_NCOHORT;
	for (ic = 0; ic < nc; ic++) {
		matrix_get_row(traits, ic, vector_to_ptr(&parse.traits));
		strata_add(&parse.strata, &parse.traits);
		assert(strata_count(&parse.strata) == ic + 1);
		actors_add_cohort(parse.actors);
	}
	
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
		actors_deinit(parse.actors);
	} else {
		const struct vector *level = strata_levels(&parse.strata);
		ssize_t il, nl = strata_count(&parse.strata);

		matrix_init(traits, nl, strata_dim(&parse.strata));
		for (il = 0; il < nl; il++) {
			matrix_set_row(traits, il, vector_to_ptr(&level[il]));
		}
	}
	strata_deinit(&parse.strata);
	
	return parse_ok;
}

bool enron_employees_init(struct actors *employees, struct matrix *traits)
{
	FILE *f = fopen(ENRON_EMPLOYEES_FILE, "r");

	if (!f) {
		fprintf(stderr, "Couldn't open employees file '%s'\n",
			ENRON_EMPLOYEES_FILE);
		return false;
	}
	
	if (!enron_employees_init_fread(employees, traits, f)) {
		fprintf(stderr, "Couldn't parse employees file '%s'\n",
			ENRON_EMPLOYEES_FILE);
		fclose(f);
		return false;
	}
	
	fclose(f);
	return true;
}
