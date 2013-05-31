#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <yajl/yajl_parse.h>
#include "coreutil.h"
#include "xalloc.h"
#include "enron/paths.h"
#include "enron/messages.h"

enum message_map_key {
	MAP_KEY_NONE = -1,
	MAP_KEY_ID,
	MAP_KEY_TIME,
	MAP_KEY_SENDER,
	MAP_KEY_RECEIVER,
	MAP_KEY_OTHER
};

struct message_parse {
	/* parse parameters */
	size_t maxrecip;

	/* total number of senders and receivers */
	size_t nsend;
	size_t nrecv;

	/* all message data */
	double *time;
	size_t *from;
	size_t **to;
	size_t *nto;
	intptr_t *attr;
	size_t nmsg, nmsg_max;

	/* current message data */
	size_t cur_id;
	double cur_time;
	size_t cur_from;
	size_t *cur_to;
	size_t cur_nto, cur_nto_max;
	intptr_t cur_attr;

	/* parse internals */
	enum message_map_key map_key;
	size_t map_depth, arr_depth;
};

static void message_parse_clear(struct message_parse *parse)
{
	parse->cur_id = SIZE_MAX;
	parse->cur_time = NAN;
	parse->cur_from = SIZE_MAX;
	parse->cur_nto = 0;
	parse->cur_attr = 0;
}

static void message_parse_init(struct message_parse *parse, size_t maxrecip)
{
	size_t nmax = 65535;
	size_t nrecv_max = 255;

	parse->maxrecip = maxrecip;

	parse->nsend = 156;
	parse->nrecv = 156;

	parse->time = xmalloc(nmax * sizeof(double));
	parse->from = xmalloc(nmax * sizeof(size_t));
	parse->to = xmalloc(nmax * sizeof(size_t *));
	parse->nto = xmalloc(nmax * sizeof(size_t));
	parse->attr = xmalloc(nmax * sizeof(intptr_t));
	parse->nmsg = 0;
	parse->nmsg_max = nmax;
	
	parse->cur_to = xmalloc(nrecv_max * sizeof(size_t));
	parse->cur_nto_max = nrecv_max;

	parse->map_key = MAP_KEY_NONE;
	parse->map_depth = 0;
	parse->arr_depth = 0;

	message_parse_clear(parse);
}

static void message_parse_deinit(struct message_parse *parse)
{
	free(parse->cur_to);
}

static int parse_integer(void *ctx, long long integerVal)
{
	struct message_parse *parse = ctx;

	if (parse->map_depth != 1)
		return 1;

	switch (parse->map_key) {
		case MAP_KEY_ID:
			if (integerVal <= 0) {
				fprintf(stderr, "non-positive message id: '%lld'\n", integerVal);
				return 0;
			}
			parse->cur_id = (size_t)(integerVal - 1);
			break;
		case MAP_KEY_TIME:
			parse->cur_time = (double)(integerVal);
			break;
		case MAP_KEY_SENDER:
			if (integerVal <= 0) {
				fprintf(stderr, "non-positive sender: '%lld'\n", integerVal);
				return 0;
			}
			parse->cur_from = (size_t)(integerVal - 1);
			break;
		case MAP_KEY_RECEIVER:
			if (integerVal <= 0) {
				fprintf(stderr, "non-positive receiver: '%lld'\n", integerVal);
				return 0;
			} else if (parse->cur_nto == parse->cur_nto_max) {
				fprintf(stderr, "number of message recipients exceeded maximum ('%zd')\b",
					parse->cur_nto_max);
				return 0;
			}

			parse->cur_to[parse->cur_nto++] = (size_t)(integerVal - 1);
			break;
		case MAP_KEY_OTHER:
			break;
		default:
			assert(0);
			return 0;
	}
	return 1;
}

static int parse_double(void *ctx, double doubleVal)
{
	struct message_parse *parse = ctx;

	if (parse->arr_depth != 1 || parse->map_depth != 1)
		return 1;

	switch (parse->map_key) {
		case MAP_KEY_ID:
			fprintf(stderr, "non-integer message id: '%g'\n", doubleVal);
			return 0;
		case MAP_KEY_TIME:
			parse->cur_time = doubleVal;
			break;
		case MAP_KEY_SENDER:
			fprintf(stderr, "non-integer sender: '%g'\n", doubleVal);
			return 0;
		case MAP_KEY_RECEIVER:
			fprintf(stderr, "non-integer receiver: '%g'\n", doubleVal);
			return 0;
		case MAP_KEY_OTHER:
			break;
		default:
			assert(0);
			return 0;
	}
	return 1;
}

static int parse_start_map(void *ctx)
{
	struct message_parse *parse = ctx;

	parse->map_depth++;

	if (parse->arr_depth != 1 || parse->map_depth != 2)
		return 1;

	switch (parse->map_key) {
		case MAP_KEY_ID:
			fprintf(stderr, "non-integer message id\n");
			return 0;
		case MAP_KEY_TIME:
			fprintf(stderr, "non-numeric time\n");
			return 0;
		case MAP_KEY_SENDER:
			fprintf(stderr, "non-integer sender\n");
			return 0;
		case MAP_KEY_RECEIVER:
			fprintf(stderr, "non-integer receiver\n");
			return 0;
		default:
			break;
	}

	return 1;
}

static int parse_map_key(void *ctx, const unsigned char *stringVal, size_t stringLen)  
{
	struct message_parse *parse = ctx;
	const char *sstringVal = (const char *)stringVal;

	if (parse->arr_depth != 1 || parse->map_depth != 1)
		return 1;

	if (strncmp("id", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_ID;
	} else if (strncmp("time", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_TIME;
	} else if (strncmp("sender", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_SENDER;
	} else if (strncmp("receiver", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_RECEIVER;
	} else {
		parse->map_key = MAP_KEY_OTHER;
	}

	return 1;
}

static int parse_end_map(void *ctx)
{
	struct message_parse *parse = ctx;
	parse->map_depth--;

	if (parse->arr_depth != 1 || parse->map_depth != 0)
		return 1;

	if (parse->cur_id == SIZE_MAX) {
		fprintf(stderr, "missing message id\n");
		return 0;
	}

	if (isnan(parse->cur_time)) {
		fprintf(stderr, "missing time for message '%zu'\n", parse->cur_id);
		return 0;
	}

	if (parse->cur_from == SIZE_MAX) {
		fprintf(stderr, "missing sender for message '%zu'\n", parse->cur_id);
		return 0;
	}

	if (!parse->cur_nto) {
		fprintf(stderr, "missing receiver for message '%zu'\n", parse->cur_id);
		return 0;
	}

	if (!(!parse->nmsg || parse->cur_time >= parse->time[parse->nmsg])) {
		fprintf(stderr, "message '%zu' time not in sorted order\n", parse->cur_id);
		return 0;
	}

	if (parse->cur_nto <= parse->maxrecip || !parse->maxrecip) {
		assert(parse->nmsg < parse->nmsg_max);

		parse->time[parse->nmsg] = parse->cur_time;
		parse->from[parse->nmsg] = parse->cur_from;
		parse->to[parse->nmsg] = xmalloc(parse->cur_nto * sizeof(size_t));
		memcpy(parse->to[parse->nmsg], parse->cur_to, parse->cur_nto * sizeof(size_t));
		parse->nto[parse->nmsg] = parse->cur_nto;
		parse->attr[parse->nmsg] = parse->cur_attr;
		parse->nmsg++;
	}

	message_parse_clear(parse);

	return 1;
}

static int parse_begin_array(void *ctx)
{
	struct message_parse *parse = ctx;
	parse->arr_depth++;

	if (parse->arr_depth != 1 || parse->map_depth != 1)
		return 1;

	switch (parse->map_key) {
		case MAP_KEY_ID:
			fprintf(stderr, "non-integer message id\n");
			return 0;
		case MAP_KEY_TIME:
			fprintf(stderr, "non-numeric time\n");
			return 0;
		case MAP_KEY_SENDER:
			fprintf(stderr, "non-integer sender\n");
			return 0;
		default:
			break;
	}

	return 1;
}

static int parse_end_array(void *ctx)
{
	struct message_parse *parse = ctx;
	parse->arr_depth--;
	return 1;
}

static yajl_callbacks parse_callbacks = {
	NULL,
	NULL,
	parse_integer,
	parse_double,
	NULL,
	NULL,
	parse_start_map,
	parse_map_key,
	parse_end_map,
	parse_begin_array,
	parse_end_array
};

int enron_messages_init_fread(size_t *nsend, size_t *nrecv, double **time,
			      size_t **from, size_t ***to, size_t **nto,
			      intptr_t **attr, size_t *nmsg, size_t maxrecip, FILE *stream)
{
	unsigned char fileData[65536];
	size_t rd;
	yajl_status stat;
	int parse_ok = 1;

	struct message_parse parse;

	message_parse_init(&parse, maxrecip);

	yajl_handle hand = yajl_alloc(&parse_callbacks, NULL, (void *) &parse);
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
		if (stat != yajl_status_ok) break;
	}

	stat = yajl_complete_parse(hand);

	if (stat != yajl_status_ok) {
		unsigned char * str = yajl_get_error(hand, 1, fileData, rd);
		fprintf(stderr, "%s", (const char *) str);
		yajl_free_error(hand, str);
		parse_ok = 0;
	}

	yajl_free(hand);
	message_parse_deinit(&parse);

	if (!parse_ok) {
		size_t i, n = parse.nmsg;
		for (i = 0; i < n; i++) {
			free(parse.to[i]);
		}
		free(parse.attr);
		free(parse.nto);
		free(parse.to);
		free(parse.from);
		free(parse.time);
	}

	*nsend = parse.nsend;
	*nrecv = parse.nrecv;
	*time = parse.time;
	*from = parse.from;
	*to = parse.to;
	*nto = parse.nto;
	*attr = parse.attr;
	*nmsg = parse.nmsg;

	return parse_ok;
}

int enron_messages_init(size_t *nsend, size_t *nrecv, double **time,
			size_t **from, size_t ***to, size_t **nto,
			intptr_t **attr, size_t *nmsg, size_t maxrecip)
{
	FILE *f = fopen(ENRON_MESSAGES_FILE, "r");

	if (!f) {
		fprintf(stderr, "Couldn't open messages file '%s'\n",
			ENRON_MESSAGES_FILE);
		return 0;
	}

	if (!enron_messages_init_fread(nsend, nrecv, time, from, to, nto, attr, nmsg, maxrecip, f)) {
		fprintf(stderr, "Couldn't parse messages file '%s'\n",
			ENRON_MESSAGES_FILE);
		fclose(f);
		return 0;
	}

	fclose(f);
	return 1;
}


