#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <yajl/yajl_parse.h>
#include "coreutil.h"
#include "xalloc.h"
#include "enron.h"

enum message_map_key {
	MAP_KEY_NONE = -1,
	MAP_KEY_ID,
	MAP_KEY_TIME,
	MAP_KEY_SENDER,
	MAP_KEY_RECEIVER,
	MAP_KEY_OTHER
};

struct message_parse {
	size_t maxrecip;
	struct messages *messages;
	size_t id;
	double time;
	size_t sender;
	size_t *receiver;
	size_t nrecv, nrecv_max;
	intptr_t attr;
	enum message_map_key map_key;
	size_t map_depth, arr_depth;
};

static void message_parse_clear(struct message_parse *parse)
{
	parse->nrecv = 0;
	parse->id = SIZE_MAX;
	parse->sender = SIZE_MAX;
	parse->attr = 0;
	parse->time = NAN;
}

static void message_parse_init(struct message_parse *parse, struct messages *messages, size_t maxrecip)
{
	parse->receiver = NULL;
	parse->nrecv_max = 0;
	parse->messages = messages;
	parse->maxrecip = maxrecip;
	parse->map_key = MAP_KEY_NONE;
	parse->map_depth = 0;
	parse->arr_depth = 0;
	message_parse_clear(parse);
}

static void message_parse_grow_recv(struct message_parse *parse, size_t delta)
{
	size_t nmax = array_grow(parse->nrecv, parse->nrecv_max, delta, SIZE_MAX);
	if (nmax > parse->nrecv_max) {
		parse->receiver = xrealloc(parse->receiver, nmax * sizeof(parse->receiver[0]));
		parse->nrecv_max = nmax;
	}
}

static void message_parse_deinit(struct message_parse *parse)
{
	free(parse->receiver);
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
			parse->id = (size_t)(integerVal - 1);
			break;
		case MAP_KEY_TIME:
			parse->time = (double)(integerVal);
			break;
		case MAP_KEY_SENDER:
			if (integerVal <= 0) {
				fprintf(stderr, "non-positive sender: '%lld'\n", integerVal);
				return 0;
			}
			parse->sender = (size_t)(integerVal - 1);
			break;
		case MAP_KEY_RECEIVER:
			if (integerVal <= 0) {
				fprintf(stderr, "non-positive receiver: '%lld'\n", integerVal);
				return 0;
			}

			message_parse_grow_recv(parse, 1);
			parse->receiver[parse->nrecv++] = (size_t)(integerVal - 1);
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
			parse->time = doubleVal;
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

	if (parse->id == SIZE_MAX) {
		fprintf(stderr, "missing message id\n");
		return 0;
	}

	if (isnan(parse->time)) {
		fprintf(stderr, "missing time for message '%zu'\n", parse->id);
		return 0;
	}

	if (parse->sender == SIZE_MAX) {
		fprintf(stderr, "missing sender for message '%zu'\n", parse->id);
		return 0;
	}

	if (!parse->nrecv) {
		fprintf(stderr, "missing receiver for message '%zu'\n", parse->id);
		return 0;
	}

	if (!(parse->time >= messages_tlast(parse->messages))) {
		fprintf(stderr, "message '%zu' time not in sorted order\n", parse->id);
		return 0;
	}

	if (parse->nrecv <= parse->maxrecip || !parse->maxrecip) {
		messages_add(parse->messages,
			     parse->time,
			     parse->sender,
			     parse->receiver,
			     parse->nrecv,
			     parse->attr);
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

int enron_messages_init_fread(struct messages *messages, size_t maxrecip, FILE *stream)
{
	unsigned char fileData[65536];
	size_t rd;
	yajl_status stat;
	int parse_ok = 1;

	struct message_parse parse;

	messages_init(messages);
	message_parse_init(&parse, messages, maxrecip);

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
		messages_deinit(messages);
	}

	return parse_ok;
}

int enron_messages_init(struct messages *messages, size_t maxrecip)
{
	FILE *f = fopen(ENRON_MESSAGES_FILE, "r");

	if (!f) {
		fprintf(stderr, "Couldn't open messages file '%s'\n",
			ENRON_MESSAGES_FILE);
		return 0;
	}

	if (!enron_messages_init_fread(messages, maxrecip, f)) {
		fprintf(stderr, "Couldn't parse messages file '%s'\n",
			ENRON_MESSAGES_FILE);
		fclose(f);
		return 0;
	}

	fclose(f);
	return 1;
}


