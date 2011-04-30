#include "port.h"
#include <assert.h>
#include <math.h>
#include <yajl/yajl_parse.h>
#include "enron.h"

enum message_map_key {
	MAP_KEY_NONE = -1,
	MAP_KEY_ID,
	MAP_KEY_TIME,
	MAP_KEY_SENDER_ID,
	MAP_KEY_RECEIVER_ID,
	MAP_KEY_OTHER
};

struct message_parse {
	struct messages *messages;
	ssize_t id;
	double time;
	ssize_t sender_id;
	struct darray receiver_id;
	enum message_map_key map_key;
	bool multiple_receivers;
};

static int parse_integer(void *ctx, long long integerVal)
{
	struct message_parse *parse = ctx;
	double doubleVal = (double)integerVal;
	ssize_t ssizeVal = (ssize_t)integerVal;
	
	switch (parse->map_key) {
		case MAP_KEY_ID:
			if (ssizeVal <= 0) {
				fprintf(stderr, "non-positive message id: '%"SSIZE_FMT"'", ssizeVal);
				return 0;
			}
			parse->id = ssizeVal;
			break;
		case MAP_KEY_TIME:
			parse->time = doubleVal;
			break;
		case MAP_KEY_SENDER_ID:
			if (ssizeVal <= 0) {
				fprintf(stderr, "non-positive sender_id: '%"SSIZE_FMT"'", ssizeVal);
				return 0;
			}
			parse->sender_id = ssizeVal;
			break;
		case MAP_KEY_RECEIVER_ID:
			if (ssizeVal <= 0) {
				fprintf(stderr, "non-positive receiver_id: '%"SSIZE_FMT"'", ssizeVal);
				return 0;
			}
			
			if (!darray_push_back(&parse->receiver_id, &ssizeVal)) {
				fprintf(stderr, "not enough memory");
				return 0;
			}
			
			if (parse->multiple_receivers)
				return 1; // exit before clearing map_key
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

static int parse_double(void *ctx, double doubleVal)
{
	struct message_parse *parse = ctx;
	switch (parse->map_key) {
		case MAP_KEY_ID:
			fprintf(stderr, "non-integer message id: '%g'", doubleVal);
			return 0;
		case MAP_KEY_TIME:
			parse->time = doubleVal;
			break;
		case MAP_KEY_SENDER_ID:
			fprintf(stderr, "non-integer sender_id: '%g'", doubleVal);
			return 0;
		case MAP_KEY_RECEIVER_ID:
			fprintf(stderr, "non-integer receiver_id: '%g'", doubleVal);
			return 0;
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
	struct message_parse *parse = ctx;
	
	parse->id = -1;
	parse->time = NAN;
	parse->sender_id = -1;
	darray_clear(&parse->receiver_id);
	parse->map_key = MAP_KEY_NONE;
	parse->multiple_receivers = false;
	return 1;
}

static int parse_map_key(void *ctx, const unsigned char *stringVal, size_t stringLen)  
{  
	struct message_parse *parse = ctx;
	const char *sstringVal = (const char *)stringVal;
	
	if (strncmp("id", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_ID;
	} else if (strncmp("time", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_TIME;
	} else if (strncmp("sender_id", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_SENDER_ID;
	} else if (strncmp("receiver_id", sstringVal, stringLen) == 0) {
		parse->map_key = MAP_KEY_RECEIVER_ID;
	} else {
		parse->map_key = MAP_KEY_OTHER;
	}
	
	return 1;
}

static int parse_end_map(void *ctx)
{
	struct message_parse *parse = ctx;
	
	if (parse->id < 0) {
		fprintf(stderr, "missing message id");
		return 0;
	}
	
	if (isnan(parse->time)) {
		fprintf(stderr, "missing time for message '%" SSIZE_FMT "'",
			parse->id);
		return 0;
	}
	
	if (parse->sender_id < 0) {
		fprintf(stderr, "missing sender_id for message '%" SSIZE_FMT "'",
			parse->id);
		return 0;
	}

	if (darray_empty(&parse->receiver_id)) {
		fprintf(stderr, "missing receiver_id for message '%" SSIZE_FMT "'",
			parse->id);
		return 0;
	}

	messages_advance_to(parse->messages, parse->time);
	if (!messages_insertm(parse->messages,
			      parse->sender_id,
			      darray_front(&parse->receiver_id),
			      darray_size(&parse->receiver_id))) {
		fprintf(stderr, "not enough memory to insert message '%" SSIZE_FMT "'",
			parse->id);
		return 0;
	}

	return 1;
}

static int parse_begin_array(void *ctx)
{
	struct message_parse *parse = ctx;
	if (parse->map_key == MAP_KEY_RECEIVER_ID)
		parse->multiple_receivers = true;
	return 1;
}

static int parse_end_array(void *ctx)
{
	struct message_parse *parse = ctx;
	parse->map_key = MAP_KEY_NONE;
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

bool enron_messages_init_fread(struct messages *messages, FILE *stream)
{
	unsigned char fileData[65536];
	size_t rd;
	yajl_status stat;
	bool parse_ok = true;
	
	struct message_parse parse;
	
	parse.messages = messages;
	if (!messages_init(messages))
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
		messages_deinit(parse.messages);
	}
	
	return parse_ok;
}

bool enron_messages_init(struct messages *messages)
{
	FILE *f = fopen(ENRON_MESSAGES_FILE, "r");

	if (!f) {
		fprintf(stderr, "Couldn't open messages file '%s'\n",
			ENRON_MESSAGES_FILE);
		return false;
	}
	
	if (!enron_messages_init_fread(messages, f)) {
		fprintf(stderr, "Couldn't parse messages file '%s'\n",
			ENRON_MESSAGES_FILE);
		fclose(f);
		return false;
	}
	
	fclose(f);
	return true;
}


