#ifndef _ENRON_H
#define _ENRON_H

#include <stdio.h>
#include "actors.h"
#include "matrix.h"
#include "messages.h"

#define PARENT_DIR "/Users/patperry/Projects/iproc/"
#define ENRON_EMPLOYEES_FILE PARENT_DIR"tests/enron/employees.json"
#define ENRON_MESSAGES_FILE  PARENT_DIR"tests/enron/messages.json"

bool enron_employees_init(struct actors *employees, struct matrix *traits);
bool enron_employees_init_fread(struct actors *employees, struct matrix *traits, FILE *stream);

bool enron_messages_init(struct messages *messages);
bool enron_messages_init_fread(struct messages *messages, FILE *stream);

#endif /* _ENRON_H */
