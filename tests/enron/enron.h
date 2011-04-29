#ifndef _ENRON_H
#define _ENRON_H

#include <stdio.h>
#include "actors.h"
#include "matrix.h"

#define ENRON_EMPLOYEES_FILE "tests/enron/employees.json"

bool enron_employees_init(struct actors *employees);
bool enron_employee_matrix_init(struct matrix *employees);
bool enron_employee_matrix_init_fread(struct matrix *employees, FILE *stream);

#endif /* _ENRON_H */
