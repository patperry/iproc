#ifndef _LA_CHECKUTILS_H
#define _LA_CHECKUTILS_H

#include <check.h>
#include <libla/ieee754.h>

#define ck_assert_feq(X, Y) \
    ck_assert_msg(la_identical(X, Y), \
                  "Assertion 'identical("#X", "#Y")' failed: " \
                  ""#X"==%.22f, "#Y"==%.22f", X, Y)



#endif /* _LA_CHECKUTILS_H */
