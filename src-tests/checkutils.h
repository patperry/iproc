#ifndef _IPROC_CHECKUTILS_H
#define _IPROC_CHECKUTILS_H

#include <check.h>
#include <iproc/ieee754.h>

#define ck_assert_feq(X, Y) \
    ck_assert_msg(iproc_identical(X, Y), \
                  "Assertion 'identical("#X", "#Y")' failed: " \
                  ""#X"==%.22f, "#Y"==%.22f", X, Y)



#endif /* _IPROC_CHECKUTILS_H */
