#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdint.h>
#include <libla/ieee754.h>


union double_uint64 {
    double d;
    uint64_t w;
};

int
la_identical (double x,
              double y)
{
    union double_uint64 ux = { x };
    union double_uint64 uy = { y };
    return ux.w == uy.w;
}
