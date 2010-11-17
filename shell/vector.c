/* vector.c
 * --------
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <glib.h>
#include <libla/index.h>
#include <libla/vector.h>

#include <lua.h>
#include <lauxlib.h>


la_size
checksize (lua_State *L, int narg)
{
    la_index n = luaL_checkinteger(L, narg);

    if (n < 0)
        luaL_argerror(L, narg, "expected non-negative");

    return (la_size)n;
}

la_size
optsize (lua_State *L, int narg, la_size def)
{
    return luaL_opt(L, checksize, narg, def);
}

void
pushsize (lua_State *L, la_size n)
{
    lua_pushinteger(L, n);
}

la_index
checkindex (lua_State *L, la_size n, int narg)
{
    la_index i = luaL_checkinteger(L, narg) - 1;
    luaL_argcheck(L, 0 <= i && i < n, narg, "index out of range");
    return i;
}

la_index
optindex (lua_State *L, la_size n, int narg, la_index def)
{
    if (lua_isnoneornil(L, narg)) {
        return def;
    } else {
        return checkindex(L, n, narg);
    }
}

void
pushindex (lua_State *L, la_index i)
{
    lua_pushinteger(L, i + 1);
}

LAVector *
checkvector (lua_State *L, int narg)
{
    LAVector *x = *((LAVector **)luaL_checkudata(L, narg, "vector"));
    if (!x) luaL_argerror(L, narg,
                              "invalid vector; likely an expired view");
    return x;
}

static int
new (lua_State *L)
{
    la_size n = checksize(L, 1);
    LAVector **px = (LAVector **)lua_newuserdata(L, sizeof(LAVector *));
    luaL_getmetatable(L, "vector");
    lua_setmetatable(L, -2);
    *px = la_vector_new(n);
    return 1;
}

static int
gc (lua_State *L)
{
    LAVector *x = *(LAVector **)lua_touserdata(L, 1);
    la_vector_free(x);
    return 0;
}

static int
dim (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    la_size n = la_vector_dim(x);
    pushsize(L, n);
    return 1;
}

static int
get (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    la_size n = la_vector_dim(x);
    la_index i = checkindex(L, n, 2);
    double e = la_vector_get(x, i);
    lua_pushnumber(L, e);
    return 1;
}

static int
set (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    la_size n = la_vector_dim(x);
    la_index i = checkindex(L, n, 2);
    double e = lua_tonumber(L, 3);
    la_vector_set(x, i, e);
    return 0;
}

static int
set_all (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    double e = lua_tonumber(L, 2);
    la_vector_set_all(x, e);
    return 0;
}

static void
checkvector_pair (lua_State *L,
                  int narg1,
                  LAVector **pv1,
                  int narg2,
                  LAVector **pv2)
{
    g_assert(pv1);
    g_assert(pv2);

    *pv1 = checkvector(L, narg1);
    *pv2 = checkvector(L, narg2);

    if (la_vector_dim(*pv1) != la_vector_dim(*pv2))
        luaL_argerror(L, narg2, "inconsistent dimension");
}

static int
copy (lua_State *L)
{
    LAVector *x, *y;
    checkvector_pair (L, 1, &x, 2, &y);
    la_vector_copy(x, y);
    return 0;
}

static int
swap (lua_State *L)
{
    LAVector *x, *y;
    checkvector_pair (L, 1, &x, 2, &y);
    la_vector_swap(x, y);
    return 0;

}

static int
swap_elems (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    la_size n = la_vector_dim(x);
    la_index i = checkindex(L, n, 2);
    la_index j = checkindex(L, n, 3);
    la_vector_swap_elems(x, i, j);
    return 0;
}

static int
reverse (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    la_vector_reverse(x);
    return 0;
}


static int
add (lua_State *L)
{
    LAVector *x, *y;
    checkvector_pair(L, 1, &x, 2, &y);
    la_vector_add(x, y);
    return 0;
}

static int
sub (lua_State *L)
{
    LAVector *x, *y;
    checkvector_pair(L, 1, &x, 2, &y);
    la_vector_sub(x, y);
    return 0;
}

static int
mul (lua_State *L)
{
    LAVector *x, *y;
    checkvector_pair (L, 1, &x, 2, &y);
    la_vector_mul(x, y);
    return 0;
}

static int
div (lua_State *L)
{
    LAVector *x, *y;
    checkvector_pair (L, 1, &x, 2, &y);
    la_vector_div(x, y);
    return 0;
}

static int
scale (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    double e = lua_tonumber(L, 2);
    la_vector_scale(x, e);
    return 0;
}

static int
shift (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    double e = lua_tonumber(L, 2);
    la_vector_shift(x, e);
    return 0;
}

static int
dot (lua_State *L)
{
    LAVector *x, *y;
    checkvector_pair (L, 1, &x, 2, &y);
    double res = la_vector_dot(x, y);
    lua_pushnumber(L, res);
    return 1;
}

static int
norm (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    double res = la_vector_norm(x);
    lua_pushnumber(L, res);
    return 1;
}

static int
sum_abs (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    double res = la_vector_sum_abs(x);
    lua_pushnumber(L, res);
    return 1;
}

static int
max_abs (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    double res = la_vector_max_abs(x);
    lua_pushnumber(L, res);
    return 1;
}

static int
max_abs_index (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    if (la_vector_dim(x) == 0)
        luaL_argerror(L, 1, "dimension-0 vector");
    la_index res = la_vector_max_abs_index(x);
    pushindex(L, res);
    return 1;
}

static int
tostring (lua_State *L)
{
    LAVector *x = checkvector(L, 1);
    lua_pushfstring(L, "vector.new(%d)", la_vector_dim(x));
    return 1;
}

static const struct luaL_Reg vectorlib_f [] = {
    { "new", new },
    { "dim", dim },
    { "set_all", set_all },
    { "get", get },
    { "set", set },
    { "copy", copy },
    { "swap", swap },
    { "swap_elements", swap_elems },
    { "reverse", reverse },
    { "add", add },
    { "sub", sub },
    { "mul", mul },
    { "div", div },
    { "scale", scale },
    { "shift", shift },
    { "dot", dot },
    { "norm", norm },
    { "sum_abs", sum_abs },
    { "max_abs", max_abs },
    { "max_abs_index", max_abs_index },
    { NULL, NULL}
};

static const struct luaL_Reg vectorlib_m [] = {
    { "__tostring", tostring },
    { "__gc", gc },
    { NULL, NULL}
};

int
luaopen_vector (lua_State *L)
{
    luaL_newmetatable(L, "vector");

    /* metatable.__index = metatable */
    lua_pushnumber(L, -1);	/* duplicates the metatable */
    lua_setfield(L, -2, "__index");

    luaL_register(L, NULL, vectorlib_m);

    luaL_register(L, "vector", vectorlib_f);
    return 1;
}

