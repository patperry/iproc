ACLOCAL_AMFLAGS        = -I m4 --install

bin_PROGRAMS           = libla-shell
lib_LTLIBRARIES        = libla.la
noinst_LTLIBRARIES     = lua/liblua.la

include_HEADERS        = libla/index.h \
                         libla/vector.h

EXTRA_DIST             = lua/COPYRIGHT


libla_la_SOURCES       = libla/blas-private.h \
                         libla/index.h \
                         libla/vector.h \
                         libla/vector.c

libla_la_CPPFLAGS      = $(GLIB_CFLAGS) \
                         $(AM_CPPFLAGS)

libla_la_LIBADD        = $(BLAS_LIBS) \
                         $(GLIB_LIBS) \
                         $(LIBS) \
                         $(FLIBS)


libla_shell_SOURCES    = shell/main.c

libla_shell_CPPFLAGS   = -Ilua/src \
                         $(GLIB_CFLAGS) \
                         $(AM_CPPFLAGS)

libla_shell_LDADD      = $(top_builddir)/lua/liblua.la \
                         $(top_builddir)/libla.la \
                         $(GLIB_LIBS) \
                         $(READLINE_LIBS)


lua_liblua_la_SOURCES  = lua/src/lapi.h \
                         lua/src/lapi.c \
                         lua/src/lauxlib.h \
                         lua/src/lauxlib.c \
                         lua/src/lbaselib.c \
                         lua/src/lcode.h \
                         lua/src/lcode.c \
                         lua/src/ldblib.c \
                         lua/src/ldebug.h \
                         lua/src/ldebug.c \
                         lua/src/ldo.h \
                         lua/src/ldo.c \
                         lua/src/ldump.c \
                         lua/src/lfunc.h \
                         lua/src/lfunc.c \
                         lua/src/lgc.h \
                         lua/src/lgc.c \
                         lua/src/linit.c \
                         lua/src/liolib.c \
                         lua/src/llex.h \
                         lua/src/llex.c \
                         lua/src/llimits.h \
                         lua/src/lmathlib.c \
                         lua/src/lmem.h \
                         lua/src/lmem.c \
                         lua/src/loadlib.c \
                         lua/src/lobject.h \
                         lua/src/lobject.c \
                         lua/src/lopcodes.h \
                         lua/src/lopcodes.c \
                         lua/src/loslib.c \
                         lua/src/lparser.h \
                         lua/src/lparser.c \
                         lua/src/lstate.h \
                         lua/src/lstate.c \
                         lua/src/lstring.h \
                         lua/src/lstring.c \
                         lua/src/lstrlib.c \
                         lua/src/ltable.h \
                         lua/src/ltable.c \
                         lua/src/ltablib.c \
                         lua/src/ltm.h \
                         lua/src/ltm.c \
                         lua/src/lua.h \
                         lua/src/luaconf.h \
                         lua/src/lualib.h \
                         lua/src/lundump.h \
                         lua/src/lundump.c \
                         lua/src/lvm.h \
                         lua/src/lvm.c \
                         lua/src/lzio.h \
                         lua/src/lzio.c

lua_liblua_la_CPPFLAGS = -Ilua/src \
                         $(AM_CPPFLAGS)
