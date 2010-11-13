ACLOCAL_AMFLAGS              = -I m4 --install

bin_PROGRAMS                 = libla-shell
lib_LTLIBRARIES              = libla.la
noinst_LTLIBRARIES           = shell/lua/liblua.la

include_HEADERS              = libla/index.h \
                               libla/vector.h

EXTRA_DIST                   = shell/lua/COPYRIGHT


libla_la_SOURCES             = libla/blas-private.h \
                               libla/index.h \
                               libla/vector.h \
                               libla/vector.c

libla_la_CPPFLAGS            = $(GLIB_CFLAGS) \
                               $(AM_CPPFLAGS)

libla_la_LIBADD              = $(BLAS_LIBS) \
                               $(GLIB_LIBS) \
                               $(LIBS) \
                               $(FLIBS)


libla_shell_SOURCES          = shell/main.c

libla_shell_CPPFLAGS         = -Ishell/lua/src \
                               $(GLIB_CFLAGS) \
                               $(AM_CPPFLAGS)

libla_shell_LDADD            = $(top_builddir)/shell/lua/liblua.la \
                               $(top_builddir)/libla.la \
                               $(GLIB_LIBS) \
                               $(READLINE_LIBS)


shell_lua_liblua_la_SOURCES  = shell/lua/src/lapi.h \
                               shell/lua/src/lapi.c \
                               shell/lua/src/lauxlib.h \
                               shell/lua/src/lauxlib.c \
                               shell/lua/src/lbaselib.c \
                               shell/lua/src/lcode.h \
                               shell/lua/src/lcode.c \
                               shell/lua/src/ldblib.c \
                               shell/lua/src/ldebug.h \
                               shell/lua/src/ldebug.c \
                               shell/lua/src/ldo.h \
                               shell/lua/src/ldo.c \
                               shell/lua/src/ldump.c \
                               shell/lua/src/lfunc.h \
                               shell/lua/src/lfunc.c \
                               shell/lua/src/lgc.h \
                               shell/lua/src/lgc.c \
                               shell/lua/src/linit.c \
                               shell/lua/src/liolib.c \
                               shell/lua/src/llex.h \
                               shell/lua/src/llex.c \
                               shell/lua/src/llimits.h \
                               shell/lua/src/lmathlib.c \
                               shell/lua/src/lmem.h \
                               shell/lua/src/lmem.c \
                               shell/lua/src/loadlib.c \
                               shell/lua/src/lobject.h \
                               shell/lua/src/lobject.c \
                               shell/lua/src/lopcodes.h \
                               shell/lua/src/lopcodes.c \
                               shell/lua/src/loslib.c \
                               shell/lua/src/lparser.h \
                               shell/lua/src/lparser.c \
                               shell/lua/src/lstate.h \
                               shell/lua/src/lstate.c \
                               shell/lua/src/lstring.h \
                               shell/lua/src/lstring.c \
                               shell/lua/src/lstrlib.c \
                               shell/lua/src/ltable.h \
                               shell/lua/src/ltable.c \
                               shell/lua/src/ltablib.c \
                               shell/lua/src/ltm.h \
                               shell/lua/src/ltm.c \
                               shell/lua/src/lua.h \
                               shell/lua/src/luaconf.h \
                               shell/lua/src/lualib.h \
                               shell/lua/src/lundump.h \
                               shell/lua/src/lundump.c \
                               shell/lua/src/lvm.h \
                               shell/lua/src/lvm.c \
                               shell/lua/src/lzio.h \
                               shell/lua/src/lzio.c

shell_lua_liblua_la_CPPFLAGS = -Ishell/lua/src \
                               $(AM_CPPFLAGS)
