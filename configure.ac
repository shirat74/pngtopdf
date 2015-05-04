#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([pngtopdf], [0.1.0], [shunsaku.hirata74@gmail.com])
AC_CONFIG_SRCDIR([pngtopdf.cc])
AC_CONFIG_HEADERS([config.h])

AM_INIT_AUTOMAKE([pngtopdf], [0.1.0])

AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([atexit getcwd memmove memset pow setlocale sqrt strtol strtoul])

dnl zlib required
AC_ARG_WITH(zlib-include-dir,
 AS_HELP_STRING([--with-zlib-include-dir=DIR],
 [location of zlib headers]),,)
AC_ARG_WITH(zlib-lib-dir,
 AS_HELP_STRING([--with-zlib-lib-dir=DIR],
 [location of zlib library binary]),,)
if test "x$with_zlib_include_dir" != "x"; then
 CPPFLAGS="-I$with_zlib_include_dir $CPPFLAGS"
fi
if test "x$with_zlib_lib_dir" != "x"; then
 LDFLAGS="-L$with_zlib_lib_dir $LDFLAGS"
fi
AC_SEARCH_LIBS([deflate], [z], [has_zlib=yes], [has_zlib=no])
if test "x${has_zlib}" != xyes; then
   AC_MSG_ERROR([Unable to find zlib library.])
fi

dnl libpng required
AC_ARG_WITH(png-include-dir,
 AS_HELP_STRING([--with-png-include-dir=DIR],
 [location of libpng headers]),,)
AC_ARG_WITH(png-lib-dir,
 AS_HELP_STRING([--with-png-lib-dir=DIR],
 [location of libpng library binary]),,)
if test "x$with_png_include_dir" != "x"; then
  CPPFLAGS="-I$with_png_include_dir $CPPFLAGS"
fi
if test "x$with_png_lib_dir" != "x"; then
  LDFLAGS="-L$with_png_lib_dir $LDFLAGS"
fi
AC_CHECK_HEADER([png.h],
  AC_SEARCH_LIBS([png_read_png], [png],
  [has_png=yes], [has_png=no]), [has_png=no],[])
if test "x${has_png}" != xyes; then
   AC_MSG_ERROR([Unable to find libpng library.])
fi

dnl libqpdf required
AC_ARG_WITH(qpdf-include-dir,
 AS_HELP_STRING([--with-qpdf-include-dir=DIR],
 [location of libpqdf headers]),,)
AC_ARG_WITH(qpdf-lib-dir,
 AS_HELP_STRING([--with-qpdf-lib-dir=DIR],
 [location of qpdf library binary]),,)
if test "x$with_qpdf_include_dir" != "x"; then
 CPPFLAGS="-I$with_qpdf_include_dir $CPPFLAGS"
fi
if test "x$with_qpdf_lib_dir" != "x"; then
 LDFLAGS="-L$with_qpdf_lib_dir $LDFLAGS"
fi
dnl don't know how to test C++ library
AC_CHECK_HEADER([qpdf/QPDF.hh])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
