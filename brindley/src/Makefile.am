AM_CPPFLAGS = -I$(top_builddir)/gl -I$(top_srcdir)/gl -DLOCALEDIR=\"$(localedir)\"

LIBS += @LTLIBMULTITHREAD@ @LTLIBINTL@

bin_PROGRAMS = brindley
brindley_SOURCES = brindley.c brindley_coordmap.c brindley_log.c
brindley_LDFLAGS = $(ZLIB_LDFLAGS) $(HTSLIB_LDFLAGS) -static
brindley_CFLAGS = $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
brindley_LDADD = $(top_srcdir)/gl/libbrindley.la

noinst_HEADERS = brindley_coordmap.h