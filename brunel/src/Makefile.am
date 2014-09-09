AM_CPPFLAGS = -I$(top_builddir)/gl -I$(top_srcdir)/gl -DLOCALEDIR=\"$(localedir)\"

LIBS += @LTLIBMULTITHREAD@ @LTLIBINTL@

bin_PROGRAMS = brunel
brunel_SOURCES = main.c
brunel_LDFLAGS = $(ZLIB_LDFLAGS) $(HTSLIB_LDFLAGS) -static
brunel_CFLAGS = $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
brunel_LDADD = $(top_srcdir)/gl/libbrunel.la
