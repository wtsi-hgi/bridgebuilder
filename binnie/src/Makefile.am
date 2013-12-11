AM_CPPFLAGS = -I$(top_builddir)/gl -I$(top_srcdir)/gl -DLOCALEDIR=\"$(localedir)\"

LIBS += @LTLIBMULTITHREAD@ @LTLIBINTL@

bin_PROGRAMS = binnie liftover
binnie_SOURCES = binnie.c binnie_files.c binnie_log.c binnie_process.c binnie_coordmap.c
binnie_LDFLAGS = $(ZLIB_LDFLAGS) $(HTSLIB_LDFLAGS) -static
binnie_CFLAGS = $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
binnie_LDADD = $(top_srcdir)/gl/libbinnie.la
liftover_SOURCES = liftover.c binnie_coordmap.c binnie_log.c
liftover_LDFLAGS = $(ZLIB_LDFLAGS) $(HTSLIB_LDFLAGS) -static
liftover_CFLAGS = $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
liftover_LDADD = $(top_srcdir)/gl/libbinnie.la

noinst_HEADERS = binnie.h binnie_files.h binnie_log.h binnie_process.h binnie_coordmap.h