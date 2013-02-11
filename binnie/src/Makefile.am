AM_CPPFLAGS = -I$(top_builddir)/gl -I$(top_srcdir)/gl -DLOCALEDIR=\"$(localedir)\"

LIBS += @LTLIBMULTITHREAD@ @LTLIBINTL@

bin_PROGRAMS = binnie
binnie_SOURCES = binnie.c binnie_files.c binnie_log.c
binnie_LDFLAGS = $(ZLIB_LDFLAGS) $(HTSLIB_LDFLAGS) -static
binnie_CFLAGS = $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
binnie_LDADD = $(top_srcdir)/gl/libbinnie.la 
