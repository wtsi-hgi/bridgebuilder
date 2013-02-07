bin_PROGRAMS = binnie
binnie_SOURCES = binnie.c
#binnie_LDFLAGS = @ZLIB_LIBS@ @HTSLIB_LIBS@ @LTLIBINTL@
binnie_LDFLAGS = @ZLIB_LIBS@ @HTSLIB_LIBS@ 
binnie_CFLAGS = @ZLIB_CFLAGS@ @HTSLIB_CFLAGS@
binnie_CPPFLAGS = -DLOCALEDIR=\"$(localedir)\"
#binnie_LDADD = $(LTLIBTHREAD) @LTLIBINTL@ $(top_srcdir)/gl/libbinnie.la
binnie_LDADD = $(LTLIBTHREAD) $(top_srcdir)/gl/libbinnie.la
