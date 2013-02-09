LIBS += @LTLIBMULTITHREAD@

bin_PROGRAMS = binnie
binnie_SOURCES = binnie.c
binnie_LDFLAGS = $(ZLIB_LDFLAGS) $(HTSLIB_LDFLAGS) 
binnie_CFLAGS = $(ZLIB_CFLAGS) $(HTSLIB_CFLAGS)
binnie_CPPFLAGS = -DLOCALEDIR=\"$(localedir)\"
binnie_LDADD = $(top_srcdir)/gl/libbinnie.la 
