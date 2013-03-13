SUBDIRS = gl src test # po
dist_doc_DATA = README.md

EXTRA_DIST = gl/m4/gnulib-cache.m4 $(top_srcdir)/.version

BUILT_SOURCES = $(top_srcdir)/.version

$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version

