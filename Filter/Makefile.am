AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(OFFLINE_MAIN)/include -I$(ROOTSYS)/include

lib_LTLIBRARIES = \
  libdnpfilter.la 

noinst_HEADERS = \
  dnpfilterLinkDef.h 

include_HEADERS = \
  dnpfilter.h

libdnpfilter_la_SOURCES = \
  dnpfilter.cxx \
  dnpfilter_Dict.cxx

libdnpfilter_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -lphool \
  -lpicodst_object \
  -lSubsysReco \
  -lmutgeom \
  -lmutoo_core \
  -lmutoo_base \
  -lmutoo_display \
  -lmutoo_interface \
  -lfvtxoo_interface \
  -lrpcoo_interface \
  -lmuioo_interface \
  -lMWG_interface \
  -lMWGOO \
  -lMWG_tools \
  -lreactionplane \
  -lMuidTriggerEmulator \
  -lfvtxgeom

# Rule for generating CINT dictionaries from class headers.
dnpfilter_Dict.cxx: \
  dnpfilter.h \
  dnpfilterLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

