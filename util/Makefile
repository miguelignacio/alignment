######################################################################
# quite generic Makefile
#
# mkw 2005-03-24

# Extension for files that go in a shared object (library)
SrcSuf        = cc
# Installation path for shared library
LIBDIR	      = $(HOME)/lib/

# Setup for ROOT libraries
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
#ROOTGLIBS   += -lX3d

# external includes / libraries
#CLHEP_DIR = /home/sebouh/CLHEP_install
#CLHEP_VERSION = 2.4.4.1
EXTINCDIR     = -I/usr/local/include/ -I$(CLHEP_DIR)/include
EXTLIBS       = -L$(CLHEP_DIR)/lib -lCLHEP-Matrix-$(CLHEP_VERSION) -lCLHEP-Vector-$(CLHEP_VERSION)

######################################################################
# You should not modify the lines below here
#

# Linux with gcc / g++ / g77
FC            = g77
FCFLAGS       = -I/cern/pro/include -fPIC -fno-automatic -fdollar-ok -fno-backslash \
		-finit-local-zero -fno-second-underscore -fugly-logint -ftypeless-boz
# -pg for gprof

CXX           = g++
CXXFLAGS      = -O2 -Wall -fPIC -fsignaling-nans -g # -DNDEBUG # -pg for gprof
CXXFLAGS     += $(EXTINCDIR)   # external header files should be here
CXXFLAGS     += $(ROOTCFLAGS)

LD            = g++
LDFLAGS       = -O2 -fsignaling-nans -g # -pg for gprof
LDFLAGS      += $(ROOTGLIBS)
LDFLAGS      += $(EXTLIBS)     # external libraries should be here
SOFLAGS       = -shared

SOURCES       = $(wildcard *.$(SrcSuf))
HEADERS       = $(SOURCES:.$(SrcSuf)=.h)
OBJECTS       = $(SOURCES:.$(SrcSuf)=.o)
DEPENDS       = $(SOURCES:.$(SrcSuf)=.d)
DEPENDS      += $(MAINSRC:.$(ExeSrcSuf)=.d)

DICT	      = Dict
SO	      = libUtilities.so

DATE	     := $(shell date "+%F")

######################################################################
# default targets

all:            $(SO)

clean:
		@rm -f $(OBJECTS) $(MAINOBJECTS) $(PROGRAM) $(SO) $(DICT).* $(DEPENDS)

distclean:      clean
		@rm -f *~ core *.def *.exp *.root *.ps .def gmon.out

tarball:
		tar cvjf kfa-$(DATE).tar.bz2 Makefile *.h *.cc *.cxx

install:	$(SO)
		cp -a $(SO) $(LIBDIR)

uninstall:
		test -x $(LIBDIR)/$(SO) && rm $(LIBDIR)/$(SO)

######################################################################
# rules

.SUFFIXES: .$(SrcSuf) .$(ExeSrcSuf) .C .f .o .so .d

.$(SrcSuf).o:
	$(CXX) $(CXXFLAGS) -c $<

.$(ExeSrcSuf).o:
	$(CXX) $(CXXFLAGS) -c $<

.C.o:
	$(CXX) $(CXXFLAGS) -c $<

.f.o:
	$(FC) $(FCFLAGS) -c $< -o $@

.$(SrcSuf).d:
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< > $@'

.$(ExeSrcSuf).d:
	$(SHELL) -ec '$(CXX) -M $(CXXFLAGS) $< > $@'

include $(DEPENDS)

######################################################################
# ROOT dictionary

$(DICT).C: $(HEADERS) LinkDef.h
	@echo "Generating dictionary $(DICT)..."
	@echo rootcint -f $@ -c $(EXTINCDIR) $+
	@rootcint -f $@ -c $(EXTINCDIR) $+

######################################################################
# targets

$(PROGRAM): 	$(SO) $(MAINOBJECTS)
	$(LD) $(LDFLAGS) $(EXTLIBS) $+ -o $@

$(SO):     $(OBJECTS) $(DICT).o
		$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@
		@echo "$@ done"

