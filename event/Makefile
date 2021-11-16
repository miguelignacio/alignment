######################################################################
# quite generic Makefile
#
# mkw 2005-03-24

# Name of the executable (code has to be in $(PROGRAM).$(ExeSrcSuf))
PROGRAM       = event
# Extension for file containing the main program
ExeSrcSuf     = cxx

# Extension for files that go in a shared object (library)
SrcSuf        = cc
# Installation path for shared library
LIBDIR	      = $(HOME)/lib/

# Setup for ROOT libraries
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
#ROOTGLIBS   += -lX3d

######################################################################
# You should not modify the lines below here
#
# main program
MAINSRC       = $(PROGRAM).$(ExeSrcSuf)
MAINOBJECTS   = $(MAINSRC:.$(ExeSrcSuf)=.o)

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
SOFLAGS       = -shared

SOURCES       = $(wildcard *.$(SrcSuf))
HEADERS       = $(SOURCES:.$(SrcSuf)=.h)
OBJECTS       = $(SOURCES:.$(SrcSuf)=.o)
DEPENDS       = $(SOURCES:.$(SrcSuf)=.d)
DEPENDS      += $(MAINSRC:.$(ExeSrcSuf)=.d)

DICT	      = Dict
SO	      = libAlignEvent.so

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

install:	$(PROGRAM)
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

align.test: $(MAINOBJECTS) $(OBJECTS) $(DICT).o
	$(LD) $(LDFLAGS) $+ -o $@

$(PROGRAM): 	$(SO) $(MAINOBJECTS)
	$(LD) $(LDFLAGS) $+ -o $@

$(SO):     $(OBJECTS) $(DICT).o
		$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@
		@echo "$@ done"

