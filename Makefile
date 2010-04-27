# This is a -*- Makefile -*- for the ublasJama Library
# Warning: this file contains tabs which cannot be converted to spaces

CC = gcc
CXX = g++
LD = g++
LDFLAGS = -g -pthread 

CFLAGS_OPT=-Wall -g -O3 -ftree-vectorize -msse3 -mssse3 -ffast-math

# Create dependencies
MAKEDEPEND = gcc -M $(CPPFLAGS) -o $*.d $<

# for gcov profiling add:
# -fprofile-arcs -ftest-coverage

BOOST_PATH = /opt/local/include
BOOST_LIBS = -L/opt/local/lib
BOOST_CPPFLAGS = -I$(BOOST_PATH)

LIBS =

CPPFLAGS = -I. $(BOOST_CPPFLAGS)

CXXFLAGS = $(CFLAGS_OPT)
CFLAGS = $(CFLAGS_OPT)

# for gcov profiling add:
#-fprofile-arcs -ftest-coverage


PROGRAMS = TestMatrix MagicSquareExample

LIBRARY=libublasJama.a

all: $(LIBRARY) $(PROGRAMS)

lib: $(LIBRARY)

ublasJama_SOURCES_CPP = \
	CholeskyDecomposition.cpp \
	EigenvalueDecomposition.cpp \
	LUDecomposition.cpp \
	QRDecomposition.cpp \
	SingularValueDecomposition.cpp

TestMatrix_SOURCES_CPP = \
	test/TestMatrix.cpp

MagicSquareExample_SOURCES_CPP = \
	examples/MagicSquareExample.cpp

ublasJama_SOURCES_C = \

ublasJama_HEADERS = \
	CholeskyDecomposition.hpp \
	EigenvalueDecomposition.hpp \
	LUDecomposition.hpp \
	QRDecomposition.hpp \
	SingularValueDecomposition.hpp

ublasJama_LIBS = $(LIBS)

ublasJama_OBJS =  $(ublasJama_SOURCES_CPP:.cpp=.o) $(ublasJama_SOURCES_C:.c=.o)
TestMatrix_OBJS = $(TestMatrix_SOURCES_CPP:.cpp=.o)
MagicSquareExample_OBJS = $(MagicSquareExample_SOURCES_CPP:.cpp=.o)

SRCS_CPP = \
	$(ublasJama_SOURCES_CPP) \
	$(TestMatrix_SOURCES_CPP) \
	$(MagicSquareExample_SOURCES_CPP)

TestMatrix:  $(TestMatrix_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(surf_LIBS) $(LDADD)

MagicSquareExample:  $(MagicSquareExample_OBJS) $(LIBRARY)
	$(LD) -o $@ $^ $(LDFLAGS) $(surf_LIBS) $(LDADD)

$(LIBRARY): $(ublasJama_OBJS)
	ar rvu $@ $^
	ranlib $@

.SUFFIXES: .c .o .cpp

## gcc-only version:
%.o : %.c
	$(COMPILE.c) -MD -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

%.o : %.cpp
	$(COMPILE.cpp) -MD -o $@ $<
	@cp $*.d $*.P; \
	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	        -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
	    rm -f $*.d

## general version:
# %.o : %.c
# 	@$(MAKEDEPEND); \
# 	    cp $*.d $*.P; \
# 	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
# 		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
# 	    rm -f $*.d
# 	$(COMPILE.c) -o $@ $<

# %.o : %.cpp
# 	@$(MAKEDEPEND); \
# 	    cp $*.d $*.P; \
# 	    sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
# 		-e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $*.P; \
# 	    rm -f $*.d
# 	$(COMPILE.cpp) -o $@ $<

clean:
	-rm -f $(PROGRAMS) $(LIBRARY) *.o */*.o *~ */*~

distclean: clean
	-rm -f $(SRCS_CPP:.cpp=.P) $(SRCS_C:.c=.P) ublasJama.xcodeproj/*.pbxuser ublasJama.xcodeproj/*.perspectivev3
	-rm -rf build Debug Release

-include $(SRCS_CPP:.cpp=.P) $(SRCS_C:.c=.P)
