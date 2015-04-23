#
# Modify HALIDE_DIR to the appropriate path on your machine.
#
# Special instructions for Mac users
# ==================================
# You need to or must have installed libpng through Macports or Homebrew.
# Assuming that the installation succeeded, you should be able to run
#
# The brew command for installing libpng is
# brew install libpng
#
# libpng-config --I_opts
# libpng-config --L_opts
#
# Please add the output of the above commands to the following variables:
# PNG_INC
# PNG_LIB
#

MKDIR	:= mkdir -p
RM		:= rm -f
CP		:= cp -f
CXX		:= g++ -std=c++11

HALIDE_DIR ?= $(HOME)/Halide
HALIDE_LIB := $(HALIDE_DIR)/bin/libHalide.a

INC  := $(wildcard  *.h)
SRC  := $(wildcard  tutorial*.cpp)
APPS := $(notdir $(SRC:.cpp=))

PNG_INC :=
PNG_LIB :=

all: a9 $(APPS)
	mkdir -p Output

CXXFLAGS := $(PNG_INC) -I$(HALIDE_DIR)/include/ -I. -g -Wall
LDFLAGS  := $(PNG_LIB) -L$(HALIDE_DIR)/bin/     -lz -lpthread -ldl -lncurses -lpng

%: %.cpp $(HALIDE_LIB)
	$(CXX) $(CXXFLAGS) $< $(HALIDE_LIB) $(LDFLAGS) -o $@

a9: a9_main.cpp a9_pset.cpp a9.h $(HALIDE_LIB)
	$(CXX) $(CXXFLAGS) a9_main.cpp a9_pset.cpp $(HALIDE_LIB) $(LDFLAGS) -o $@

run:
	./a9

clean:
	$(RM) $(APPS) Output/*
	rm -rf Output
