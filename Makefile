# This makefile will build the C++ components of Porechop.

# Example commands:
#   make (build in release mode)
#   make debug (build in debug mode)
#   make clean (deletes *.o files, which aren't required to run the aligner)
#   make distclean (deletes *.o files and the *.so file, which is required to run the aligner)
#   make CXX=g++-5 (build with a particular compiler)
#   make CXXFLAGS="-Werror -g3" (build with particular compiler flags)


# CXX and CXXFLAGS can be overridden by the user.
CXX         ?= g++
CXXFLAGS    ?= -Wall -Wextra -pedantic -mtune=native

# These flags are required for the build to work.
FLAGS        = -std=c++14 -Iporechop/include -fPIC
LDFLAGS      = -shared

# Different debug/optimisation levels for debug/release builds.
DEBUGFLAGS   = -DSEQAN_ENABLE_DEBUG=1 -g
RELEASEFLAGS = -O3 -D NDEBUG

TARGET       = porechop/cpp_functions.so
SHELL        = /bin/sh
SOURCES      = $(shell find porechop -name "*.cpp")
HEADERS      = $(shell find porechop -name "*.h")
OBJECTS      = $(SOURCES:.cpp=.o)

# Linux needs '-soname' while Mac needs '-install_name'
PLATFORM     = $(shell uname)
ifeq ($(PLATFORM), Darwin)
SONAME       = -install_name
else
SONAME       = -soname
endif

.PHONY: release
release: FLAGS+=$(RELEASEFLAGS)
release: $(TARGET)

.PHONY: debug
debug: FLAGS+=$(DEBUGFLAGS)
debug: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(FLAGS) $(CXXFLAGS) $(LDFLAGS) -Wl,$(SONAME),$(TARGET) -o $(TARGET) $(OBJECTS)

clean:
	$(RM) $(OBJECTS)

distclean: clean
	$(RM) $(TARGET)

%.o: %.cpp $(HEADERS)
	$(CXX) $(FLAGS) $(CXXFLAGS) -c -o $@ $<
