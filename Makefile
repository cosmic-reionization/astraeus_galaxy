PYTHON = python3.11
CXX = g++

PYTHONINC := $(shell $(PYTHON)-config --includes)
NUMPYINC := $(shell $(PYTHON) -c "import numpy; print(numpy.get_include())")
EXT_SUFFIX := $(shell $(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX') or '.so')")

CPPFLAGS += -I$(NUMPYINC) $(PYTHONINC)
CXXFLAGS += -O3 -std=c++17 -fPIC
LDFLAGS += -shared

TARGET := astraeus$(EXT_SUFFIX)
SOURCES := python_module.cpp evol_gal.cpp radiative_feedback.cpp starformation_SNfeedback.cpp
OBJECTS := $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(TARGET)
