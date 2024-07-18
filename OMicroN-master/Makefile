appname := OMicroN
# choose compiler
CXX := g++
# compiler flags
CXXFLAGS := -std=c++17 -mno-avx -ffast-math -Wno-unknown-pragmas -march=native -fopenmp -fPIC -fmax-include-depth=300
BUILD_DIR := build
BIN_DIR := bin
SRC_DIR := src
UTILS_DIR := $(SRC_DIR)/utils
LIBS := -lm -ltiff -lpthread -fopenmp

# Definitions
DEFINES := -D_GL_LINUX_
# Includes source directories
INCLUDE := -I$(SRC_DIR) -I$(UTILS_DIR)

# Include directories for HDF5
INCLUDES := -I/usr/include/hdf5/serial

# Libraries for HDF5
LIBS += -lrt -lhdf5_serial -lhdf5_cpp -ldl -L/usr/lib/x86_64-linux-gnu/hdf5/serial
DEFINES += -DUSE_PTHREADS
SWIG_DEFINES :=

# Source files
SOURCES := $(wildcard $(SRC_DIR)/*.cpp) $(wildcard $(UTILS_DIR)/*.cpp)

# Temporary sources without the directory prefix
SOURCES_TMP := $(patsubst $(SRC_DIR)/%,%,$(SOURCES))
SOURCES_TMP := $(patsubst $(UTILS_DIR)/%,%,$(SOURCES_TMP))

# Object files
OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(SOURCES_TMP))

OMicroN: $(BIN_DIR)/OMicroN

# Main executable
$(BIN_DIR)/OMicroN: $(OBJECTS) $(BUILD_DIR)/main.o | $(BIN_DIR)
	$(CXX) -o $@ $^ $(LIBS)

## Targets for the object files

# .cpp files in the source directory
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(DEFINES) $(INCLUDE) $(INCLUDES) -o $@ $<

# .cpp files in the utils directory
$(BUILD_DIR)/%.o: $(UTILS_DIR)/%.cpp | $(BUILD_DIR)
	@mkdir -p $(dir $@)
	$(CXX) -c $(CXXFLAGS) $(DEFINES) $(INCLUDE) $(INCLUDES) -o $@ $<

## Target for making directories
$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

depend: .depend

.depend: $(SOURCES)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) $(DEFINES) $(INCLUDE) $(INCLUDES) -MM $^ >>./.depend;

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)

dist-clean: clean
	rm -f *~ .depend

include .depend

