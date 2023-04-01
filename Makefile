# Description: Makefile for the project

# Compiler and flags
CXX := g++
CC := gcc
CXXFLAGS := -Wall -Wextra -pedantic -std=c++17
CFLAGS := -Wall -Wextra -pedantic -std=c11

# Libraries
LIBS := -lm

# Folders
SRC := src
INCLUDE := include
BIN := bin

# Executable name
TARGET := main

# We create a list of all the sources by looking for all the .cpp and .c files
SOURCES := $(wildcard $(SRC)/*.c) $(wildcard $(SRC)/*.cpp)

# We create a list of object files by replacing the .cpp or .c extension with .o in the list of sources
OBJECTS := $(patsubst $(SRC)/%.cpp, $(BIN)/%.o, $(filter %.cpp, $(SOURCES))) $(patsubst $(SRC)/%.c, $(BIN)/%.o, $(filter %.c, $(SOURCES))) 

# We need to tell the compiler where to find the headers
HEADERS := $(wildcard $(INCLUDE)/*.h)

#  .PHONY target specifies that all and clean are not real files, but are just targets that don't produce output files.
.PHONY: all clean

all: $(BIN)/$(TARGET)

# We link all the object files together to create the executable
$(BIN)/$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LIBS)

# We compile the .cpp files
$(BIN)/%.o: $(SRC)/%.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -c $< -o $@ $(LIBS)

# We compile the .c files
$(BIN)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) -I$(INCLUDE) -c $< -o $@ $(LIBS)

clean:
	rm -f $(BIN)/*.o $(BIN)/$(TARGET)
