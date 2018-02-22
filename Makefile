DEP := .d
BUILD := build
BIN := bin
DATA := data
INC := -Iinclude -Isolver

CC := c++
CFLAGS := -Ofast -Wall -std=c++14
DFLAGS := -MMD -MP -MF $(DEP)

SRCS := $(wildcard tester/*.cpp)
TARGET := $(patsubst %.cpp, $(BIN)/%, $(notdir $(SRCS)))

$(shell mkdir -p $(DEP))
$(shell mkdir -p $(DATA))

all: $(TARGET)

$(BIN)/%: tester/%.cpp
	@mkdir -p $(BIN)
	$(CC) $(CFLAGS) $(INC) -o $@ $< $(DFLAGS)/$*.d

clean:
	@echo " Cleaning..."
	$(RM) $(BIN)/* $(DEP)/*

.PHONY: all clean
-include $(patsubst %.cpp, $(DEP)/%.d, $(notdir $(SRCS)))