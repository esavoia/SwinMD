#=======================================================
# Makefile for SwinMD C++ MPI(SGI MPT)  program        =
# Created by Jianhui LI, CMS, Swinburne Univ. of Tech. =
# @Copyright Jianhui LI, jli@it.swin.edu.au            =
# Usage: make                                          =
#=======================================================

## DIRECTORIES

SRC_BASE_DIR := ./src

OBJS_DIR := objs

PROD_DIR := products

# C++ files
SRCS := $(wildcard $(SRC_BASE_DIR)/*.cpp)
#$(shell find . -name *.cpp)

# Object files
OBJS := $(SRCS:$(SRC_BASE_DIR)/%.cpp=$(OBJS_DIR)/%.o)

# Lib Dependencies files
DEPS := $(OBJS:.o=.d)

# Default target name
TARGET := $(PROD_DIR)/SwinMD

RUN_DIR := $(PROD_DIR)/_run


# Forces
F_SRCS := $(wildcard $(SRC_BASE_DIR)/forces/*.cpp)
F_OBJS := $(F_SRCS:$(SRC_BASE_DIR)/forces/%.cpp=$(OBJS_DIR)/%.o)
# Integrators
I_SRCS := $(wildcard $(SRC_BASE_DIR)/integrators/*.cpp)
I_OBJS := $(I_SRCS:$(SRC_BASE_DIR)/integrators/%.cpp=$(OBJS_DIR)/%.o)
# Observables
O_SRCS := $(wildcard $(SRC_BASE_DIR)/observables/*.cpp)
O_OBJS := $(O_SRCS:$(SRC_BASE_DIR)/observables/%.cpp=$(OBJS_DIR)/%.o)
# Utils
U_SRCS := $(wildcard $(SRC_BASE_DIR)/utils/*.cpp)
U_OBJS := $(U_SRCS:$(SRC_BASE_DIR)/utils/%.cpp=$(OBJS_DIR)/%.o)

## COMPILATION STUFFS

#CC ?= mpicxx
CC = mpicxx
## FLAGS

# CPPFLAGS: C PreProcessor Flags
# The -MMD and -MP flags together generate Makefiles for us!
# These files will have .d instead of .o as the output.
CPPFLAGS := $(INC_FLAGS) -MMD -MP

# Warnings flags
WARNINGS := -Wall -Wextra

MORE_WARNINGS := -Wshadow -Wpointer-arith -Wcast-align -Wwrite-strings -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls -Wnested-externs -Winline -Wno-long-long -Wuninitialized -Wstrict-prototypes

# C++ flags
CFLAGS ?= -O2 -g $(WARNINGS)
#CFLAGS ?= -O2 -g3 $(WARNINGS) $(MORE_WARNINGS)
#CFLAGS := $(CFLAGS) -std=c99 -U__STRICT_ANSI__


$(TARGET): $(OBJS_DIR) $(PROD_DIR) $(OBJS) $(U_OBJS) $(F_OBJS) $(I_OBJS) $(O_OBJS)
	@echo $(OBJS)
	@echo $(U_OBJS)
	@echo "\x1B[31mMake SwinMD\x1B[0m"
	$(CC) -o $(TARGET) $(CFLAGS) $(OBJS) $(U_OBJS) $(F_OBJS) $(I_OBJS) $(O_OBJS)
	@echo "Make SwinMD \x1B[31mDone\x1B[0m"

$(OBJS_DIR):
	@echo "\x1B[31m Build <objs> dir\x1B[0m: \"$@\""
	mkdir -p $@

$(PROD_DIR):
	@echo "\x1B[31m Build <products> dir\x1B[0m: \"$@\""
	mkdir -p $@

$(OBJS): $(OBJS_DIR)/%.o : $(SRC_BASE_DIR)/%.cpp
	@echo "\x1B[31m Make obj \x1B[0m: \"$@\" "
	@echo " from src: \"$<\""
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(U_OBJS): $(OBJS_DIR)/%.o : $(SRC_BASE_DIR)/utils/%.cpp
	@echo "\x1B[31m Make Util Obj\x1B[0m: \"$@\" "
	@echo " from src: \"$<\""
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(F_OBJS): $(OBJS_DIR)/%.o : $(SRC_BASE_DIR)/forces/%.cpp
	@echo "\x1B[31m Make Force Obj\x1B[0m: \"$@\" "
	@echo " from src: \"$<\""
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(I_OBJS): $(OBJS_DIR)/%.o : $(SRC_BASE_DIR)/integrators/%.cpp
	@echo "\x1B[31m Make Integrator Obj\x1B[0m: \"$@\" "
	@echo " from src: \"$<\""
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(O_OBJS): $(OBJS_DIR)/%.o : $(SRC_BASE_DIR)/observables/%.cpp
	@echo "\x1B[31m Make Observable Obj\x1B[0m: \"$@\" "
	@echo " from src: \"$<\""
	mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

all: $(TARGET)
	@echo "Make all:\x1B[31m DONE \x1B[0m"

clean:
	@echo "\x1B[31mMake clean\x1B[0m"
	$(RM) -r $(OBJS_DIR)
	@echo "Make clean:\x1B[31m DONE \x1B[0m"

cleanAll:
	@echo "\x1B[31mMake cleanAll\x1B[0m"
	$(RM) -r $(OBJS_DIR) $(PROD_DIR)
	@echo "Make clean all:\x1B[31m DONE \x1B[0m"

run: $(TARGET)
	mkdir -p $(RUN_DIR)
	cp $(TARGET) $(RUN_DIR)/
#	mpirun -n 4 ...
	@echo "Make run:\x1B[31m DONE \x1B[0m"

cleanRun:
	@echo "\x1B[31mMake cleanRun\x1B[0m"
	$(RM) -r $(RUN_DIR)/*.out
	@echo "Make clean run:\x1B[31m DONE \x1B[0m"

-include $(DEPS)
