SRC_DIRS = strand sequencer tr_system
INC_DIRS = $(patsubst %, -I./%, ${SRC_DIRS}) -I./KiteSTL/inc

CC 			= gcc
CFLAGS 	= ${INC_DIRS} -g -fsanitize=address -Wall -Werror
LDFLAGS = -Wl,-rpath,./KiteSTL/lib -L./KiteSTL/lib
LDLIBS  = -lKiteSTL
SOURCES = $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.c))

SOURCES_2018_Patent  = $(foreach dir, $(2018_Patent), $(wildcard $(dir)/*.c))
SOURCES_2020_Patent	 = $(foreach dir, $(2020_Patent), $(wildcard $(dir)/*.c))

.PHONY: all test clean

all: test

test: test_vector test_strand test_sequencer test_tr_system_normal

test_vector:
	${CC} ${CFLAGS} test/test_vector.c ${SOURCES} ${LDFLAGS} ${LDLIBS}

test_strand:
	${CC} ${CFLAGS} test/test_strand.c ${SOURCES} ${LDFLAGS} ${LDLIBS}

test_sequencer:
	${CC} ${CFLAGS} test/test_sequencer.c ${SOURCES} ${LDFLAGS} ${LDLIBS}

test_tr_system_normal:
	${CC} ${CFLAGS} test/test_tr_system.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=1

test_tr_system_indeterminant:
	${CC} ${CFLAGS} test/test_tr_system.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=0 -DTR_INDETERMINATED=1 -DTR_ALIGNMENT=0

test_tr_system_alignment:
	${CC} ${CFLAGS} test/test_tr_system.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=0 -DTR_INDETERMINATED=0 -DTR_ALIGNMENT=1

test_tr_system_indeterminant_alignment:
	${CC} ${CFLAGS} test/test_tr_system.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=0 -DTR_INDETERMINATED=1 -DTR_ALIGNMENT=1

test_edit_distance_3:
	${CC} ${CFLAGS} test/test_edit_distance_3.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=1

test_edit_distance_41:
	${CC} ${CFLAGS} test/test_edit_distance_4.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=0 -DTR_INDETERMINATED=1 -DTR_ALIGNMENT=0

test_edit_distance_42:
	${CC} ${CFLAGS} test/test_edit_distance_4.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=0 -DTR_INDETERMINATED=0 -DTR_ALIGNMENT=1

test_edit_distance_43:
	${CC} ${CFLAGS} test/test_edit_distance_4.c ${SOURCES} ${LDFLAGS} ${LDLIBS} -DTR_NORMAL=0 -DTR_INDETERMINATED=1 -DTR_ALIGNMENT=1

clean:
	rm -f *.out
