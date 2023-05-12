SRC_DIRS:=strand sequencer tr_system
INC_DIRS:=-I./strand -I./sequencer -I./tr_system

CC:=gcc
CFLAGS:=${INC_DIRS} -g -fsanitize=address -Wall -Werror
SOURCES:=$(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.c))

all: test

test: test_strand test_sequencer test_tr_system

test_strand:
	${CC} ${CFLAGS} test/test_strand.c ${SOURCES}

test_sequencer:
	${CC} ${CFLAGS} test/test_sequencer.c ${SOURCES}

test_tr_system:
	${CC} ${CFLAGS} test/test_tr_system.c ${SOURCES}

clean:
	rm -f *.out
