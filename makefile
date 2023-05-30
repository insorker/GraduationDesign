SRC_DIRS:=vector strand sequencer
INC_DIRS:=-I./vector -I./strand -I./sequencer

CC:=gcc
CFLAGS:=${INC_DIRS} -g -fsanitize=address -Wall -Werror
SOURCES:=$(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.c))

SRC_DIRS_2018_Patent:=2018-Patent
SRC_DIRS_2020_Biorxiv:=2020_Biorxiv
SRC_DIRS_2020_Patent:=2020-Patent
INC_DIRS_2018_Patent:=-I./2018-Patent
INC_DIRS_2020_Biorxiv:=-I./2020_Biorxiv
INC_DIRS_2020_Patent:=-I./2020-Patent
SOURCES_2018_Patent:=$(foreach dir, $(SRC_DIRS_2018_Patent), $(wildcard $(dir)/*.c))
SOURCES_2020_Biorxiv:=$(foreach dir, $(SRC_DIRS_2020_Biorxiv), $(wildcard $(dir)/*.c))
SOURCES_2020_Patent:=$(foreach dir, $(SRC_DIRS_2020_Patent), $(wildcard $(dir)/*.c))

all: test

test: test_vector test_strand test_sequencer test_2018

test_vector:
	${CC} ${CFLAGS} test/test_vector.c ${SOURCES}

test_strand:
	${CC} ${CFLAGS} test/test_strand.c ${SOURCES}

test_sequencer:
	${CC} ${CFLAGS} test/test_sequencer.c ${SOURCES}

test_2018:
	${CC} ${CFLAGS} test/test_2018_Patent.c ${SOURCES} ${INC_DIRS_2018_Patent} ${SOURCES_2018_Patent}

test_2020:
	${CC} ${CFLAGS} test/test_2020_Patent.c ${SOURCES} ${INC_DIRS_2020_Patent} ${SOURCES_2020_Patent}

testbench: edit_distance

edit_distance_18:
	${CC} ${CFLAGS} testbench/edit_distance.c ${SOURCES} ${INC_DIRS_2018_Patent} ${SOURCES_2018_Patent}

edit_distance_20:
	${CC} ${CFLAGS} testbench/edit_distance.c ${SOURCES} ${INC_DIRS_2020_Patent} ${SOURCES_2020_Patent}

clean:
	rm -f *.out
