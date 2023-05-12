#include "read.h"
#include <malloc.h>
#include <stdio.h>

int         read_size(Read *);
Nucleotide  read_at(Read *, int index);
void        read_push_back(Read *, Nucleotide n);
Nucleotide  read_pop_back(Read *);
void        read_clear(Read *);

Read *new_read() {
  Read *read = (Read *)malloc(sizeof(Read));

/* private */
  read->super = new_strand();

/* public extends */
  read->size = read_size;
  read->at = read_at;
  read->push_back = read_push_back;
  read->pop_back = read_pop_back;
  read->clear = read_clear;

  return read;
}

void free_read(Read *read) {
  free_strand(read->super);
  free(read);
}

void print_read(Read *read) {
  printf("Read: \n");
  printf("- "), print_strand(read->super);
}

int read_size(Read *read) {
  return read->super->size(read->super);
}

Nucleotide read_at(Read *read, int index) {
  return read->super->at(read->super, index);
}

void read_push_back(Read *read, Nucleotide n) {
  return read->super->push_back(read->super, n);
}

Nucleotide read_pop_back(Read *read) {
  return read->super->pop_back(read->super);
}

void read_clear(Read *read) {
  return read->super->clear(read->super);
}
