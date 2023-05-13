#include "read.h"
#include <malloc.h>
#include <stdio.h>

/* public extends */
int           read_size(Read *);
Nucleotide    read_at(Read *, int index);
void          read_push_back(Read *, Nucleotide);
Nucleotide    read_pop_back(Read *);
void          read_clear(Read *);


Read *new_read() {
  Read *read = (Read *)malloc(sizeof(Read));

/* super */
  read->super = new_strand();

/* public extends */
  read->size = read_size;
  read->at = read_at;
  read->push_back = read_push_back;
  read->pop_back = read_pop_back;
  read->clear = read_clear;

/* public */
  read->error = new_vector();

  return read;
}

void free_read(Read *read) {
  free_strand(read->super);
  free_vector(read->error);
  free(read);
}

void print_read(Read *read) {
  printf("Read:\n");
  printf("- "), print_strand(read->super);
  printf("- "), print_read_error(read);
}

void print_read_error(Read *read) {
  Vector *error = read->error;

  printf("error: ");
  for (int i = 0; i < error->size(error); i++) {
    switch (error->at(error, i)) {
      case READ_ERROR_NONE: printf("n"); break;
      case READ_ERROR_SUB:  printf("s"); break;
      case READ_ERROR_DEL:  printf("d"); break;
      case READ_ERROR_INS:  printf("i"); break;
    }
  }
}


int read_size(Read *read) {
  return read->super->size(read->super);
}

Nucleotide read_at(Read *read, int index) {
  return read->super->at(read->super, index);
}

void read_push_back(Read *read, Nucleotide n) {
  read->super->push_back(read->super, n);
}

Nucleotide read_pop_back(Read *read) {
  return read->super->pop_back(read->super);
}

void read_clear(Read *read) {
  return read->super->clear(read->super);
}
