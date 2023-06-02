#include "tr_read.h"
#include <malloc.h>
#include <stdio.h>

int         tr_read_size(TRRead *);
Nucleotide  tr_read_at(TRRead *, int index);
void        tr_read_push_back(TRRead *, Nucleotide n);
void        tr_read_pop_back(TRRead *);
void        tr_read_clear(TRRead *);


TRRead *new_tr_read() {
  TRRead *tr_read = (TRRead *)malloc(sizeof(TRRead));

/* super */
  tr_read->super = new_read();

/* public extends */
  tr_read->size = tr_read_size;
  tr_read->at = tr_read_at;
  tr_read->push_back = tr_read_push_back;
  tr_read->pop_back = tr_read_pop_back;
  tr_read->clear = tr_read_clear;

/* public */
  tr_read->state = TR_READ_STATE_ACTIVE;

  return tr_read;
}

void free_tr_read(TRRead *tr_read) {
  free_read(tr_read->super);
  free(tr_read);
}

void copy_tr_read(TRRead *tr_read, Read *read) {
  tr_read->clear(tr_read);

  for (int i = 0; i < read->size(read); i++) {
    tr_read->push_back(tr_read, read->at(read, i));
  }

  for (int j = 0; j < read->error->size(read->error); j++) {
    tr_read->super->error->push_back(
        tr_read->super->error,
        read->error->at(read->error, j));
  }

  for (int j = 0; j < read->anchor->size(read->anchor); j++) {
    tr_read->super->anchor->push_back(
        tr_read->super->anchor,
        read->anchor->at(read->anchor, j));
  }
}

void print_tr_read(TRRead *tr_read) {
  printf("TRRead:\n");

  printf("  ");
  for (int i = 0; i < tr_read->size(tr_read); i++) {
    print_nucleotide(tr_read->at(tr_read, i));
  }
  printf("\n");

  printf("  "), print_read_error(tr_read->super);

  printf("  "), print_read_anchor(tr_read->super);

  printf("  "), print_tr_read_state(tr_read);
}

void print_tr_read_state(TRRead *tr_read) {
  printf("state: ");

  switch (tr_read->state) {
    case TR_READ_STATE_ACTIVE:    printf("active"); break;
    case TR_READ_STATE_INACTIVE:  printf("inactive"); break;
    case TR_READ_STATE_VARIANT:   printf("variant"); break;
    case TR_READ_STATE_OMITTED:   printf("omitted"); break;
  }
}


int tr_read_size(TRRead *tr_read) {
  return tr_read->super->size(tr_read->super);
}

Nucleotide tr_read_at(TRRead *tr_read, int index) {
  return tr_read->super->at(tr_read->super, index);
}

void tr_read_push_back(TRRead *tr_read, Nucleotide n) {
  return tr_read->super->push_back(tr_read->super, n);
}

void tr_read_pop_back(TRRead *tr_read) {
  tr_read->super->pop_back(tr_read->super);
}

void tr_read_clear(TRRead *tr_read) {
  return tr_read->super->clear(tr_read->super);
}
