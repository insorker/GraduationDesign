#include "tr_read.h"
#include <malloc.h>
#include <stdio.h>

int         tr_read_size(TRRead *);
Nucleotide  tr_read_at(TRRead *, int index);
void        tr_read_push_back(TRRead *, Nucleotide n);
Nucleotide  tr_read_pop_back(TRRead *);
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
  tr_read->state = TR_READ_STATE_CREDIBLE;

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
}

void print_tr_read(TRRead *tr_read) {
  printf("TRRead: \n");
  printf("- "), print_read(tr_read->super);
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

Nucleotide tr_read_pop_back(TRRead *tr_read) {
  return tr_read->super->pop_back(tr_read->super);
}

void tr_read_clear(TRRead *tr_read) {
  return tr_read->super->clear(tr_read->super);
}
