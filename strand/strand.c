#include "strand.h"
#include <malloc.h>
#include <assert.h>

int         strand_size(Strand *);
Nucleotide  strand_at(Strand *, int index);
void        strand_push_back(Strand *, Nucleotide n);
Nucleotide  strand_pop_back(Strand *);
void        strand_clear(Strand *);


Strand *new_strand() {
  Strand *s = (Strand *)malloc(sizeof(Strand));

/* public */
  s->size = strand_size;
  s->at = strand_at;
  s->push_back = strand_push_back;
  s->pop_back = strand_pop_back;
  s->clear = strand_clear;

/* private */
  s->_nucleotide = new_vector();

  return s;
}

void free_strand(Strand *s) {
  free_vector(s->_nucleotide);
  free(s);
}

int compare_strand(Strand *sa, Strand *sb) {
  return compare_vector(sa->_nucleotide, sb->_nucleotide);
}

void print_strand(Strand *s) {
  printf("strand: ");
  for (int i = 0; i < s->size(s); i++) {
    print_nucleotide(s->at(s, i));
  }
  printf("\n");
}


int strand_size(Strand *s) {
  return s->_nucleotide->size(s->_nucleotide);
}

Nucleotide strand_at(Strand *s, int index) {
  return s->_nucleotide->at(s->_nucleotide, index);
}

void strand_push_back(Strand *s, Nucleotide n) {
  s->_nucleotide->push_back(s->_nucleotide, n);
}

Nucleotide strand_pop_back(Strand *s) {
  return s->_nucleotide->pop_back(s->_nucleotide);
}

void strand_clear(Strand *s) {
  s->_nucleotide->clear(s->_nucleotide);
}
