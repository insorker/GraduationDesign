#include "strand.h"
#include <malloc.h>
#include <assert.h>

int         strand_size(Strand *);
Nucleotide  strand_at(Strand *, int index);
void        strand_push_back(Strand *, Nucleotide n);
void        strand_pop_back(Strand *);
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
  s->_nucleotides = new_vector(sizeof(Nucleotide));

  return s;
}

void free_strand(Strand *s) {
  free_vector(s->_nucleotides);
  free(s);
}

bool compare_strand(Strand *sa, Strand *sb) {
  if (sa->size(sa) != sb->size(sb)) {
    return false;
  }

  for (int i = 0; i < sa->size(sa); i++) {
    if (sa->at(sa, i) != sb->at(sb, i)) {
      return false;
    }
  }

  return true;
}

void print_strand(Strand *s) {
  printf("strand:\n");

  printf("  ");
  for (int i = 0; i < s->size(s); i++) {
    print_nucleotide(s->at(s, i));
  }
  printf("\n");
}


int strand_size(Strand *s) {
  return s->_nucleotides->size(s->_nucleotides);
}

Nucleotide strand_at(Strand *s, int index) {
  return *(Nucleotide *)s->_nucleotides->at(s->_nucleotides, index);
}

void strand_push_back(Strand *s, Nucleotide n) {
  s->_nucleotides->push_back(s->_nucleotides, &n);
}

void strand_pop_back(Strand *s) {
  s->_nucleotides->pop_back(s->_nucleotides);
}

void strand_clear(Strand *s) {
  s->_nucleotides->clear(s->_nucleotides);
}
