#include "strand.h"
#include <malloc.h>
#include <assert.h>
#include <limits.h>
#include <string.h>

#define STRAND_MIN_SIZE 3
#define STRAND_MAX_SIZE INT_MAX
#define STRAND_MAX_LENGTH (STRAND_MAX_SIZE / 4) // 4 是大于 3 的最大二次幂

int         strand_size(Strand *);
Nucleotide  strand_at(Strand *, int index);
void        strand_push_back(Strand *, Nucleotide n);
Nucleotide  strand_pop_back(Strand *);
void        strand_clear(Strand *);

void strand_capacity(Strand *, int newsize);
void strand_expand(Strand *);
void strand_shrink(Strand *);

Strand *new_strand() {
  Strand *s = (Strand *)malloc(sizeof(Strand));

/* public */
  s->size = strand_size;
  s->at = strand_at;
  s->push_back = strand_push_back;
  s->pop_back = strand_pop_back;
  s->clear = strand_clear;

/* private */
  s->length = 0;
  s->_size = STRAND_MIN_SIZE;
  s->items = (Nucleotide *)malloc(s->_size * sizeof(Nucleotide));

  s->capacity = strand_capacity;
  s->expand = strand_expand;
  s->shrink = strand_shrink;

  return s;
}

void free_strand(Strand *s) {
  free(s->items);
  free(s);
}

int compare_strand(Strand *sa, Strand *sb) {
  if (sa->length != sb->length) {
    return 1;
  }

  for (int i = 0; i < sa->length; i++) {
    if (sa->at(sa, i) != sb->at(sb, i)) {
      return 1;
    }
  }
  
  return 0;
}

void print_strand(Strand *s) {
  printf("strand: ");
  for (int i = 0; i < s->length; i++) {
    print_nucleotide(s->items[i]);
  }
  printf("\n");
}

int strand_size(Strand *s) {
  return s->length;
}

Nucleotide strand_at(Strand *s, int index) {
  assert(index >= 0);
  assert(index < s->length);

  return s->items[index];
}

void strand_push_back(Strand *s, Nucleotide n) {
  s->capacity(s, s->length + 1);
  s->length += 1;
  s->items[s->length - 1] = n;
}

Nucleotide strand_pop_back(Strand *s) {
  Nucleotide item = s->items[s->length - 1];

  s->capacity(s, s->length - 1);
  s->length -= 1;

  return item;
}

void strand_clear(Strand *s) {
  while (s->length != 0) {
    s->pop_back(s);
  }
}

void strand_capacity(Strand *s, int nlength) {
  assert(nlength >= 0);
  assert(nlength <= STRAND_MAX_LENGTH);

  while (nlength >= s->_size)
    s->expand(s);
  while (nlength < s->_size / 2 && s->_size > STRAND_MIN_SIZE)
    s->shrink(s);
}

void strand_expand(Strand *s) {
  int nsize = s->_size * 2;
  Nucleotide *nitems = (Nucleotide *)malloc(nsize * sizeof(Nucleotide));
  
  for (int i = 0; i < s->_size; i++) {
    nitems[i] = s->items[i];
  }

  free(s->items);
  s->items = nitems;
  s->_size = nsize;
}

void strand_shrink(Strand *s) {
  int nsize = s->_size / 2;
  Nucleotide *nitems = (Nucleotide *)malloc(nsize * sizeof(Nucleotide));
  
  for (int i = 0; i < nsize; i++) {
    nitems[i] = s->items[i];
  }

  free(s->items);
  s->items = nitems;
  s->_size = nsize;
}
