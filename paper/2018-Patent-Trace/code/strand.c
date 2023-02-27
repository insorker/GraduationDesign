#include "strand.h"
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

Strand *strand_new(StrandType type) {
  Strand *s = (Strand *)malloc(sizeof(Strand));
  s->len = 0;
  s->type = type;
  s->size = STRAND_INIT_SIZE;
  s->seq = (Nucleotide *)malloc(sizeof(Nucleotide) * s->size);
  return s;
}

void strand_free(Strand *s) {
  free(s->seq);
  free(s);
}

Strand *strand_copy(Strand *s) {
  Strand *sp = (Strand *)malloc(sizeof(Strand));
  sp->len = s->len;
  sp->type = s->type;
  sp->size = s->size;
  sp->seq = (Nucleotide *)malloc(sizeof(Nucleotide) * s->size);
  for (int i = 0; i < s->len; i++) {
    sp->seq[i] = s->seq[i];
  }
  return sp;
}

int strand_resize(Strand *s, int len) {
  if (s == NULL   ||
      len < 0     ||
      len * 2 < 0)
  {
    return 1;
  }

  while (len < s->size / 2 && len >= STRAND_INIT_SIZE) {
    if (strand_shrink(s)) { assert(0); }
  }
  while (len >= s->size) {
    if (strand_expand(s)) { assert(0); }
  }

  s->len = len;
  return 0;
}

int strand_append(Strand *s, Nucleotide n) {
  assert(0 <= n && n <= 3);
  s->seq[s->len++] = n;
  if (s->len == s->size) {
    if (strand_resize(s, s->len)) {
      return 1;
    }
  }
  return 0;
}

void strand_print(Strand *s) {
  switch (s->type) {
    case STRAND: printf("strand: "); break;
    case READ:   printf("read:   "); break;
  }
  printf("%3d ", s->len);
  for (int i = 0; i < s->len; i++) {
    switch (s->seq[i]) {
      case A: printf("A"); break;
      case C: printf("C"); break;
      case G: printf("G"); break;
      case T: printf("T"); break;
      default: printf("#"); break;
    }
  }
  printf("\n");
}

int strand_expand(Strand *s) {
  int newSize = s->size * 2;
  if (newSize <= 0) {
    return 1;
  }
  Nucleotide *newSeq = (Nucleotide *)malloc(sizeof(Nucleotide) * newSize);

  for (int i = 0; i < s->len; i++) {
    newSeq[i] = s->seq[i];
  }

  free(s->seq);
  s->seq = newSeq;
  s->size = newSize;
  return 0;
}

int strand_shrink(Strand *s) {
  int newSize = s->size / 2;
  if (newSize < STRAND_INIT_SIZE) {
    return 1;
  }
  Nucleotide *newSeq = (Nucleotide *)malloc(sizeof(Nucleotide) * newSize);

  for (int i = 0; i < s->len; i++) {
    newSeq[i] = s->seq[i];
  }

  free(s->seq);
  s->seq = newSeq;
  s->size = newSize;
  return 0;
}
