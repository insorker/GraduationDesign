#include "strand.h"
#include <malloc.h>
#include <stdio.h>
#include <assert.h>

Strand *strand_new(StrandType type) {
  Strand *s = (Strand *)malloc(sizeof(Strand));
  s->len = 0;
  s->type = type;
  s->size = STRAND_INIT_SIZE;
  s->state = CREDIBLE;
  s->seq = (Nucleotide *)malloc(sizeof(Nucleotide) * s->size);
  return s;
}

void strand_free(Strand *s) {
  free(s->seq);
  free(s);
}

void strand_clear(Strand *s) {
  strand_resize(s, 0);
}

Strand *strand_copy(Strand *s) {
  Strand *sp = (Strand *)malloc(sizeof(Strand));
  sp->len = s->len;
  sp->type = s->type;
  sp->size = s->size;
  sp->state = s->state;
  sp->seq = (Nucleotide *)malloc(sizeof(Nucleotide) * s->size);
  for (int i = 0; i < s->len; i++) {
    sp->seq[i] = s->seq[i];
  }
  return sp;
}

int strand_resize(Strand *s, int len) {
  if (s == NULL || len < 0 || len * 2 < 0) {
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
  assert(0 <= n && n <= 4);

  s->seq[s->len++] = n;
  if (s->len == s->size) {
    if (strand_resize(s, s->len)) {
      return 1;
    }
  }

  return 0;
}

static int min(int a, int b) {
  if (a <= b)
    return a;
  return b;
}

int strand_cmp_editdistance(Strand *s1, Strand *s2) {
  int f[200][200] = {};
  int a[200] = {}, b[200] = {};
  int n = s1->len, m = s2->len;
  for (int i = 0; i < n; i++) a[i + 1] = s1->seq[i];
  for (int i = 0; i < m; i++) b[i + 1] = s2->seq[i];

  for (int i = 0; i <= n; i++) f[i][0] = i;
  for (int i = 0; i <= m; i++) f[0][i] = i;

  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      f[i][j] = min(f[i - 1][j] + 1, f[i][j - 1] + 1);
      f[i][j] = min(f[i][j], f[i - 1][j - 1] + (a[i] != b[j]));
    }
  }

  return f[n][m];
}

void strand_print(Strand *s) {
  switch (s->type) {
    case STRAND:    printf("strand:    "); break;
    case READ:      printf("read:      "); break;
    case CONSENSUS: printf("consensus: "); break;
    default:        printf("errortype: "); break;
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

int nucleotide_cmp(Nucleotide *n1, Nucleotide *n2, int l1, int l2) {
  for (int i = 0; i < l1 || i < l2; i++) {
    Nucleotide a = i < l1 ? n1[i] : N;
    Nucleotide b = i < l2 ? n2[i] : N;

    if (a != b) {
      return 1;
    }
  }
  
  return 0;
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
