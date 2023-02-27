#include "generator.h"
#include "util.h"
#include "stdio.h"
#include "assert.h"
#define ERROR_RATE_CHECK 0

ErrorRate gen_er() {
  FILE *fp = fopen("config", "r");
  ErrorRate er = { 0, 0, 0 };
  if (fp == NULL) {
    assert(0);
  }

  fscanf(fp, "%*s%lf", &er.sub);
  fscanf(fp, "%*s%lf", &er.ins);
  fscanf(fp, "%*s%lf", &er.del);
  er.ins = er.sub + er.ins;
  er.del = er.ins + er.del;
  fclose(fp);

#if ERROR_RATE_CHECK
  printf("sub rate: %f\n", er.sub);
  printf("ins rate: %f\n", er.ins);
  printf("del rate: %f\n", er.del);
#endif

  return er;
}

Strand *gen_strand_rand(int len) {
  Strand *strand = strand_new(STRAND);
  for (int i = 0; i < len; i++) {
    strand_append(strand, nrand(4));
  }
  return strand;
}

Strand *gen_strand_str(char *str) {
  return NULL;
}

Strand *gen_strand_file(char *file) {
  return NULL;
}

Strand *gen_read(Strand *strand, ErrorRate er) {
  Strand *read = strand_new(READ);
  for (int i = 0; i < strand->len; ) {
    double x = roll();
#if ERROR_RATE_CHECK
    printf("read error rate: %f\n", x);
#endif

    if (x < er.sub) {
      strand_append(read, nrand(4));
      i++;
    }
    else if (x < er.ins) {
      strand_append(read, nrand(4));
    }
    else if (x < er.del) {
      i++;
    }
    else {
      strand_append(read, strand->seq[i]);
      i++;
    }
  }

  return read;
}
