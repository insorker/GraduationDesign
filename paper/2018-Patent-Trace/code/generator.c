#include "generator.h"
#include "util.h"
#include "stdio.h"
#include "assert.h"
#define ERROR_RATE_CHECK 0

/* deprecated, use gen_trc instead */
ErrorRate gen_er() {
  FILE *fp = fopen("config", "r");
  ErrorRate er = { 0, 0, 0 };
  if (fp == NULL) {
    printf("No input config\n");
    return er;
  }

  fscanf(fp, "%*s%lf", &er.sub);
  fscanf(fp, "%*s%lf", &er.del);
  fscanf(fp, "%*s%lf", &er.ins);
  er.del = er.sub + er.del;
  er.ins = er.del + er.ins;
  fclose(fp);

#if ERROR_RATE_CHECK
  printf("sub rate: %f\n", er.sub);
  printf("del rate: %f\n", er.del);
  printf("ins rate: %f\n", er.ins);
#endif

  return er;
}

TRConfig gen_trc() {
  FILE *fp = fopen("config", "r");
  TRConfig trc = { 0, 0, 0, { 0, 0, 0 } };

  if (fp == NULL) {
    printf("No input config\n");
    assert(0);
  }

  fscanf(fp, "%*s%lf", &trc.er.sub);
  fscanf(fp, "%*s%lf", &trc.er.del);
  fscanf(fp, "%*s%lf", &trc.er.ins);
  fscanf(fp, "%*s%d", &trc.num);
  fscanf(fp, "%*s%d", &trc.len);
  fscanf(fp, "%*s%d", &trc.width);
  trc.er.del = trc.er.sub + trc.er.del;
  trc.er.ins = trc.er.del + trc.er.ins;
  fclose(fp);

#if ERROR_RATE_CHECK
  printf("sub rate: %f\n", trc.er.sub);
  printf("del rate: %f\n", trc.er.del);
  printf("ins rate: %f\n", trc.er.ins);
#endif

  return trc;
}

Strand *gen_strand(int len) {
  Strand *strand = strand_new(STRAND);
  for (int i = 0; i < len; i++) {
    strand_append(strand, nrand(4));
  }
  return strand;
}

Strand *gen_read(Strand *strand, ErrorRate er) {
  Strand *read = strand_new(READ);

  for (int i = 0; i < strand->len; ) {
    double x = roll();

    if (x < er.sub) {
      strand_append(read, nrand(4));
      i++;
    }
    else if (x < er.del) {
      i++;
    }
    else if (x < er.ins) {
      strand_append(read, nrand(4));
    }
    else {
      strand_append(read, strand->seq[i]);
      i++;
    }
  }

  return read;
}
