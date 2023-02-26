#include "generator.h"
#include "util.h"
#include "stdio.h"
#include "assert.h"

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
