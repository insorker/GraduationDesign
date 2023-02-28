#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include "generator.h"
#include "strand.h"
#include "trsystem.h"
#include "config.h"
#include "util.h"

int main() {
  roll_seed();

  TRConfig trc = gen_trc();
  Strand *raw = gen_strand(trc.len);
  Strand *reads[trc.num];

  strand_print(raw);
  for (int i = 0; i < trc.num; i++) {
    reads[i] = gen_read(raw, trc.er);
    printf("%2d ", i);
    strand_print(reads[i]);
  }

  TRSystem *trs = trs_new(trc, reads);
  trs_run(trs);
  printf("\n");
  trs_print(trs);
  printf("edit distance: %d\n", strand_cmp_editdistance(raw, trs->consensus));

  trs_free(trs);
  strand_free(raw);
  for (int i = 0; i < trc.num; i++) {
    strand_free(reads[i]);
  }
}
