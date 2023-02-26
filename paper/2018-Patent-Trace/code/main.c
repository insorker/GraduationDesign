#include <stdio.h>
#include <assert.h>
#include "generator.h"
#include "strand.h"
#include "util.h"

int main() {
  roll_seed();

  ErrorRate er = gen_er();
  Strand *raw = gen_strand_rand(40);
  strand_print(raw);

  Strand *read;
  for (int i = 0; i < 10; i++) {
    read = gen_read(raw, er);
    strand_print(read);
    strand_free(read);
  }

  strand_free(raw);
}
