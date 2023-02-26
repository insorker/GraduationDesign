#include <stdio.h>
#include <assert.h>
#include "../generator.h"
#include "../strand.h"
#include "../util.h"

void test_strand();
void test_generator();

int main() {
  roll_seed();

  test_strand();
  test_generator();
}

void test_strand() {
#define N "new: "
#define R "resize: "
#define L "len: "
#define S "size: "

  Strand *s = strand_new(STRAND);
  printf(N L "%d\n", s->len);

  if (strand_resize(s, 100)) {
    printf("error\n");
    assert(0);
  }
  printf(R L "%d\n", s->len);
  printf(R S "%d\n", s->size);
  if (strand_resize(s, 10)) {
    printf("error\n");
    assert(0);
  }
  printf(R L "%d\n", s->len);
  printf(R S "%d\n", s->size);
}

void test_generator() {
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
