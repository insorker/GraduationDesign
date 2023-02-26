#include <stdio.h>
#include "../strand.h"

#define N "new: "
#define R "resize: "
#define L "len: "
#define S "size: "

int main() {
  Strand *s = strand_new(STRAND);
  printf(N L "%d\n", s->len);

  if (strand_resize(&s, 100)) {
    printf("error\n");
    return 0;
  }
  printf(R L "%d\n", s->len);
  printf(R S "%d\n", s->size);
  if (strand_resize(&s, 10)) {
    printf("error\n");
    return 0;
  }
  printf(R L "%d\n", s->len);
  printf(R S "%d\n", s->size);
  // strand_free(s);
}
