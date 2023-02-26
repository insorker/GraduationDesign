#include <stdlib.h>
#include <time.h>

void roll_seed() {
  srand((unsigned int)time(NULL));
}

double roll() {
  return (double)(rand() % RAND_MAX) / RAND_MAX;
}

int nrand(int n) {
  return rand() % n;
}
