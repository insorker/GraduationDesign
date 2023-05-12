#include "nucleotide.h"
#include <stdio.h>
#include <stdlib.h>

Nucleotide rand_nucleotide() {
  return rand() % (NUCLEOTIDE_SIZE - 1);
}

void print_nucleotide(Nucleotide n) {
  switch (n) {
    case NUCLEOTIDE_A:  printf("A"); break;
    case NUCLEOTIDE_C:  printf("C"); break;
    case NUCLEOTIDE_G:  printf("G"); break;
    case NUCLEOTIDE_T:  printf("T"); break;
    default:            printf("N"); break;
  }
}
