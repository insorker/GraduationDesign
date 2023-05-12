#ifndef NUCLEOTIDE_H
#define NUCLEOTIDE_H

#define NUCLEOTIDE_SIZE 5
typedef enum Nucleotide {
  NUCLEOTIDE_A,
  NUCLEOTIDE_C,
  NUCLEOTIDE_G,
  NUCLEOTIDE_T,
  NUCLEOTIDE_N
} Nucleotide;

Nucleotide rand_nucleotide();
void print_nucleotide(Nucleotide n);

#endif
