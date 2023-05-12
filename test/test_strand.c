#include <stdio.h>
#include "test.h"
#include "strand.h"
#include "nucleotide.h"

void test_strand_push_pop() {
  PRINT_TEST_FUNC();

  Strand *s = new_strand();

  s->push_back(s, NUCLEOTIDE_A);
  s->push_back(s, NUCLEOTIDE_C);
  s->push_back(s, NUCLEOTIDE_C);
  s->push_back(s, NUCLEOTIDE_C);
  s->push_back(s, NUCLEOTIDE_C);
  s->push_back(s, NUCLEOTIDE_C);
  print_strand(s);
  s->push_back(s, NUCLEOTIDE_G);
  s->push_back(s, NUCLEOTIDE_C);
  print_strand(s);

  s->pop_back(s);
  print_strand(s);
  s->clear(s);
  print_strand(s);

  s->push_back(s, NUCLEOTIDE_T);
  print_strand(s);

  free_strand(s);

  printf("\n");
}

void test_strand_size() {
  PRINT_TEST_FUNC();

  Strand *s = new_strand();

  for (int i = 0; i < 100; i++) {
    s->push_back(s, NUCLEOTIDE_A);
  }
  print_strand(s);
  printf("%d is equal to 100\n", s->length);

  free_strand(s);

  printf("\n");
}

void test_strand_compare() {
  PRINT_TEST_FUNC();

  Strand *sa = new_strand();
  Strand *sb = new_strand();

  for (int i = 0; i < 10; i++) {
    sa->push_back(sa, NUCLEOTIDE_A);
  }
  for (int i = 0; i < 10; i++) {
    sb->push_back(sb, NUCLEOTIDE_A);
  }

  print_strand(sa);
  print_strand(sb);
  printf("%d\n", compare_strand(sa, sb));

  sb->pop_back(sb);
  print_strand(sa);
  print_strand(sb);
  printf("%d\n", compare_strand(sa, sb));

  sa->pop_back(sa);
  sa->pop_back(sa);
  sa->push_back(sa, NUCLEOTIDE_C);
  print_strand(sa);
  print_strand(sb);
  printf("%d\n", compare_strand(sa, sb));
  
  free_strand(sa);
  free_strand(sb);

  printf("\n");
}

int main() {
  PRINT_TEST_FILE();

  test_strand_push_pop();
  test_strand_size();
  test_strand_compare();
}
