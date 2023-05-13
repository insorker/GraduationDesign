#include "sequencer.h"
#include "strand.h"
#include "nucleotide.h"
#include "read.h"
#include "test.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

SequencerErrorRate ser_null = { 0, 0, 0 };
SequencerErrorRate ser_sm = { 0.05, 0.03, 0.01 };
SequencerErrorRate ser_mid = { 0.05, 0.05, 0.05 };
SequencerErrorRate ser_lg = { 0.1, 0.2, 0.35 };
SequencerErrorRate ser_full = { 0.3, 0.3, 0.4 };

Strand *test_create_strand() {
  Strand *s = new_strand();

  for (int i = 0; i < 30; i++) {
    s->push_back(s, NUCLEOTIDE_A);
  }
  
  return s;
}

void test_create() {
  PRINT_TEST_FUNC();

  Sequencer *seq = new_sequencer(ser_null);
  Strand *s = test_create_strand();
  Read *read = seq->process(seq, s);

  printf("- "), print_strand(s);
  print_read(read);

  free_read(read);
  free_strand(s);
  free_sequencer(seq);

  printf("\n");
}

void test_create_sm() {
  PRINT_TEST_FUNC();

  Sequencer *seq = new_sequencer(ser_sm);
  Strand *s = test_create_strand();
  Read *read = seq->process(seq, s);

  printf("- "), print_strand(s);
  print_read(read);

  free_read(read);
  free_strand(s);
  free_sequencer(seq);

  printf("\n");
}

void test_create_mid() {
  PRINT_TEST_FUNC();

  Sequencer *seq = new_sequencer(ser_mid);
  Strand *s = test_create_strand();
  Read *read = seq->process(seq, s);

  printf("- "), print_strand(s);
  print_read(read);

  free_read(read);
  free_strand(s);
  free_sequencer(seq);

  printf("\n");
}

void test_create_lg() {
  PRINT_TEST_FUNC();

  Sequencer *seq = new_sequencer(ser_lg);
  Strand *s = test_create_strand();
  Read *read = seq->process(seq, s);

  printf("- "), print_strand(s);
  print_read(read);

  free_read(read);
  free_strand(s);
  free_sequencer(seq);

  printf("\n");
}

void test_create_full() {
  PRINT_TEST_FUNC();

  Sequencer *seq = new_sequencer(ser_full);
  Strand *s = test_create_strand();
  Read *read = seq->process(seq, s);

  printf("- "), print_strand(s);
  print_read(read);

  free_read(read);
  free_strand(s);
  free_sequencer(seq);

  printf("\n");
}

int main() {
  PRINT_TEST_FILE();

  srand((unsigned int)time(NULL));

  test_create();
  test_create_sm();
  test_create_mid();
  test_create_lg();
  test_create_full();
}
