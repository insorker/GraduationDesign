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

void test_create(SequencerErrorRate ser) {
  PRINT_TEST_FUNC();

  Sequencer *seq = new_sequencer(ser);
  Strand *s = test_create_strand();

  seq->set_strand(seq, s);

  Read *read = seq->process(seq);

  printf("- "), print_strand(s);
  print_read(read);

  free_read(read);
  free_strand(s);
  free_sequencer(seq);

  PRINT_TEST_FUNC();
}

void test_anchor() {
  Sequencer *seq = new_sequencer(ser_sm);
  Strand *s = test_create_strand();
  vector_t *anchor = new_vector(sizeof(double));

  anchor->push_back(anchor, &(double){0.1});
  anchor->push_back(anchor, &(double){0.4});
  anchor->push_back(anchor, &(double){0.5});

  seq->set_strand(seq, s);
  seq->set_anchor(seq, anchor);

  Read *read = seq->process(seq);

  printf("- "), print_strand(s);
  print_read(read);

  free_read(read);
  free_strand(s);
  free_vector(anchor);
  free_sequencer(seq);
}

int main() {
  PRINT_TEST_FILE();

  srand((unsigned int)time(NULL));

  test_create(ser_null);
  test_create(ser_sm);
  test_create(ser_mid);
  test_create(ser_lg);
  test_create(ser_full);

  test_anchor();
}
