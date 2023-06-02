#include "test.h"
#include "tr_system.h"
#include "sequencer.h"
#include <stdlib.h>
#include <time.h>

SequencerErrorRate ser_null = { 0, 0, 0 };
SequencerErrorRate ser_sm = { 0.05, 0.03, 0.01 };
SequencerErrorRate ser_mid = { 0.05, 0.05, 0.05 };
SequencerErrorRate ser_lg = { 0.1, 0.2, 0.35 };
SequencerErrorRate ser_full = { 0.3, 0.3, 0.4 };

Strand *test_create_rand_strand(int length) {
  Strand *s = new_strand();

  for (int i = 0; i < length; i++) {
    s->push_back(s, rand() % (NUCLEOTIDE_SIZE - 1));
  }

  return s;
}

void test_tr_system_consensus() {
  PRINT_TEST_FUNC();

  Sequencer *seq = new_sequencer(ser_mid);
  Strand *s = test_create_rand_strand(60);
  vector_t *anchor = new_vector(sizeof(double));

  anchor->push_back(anchor, &(double){0.1});
  anchor->push_back(anchor, &(double){0.2});
  anchor->push_back(anchor, &(double){0.9});
  seq->set_strand(seq, s);
  seq->set_anchor(seq, anchor);

  int nreads = 20;
  Read *reads[20];
  TRSystem *trs = new_tr_system(2);
  
  printf("original strand:\n");
  printf("  "), print_strand(s);
  printf("\n");

  // int err_read_cnt = 0;
  for (int i = 0; i < nreads; i++) {
    reads[i] = seq->process(seq);
    // print_read(reads[i]);
    // vector_t *err = reads[i]->error;
    // for (int j = 0; j < err->size(err) - 2; j++) {
    //   if ((err->at(err, j) != READ_ERROR_NONE &&
    //       err->at(err, j + 1) != READ_ERROR_NONE) ||
    //       (err->at(err, j) != READ_ERROR_NONE &&
    //       err->at(err, j + 2) != READ_ERROR_NONE) ||
    //       (err->at(err, j + 1) != READ_ERROR_NONE &&
    //       err->at(err, j + 2) != READ_ERROR_NONE) ||
    //       (err->at(err, j) != READ_ERROR_NONE &&
    //       err->at(err, j + 1) != READ_ERROR_NONE &&
    //       err->at(err, j + 2) != READ_ERROR_NONE))
    //   {
    //     err_read_cnt++;
    //     break;
    //   }
    // }
  }
  // printf("%d \n", err_read_cnt);
  // printf("\n");

  trs->set_reads(trs, reads, nreads);

  Strand *consensus = trs->get_consensus(trs);
  printf("consensus strand:\n");
  printf("  "), print_strand(consensus);

  free_sequencer(seq);
  free_strand(s);
  free_vector(anchor);
  for (int i = 0; i < nreads; i++) {
    free_read(reads[i]);
  }
  free_tr_system(trs);
  free_strand(consensus);

  printf("\n");
}

int main() {
  PRINT_TEST_FILE();

  srand((unsigned int)time(NULL));

  test_tr_system_consensus();
}
