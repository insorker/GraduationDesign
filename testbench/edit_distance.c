#include "../test/test.h"
#include "tr_system.h"
#include "sequencer.h"
#include <stdlib.h>
#include <time.h>

SequencerErrorRate ser_null = { 0, 0, 0 };
SequencerErrorRate ser_sm = { 0.05, 0.03, 0.01 };
SequencerErrorRate ser_mid = { 0.05, 0.05, 0.05 };
SequencerErrorRate ser_lg = { 0.1, 0.2, 0.35 };
SequencerErrorRate ser_full = { 0.3, 0.3, 0.4 };

int min(int a, int b, int c) {
    int min = a;
    if (b < min) {
        min = b;
    }
    if (c < min) {
        min = c;
    }
    return min;
}

int edit_distance(Strand *a, Strand *b) {
  int m = a->size(a), n = b->size(b);
  int dp[200][200];

  for (int i = 0; i <= m; i++) {
    for (int j = 0; j <= n; j++) {
      if (i == 0) {
        dp[i][j] = j;
      } else if (j == 0) {
        dp[i][j] = i;
      } else if (a->at(a, i - 1) == b->at(b, j - 1)) {
        dp[i][j] = dp[i - 1][j - 1];
      } else {
        dp[i][j] = 1 + min(dp[i - 1][j - 1], dp[i - 1][j], dp[i][j - 1]);
      }
    }
  }

  return dp[m][n];
}

Strand *test_create_strand(char *str, int length) {
  Strand *s = new_strand();

  for (int i = 0; i < length; i++) {
    switch (str[i]) {
      case 'A': s->push_back(s, NUCLEOTIDE_A); break;
      case 'C': s->push_back(s, NUCLEOTIDE_C); break;
      case 'G': s->push_back(s, NUCLEOTIDE_G); break;
      case 'T': s->push_back(s, NUCLEOTIDE_T); break;
    }
  }

  return s;
}

void test_tr_system_consensus() {
  PRINT_TEST_FUNC();

  SequencerErrorRate ser = { 0.1, 0.1, 0.1 };
  Sequencer *seq = new_sequencer(ser_mid);
  Strand *s = test_create_strand("GTTGACCCCACGTCGTAGATCGACCGATCTGGAGCTAGCTAGTCGATTTCACATTCAGTGTGCTAGTACACTTTTTCTATCTAAGATCATACGGCTATCGGATTCGATTTCAA", 100);
  TRSystem *trs = new_tr_system(5);
  
  printf("original strand:\n");
  printf("- "), print_strand(s);
  printf("\n");
  
  for (double r = 0.01; r <= 0.15; r += 0.01) {
    // average
    ser.del = r / 3;
    ser.ins = r / 3;
    ser.sub = r / 3;

    // ser.sub = r / 10 * 8;
    // ser.del = r / 10 * 1;
    // ser.ins = r / 10 * 1;

    // ser.sub = r / 10 * 1;
    // ser.del = r / 10 * 8;
    // ser.ins = r / 10 * 1;

    // ser.sub = r / 10 * 1;
    // ser.del = r / 10 * 1;
    // ser.ins = r / 10 * 8;

    int nreads = 20;
    Read *reads[20];

    seq->set_error_rate(seq, ser);
    for (int i = 0; i < nreads; i++) {
      reads[i] = seq->process(seq, s);
      // print_read(reads[i]);
    }
    // printf("\n");

    trs->set_reads(trs, reads, nreads);
    Strand *consensus = trs->get_consensus(trs);
    // printf("consensus strand:\n");
    // printf("- "), print_strand(consensus);
    // printf("\n");
    // printf("edit distance: %d\n", edit_distance(s, consensus));
    printf("%d ", edit_distance(s, consensus));

    for (int i = 0; i < nreads; i++) {
      free_read(reads[i]);
    }
    free_strand(consensus);
  }

  free_sequencer(seq);
  free_strand(s);
  free_tr_system(trs);

  printf("\n");
}

int main() {
  PRINT_TEST_FILE();

  srand((unsigned int)time(NULL));

  test_tr_system_consensus();
}
