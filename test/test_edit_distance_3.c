#include "test.h"
#include "tr_system.h"
#include "sequencer.h"
#include <stdlib.h>
#include <time.h>

#define WIDTH 2
#define TYPE  0
#define TIMES 100
#define MAX_RATE 15

int min(int a, int b, int c);
int edit_distance(Strand *a, Strand *b);
void test_tr_system_consensus();

int main() {
  PRINT_TEST_FILE();

  srand((unsigned int)time(NULL));

  test_tr_system_consensus();
}

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

Strand *test_create_rand_strand(int length) {
  Strand *s = new_strand();

  for (int i = 0; i < length; i++) {
    s->push_back(s, rand() % (NUCLEOTIDE_SIZE - 1));
  }

  return s;
}

void test_tr_system_consensus() {
  PRINT_TEST_FUNC();

  int result[MAX_RATE + 1] = {};

  SequencerErrorRate ser = {};
  Sequencer *seq = new_sequencer(ser);
  vector_t *anchor = new_vector(sizeof(double));

  for (int i = 0; i < TIMES; i++) {
    for (int j = 1; j <= MAX_RATE; j++) {
      double r = (double)j / 100;
#if TYPE == 0
      ser.sub = r / 3;
      ser.del = r / 3;
      ser.ins = r / 3;
#elif TYPE == 1
      ser.sub = r / 10 * 8;
      ser.del = r / 10 * 2;
      ser.ins = r / 10 * 2;
#elif TYPE == 2
      ser.sub = r / 10 * 2;
      ser.del = r / 10 * 8;
      ser.ins = r / 10 * 2;
#elif TYPE == 3
      ser.sub = r / 10 * 2;
      ser.del = r / 10 * 2;
      ser.ins = r / 10 * 8;
#endif
    
      Strand *s = test_create_rand_strand(100);
      seq->set_strand(seq, s);
      anchor->push_back(anchor, &(double){0.1});
      anchor->push_back(anchor, &(double){0.2});
      anchor->push_back(anchor, &(double){0.9});
      seq->set_anchor(seq, anchor);
      seq->set_error_rate(seq, ser);

      int nreads = 20;
      Read *reads[20];
      TRSystem *trs = new_tr_system(WIDTH);
      
      // printf("original strand:\n");
      // printf("  "), print_strand(s);
      // printf("\n");

      for (int i = 0; i < nreads; i++) {
        reads[i] = seq->process(seq);
      }

      trs->set_reads(trs, reads, nreads);

      Strand *consensus = trs->get_consensus(trs);
      // printf("consensus strand:\n");
      // printf("  "), print_strand(consensus);
      // printf("edit_distance: %d\n", edit_distance(s, consensus));
      // printf("%d\n", edit_distance(s, consensus));
      result[j] += edit_distance(s, consensus);

      for (int i = 0; i < nreads; i++) {
        free_read(reads[i]);
      }
      free_tr_system(trs);
      free_strand(consensus);
      free_strand(s);
    }
  }

  for (int i = 1; i <= MAX_RATE; i++) {
    printf("%f ", (double)result[i] / TIMES);
  }

  free_sequencer(seq);
  free_vector(anchor);

  PRINT_TEST_FUNC();
}

