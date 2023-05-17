#include "ml_scs.h"
#include "read.h"
#include "read_tuple.h"
#include "nucleotide.h"
#include <malloc.h>

Read *scs_helper(ReadTuple tuple);
Read *scs_computer(Read *r1, Read *r2);
Read *scs(ReadTuple tuple) {
  return scs_helper(tuple);
}

Read *scs_helper(ReadTuple tuple) {
  if (tuple.size == 1) {
    return copy_read(tuple.elem[0]);
  }
  else if (tuple.size == 2) {
    return scs_computer(tuple.elem[0], tuple.elem[1]);
  }
  else {
    int mid = tuple.size >> 1;
    ReadTuple le_tuple = { mid };
    ReadTuple ri_tuple = { tuple.size - mid };

    for (int i = 0; i <= mid; i++) {
      le_tuple.elem[i] = tuple.elem[i];
    }
    for (int i = mid + 1; i < tuple.size; i++) {
      ri_tuple.elem[i - mid - 1] = tuple.elem[i];
    }

    Read *le_result = scs_helper(le_tuple);
    Read *ri_result = scs_helper(ri_tuple);
    Read *result = scs_computer(le_result, ri_result);

    free_read(le_result);
    free_read(ri_result);
    
    return result;
  }
}

Read *scs_computer(Read *r1, Read *r2) {
  int m = r1->size(r1);
  int n = r2->size(r2);
  Nucleotide **dp = (Nucleotide **)malloc((m + 1) * sizeof(Nucleotide *));
  for (int i = 0; i <= m; i++) {
    dp[i] = (Nucleotide *)malloc((n + 1) * sizeof(Nucleotide));
  }

  for (int i = 0; i <= m; i++) dp[i][0] = i;
  for (int i = 0; i <= n; i++) dp[0][i] = i;
  
  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= n; j++) {
      if (r1->at(r1, i - 1) == r2->at(r2, j - 1)) {
        dp[i][j] = dp[i - 1][j - 1] + 1;
      }
      else {
        dp[i][j] = (dp[i - 1][j] < dp[i][j - 1])
          ? (dp[i - 1][j] + 1) : (dp[i][j - 1] + 1);
      }
    }
  }

  Read *result = new_read();
  Read *result_reverse = new_read();

  int i = m, j = n;
  while (i > 0 && j > 0) {
    if (r1->at(r1, i - 1) == r2->at(r2, j - 1)) {
      result_reverse->push_back(result_reverse, r1->at(r1, i - 1));
      i--, j--;
    }
    else if (dp[i - 1][j] < dp[i][j - 1]) {
      result_reverse->push_back(result_reverse, r1->at(r1, i - 1));
      i--;
    }
    else {
      result_reverse->push_back(result_reverse, r2->at(r2, j - 1));
      j--;
    }
  }
  while (i > 0) {
    result_reverse->push_back(result_reverse, r1->at(r1, i - 1));
    i--;
  }
  while (j > 0) {
    result_reverse->push_back(result_reverse, r2->at(r2, j - 1));
    j--;
  }

  for (int i = result_reverse->size(result_reverse) - 1; i >= 0; i--) {
    result->push_back(result, result_reverse->at(result_reverse, i));
  }
  
  free_read(result_reverse);
  for (int i = 0; i <= m; i++) {
    free(dp[i]);
  }
  free(dp);

  return result;
}
