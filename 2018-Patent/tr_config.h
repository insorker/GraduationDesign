#ifndef TR_CONFIG_H
#define TR_CONFIG_H

typedef struct TRErrorRate {
  double sub; // substitution
  double del; // deletion
  double ins; // insertion
} TRErrorRate;

typedef struct TRConfig {
  int num_reads;     // num of reads
  int len_strand;    // length of strand
  int width_window;  // width of look-ahead window
  TRErrorRate er;    // error rate
} TRConfig;

void read_tr_config(TRConfig *tr_config);

#endif
