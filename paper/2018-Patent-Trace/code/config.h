#ifndef CONFIG_H
#define CONFIG_H

struct ErrorRate {
  double sub; // substitution
  double del; // deletion
  double ins; // insertion
};
typedef struct ErrorRate ErrorRate;

struct TRConfig {
  int num;        // number of reads
  int len;        // length of strand
  int width;      // width of look-ahead window
  ErrorRate er;   // error rate
};
typedef struct TRConfig TRConfig;

#endif
