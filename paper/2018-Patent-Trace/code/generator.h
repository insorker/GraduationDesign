#ifndef GENERATOR_H
#define GENERATOR_H

#include "strand.h"

struct ErrorRate {
  double sub;
  double ins;
  double del;
};
typedef struct ErrorRate ErrorRate;

/* public */
ErrorRate gen_er();
Strand *gen_strand_rand(int len);
Strand *gen_strand_str(char *str);
Strand *gen_strand_file(char *file);
Strand *gen_read(Strand *strand, ErrorRate er);

/* private */

#endif
