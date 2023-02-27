#ifndef GENERATOR_H
#define GENERATOR_H

#include "strand.h"
#include "config.h"

/* public */
ErrorRate gen_er();
TRConfig gen_trc();
Strand *gen_strand(int len);
Strand *gen_read(Strand *strand, ErrorRate er);

/* private */

#endif
