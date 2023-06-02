#ifndef SEQUENCER_H
#define SEQUENCER_H

#include "strand.h"
#include "read.h"

typedef struct SequencerErrorRate {
  double sub; // substitution
  double del; // deleltion
  double ins; // insertion
} SequencerErrorRate;

typedef struct Sequencer {
/* public */
  struct SequencerErrorRate er;

  void (*set_error_rate)(struct Sequencer *, SequencerErrorRate);
  void (*set_strand)(struct Sequencer *, Strand *strand);
  void (*set_anchor)(struct Sequencer *, vector_t *anchor);
  Read *(*process)(struct Sequencer *);

/* private */
  Strand *_strand;
  vector_t *_anchor;

} Sequencer;

Sequencer *new_sequencer(SequencerErrorRate);
void free_sequencer(Sequencer *);

#endif
