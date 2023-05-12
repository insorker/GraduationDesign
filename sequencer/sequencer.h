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

  Read *(*process)(struct Sequencer *c, Strand *s);

} Sequencer;

Sequencer *new_sequencer(SequencerErrorRate);
void free_sequencer(Sequencer *);

#endif
