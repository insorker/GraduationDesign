#include "sequencer.h"
#include <malloc.h>
#include <stdlib.h>

Read *sequencer_process(Sequencer *seq, Strand *s);

Sequencer *new_sequencer(SequencerErrorRate er) {
  Sequencer *seq = (Sequencer *)malloc(sizeof(Sequencer));

  seq->er.sub = er.sub;
  seq->er.del = er.del + seq->er.sub;
  seq->er.ins = er.ins + seq->er.del;

  seq->process = sequencer_process;

  return seq;
}

void free_sequencer(Sequencer *c) {
  free(c);
}

Read *sequencer_process(Sequencer *seq, Strand *s) {
  Read *read = new_read();

  for (int i = 0; i < s->size(s); i++) {
    double rate = (double)rand() / RAND_MAX;

    // substitution
    if (seq->er.sub >= rate) {
      Nucleotide n = rand_nucleotide();
      
      while (n == s->at(s, i)) {
        n = rand_nucleotide();
      }

      read->push_back(read, n);
      read->error->push_back(read->error, READ_ERROR_SUB);
    }
    // deletion
    else if (seq->er.del >= rate) {
      read->error->push_back(read->error, READ_ERROR_DEL);
      continue;
    }
    // insertion
    else if (seq->er.ins >= rate) {
      read->error->push_back(read->error, READ_ERROR_INS);
      read->push_back(read, rand_nucleotide());
      i--;
    }
    // no error
    else {
      read->error->push_back(read->error, READ_ERROR_NONE);
      read->push_back(read, s->at(s, i));
    }
  }

  return read;
}
