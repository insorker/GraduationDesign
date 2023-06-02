#include "sequencer.h"
#include <malloc.h>
#include <stdlib.h>

void sequencer_set_error_rate(Sequencer *seq, SequencerErrorRate er);
void sequencer_set_strand(Sequencer *seq, Strand *strand);
void sequencer_set_anchor(Sequencer *seq, vector_t *anchor);
Read *sequencer_process(Sequencer *seq);

Sequencer *new_sequencer(SequencerErrorRate er) {
  Sequencer *seq = (Sequencer *)malloc(sizeof(Sequencer));

  seq->er.sub = er.sub;
  seq->er.del = er.del + seq->er.sub;
  seq->er.ins = er.ins + seq->er.del;

  seq->set_error_rate = sequencer_set_error_rate;
  seq->set_strand = sequencer_set_strand;
  seq->set_anchor = sequencer_set_anchor;
  seq->process = sequencer_process;

  seq->_strand = new_strand();
  seq->_anchor = new_vector(sizeof(double));

  return seq;
}

void free_sequencer(Sequencer *seq) {
  free_strand(seq->_strand);
  free_vector(seq->_anchor);
  free(seq);
}

void sequencer_set_error_rate(Sequencer *seq, SequencerErrorRate er) {
  seq->er.sub = er.sub;
  seq->er.del = er.del + seq->er.sub;
  seq->er.ins = er.ins + seq->er.del;
}

void sequencer_set_strand(Sequencer *seq, Strand *strand) {
  seq->_strand->clear(seq->_strand);

  for (int i = 0; i < strand->size(strand); i++) {
    seq->_strand->push_back(
        seq->_strand,
        strand->at(strand, i));
  }
}

void sequencer_set_anchor(Sequencer *seq, vector_t *anchor) {
  seq->_anchor->clear(seq->_anchor);

  for (int i = 0; i < anchor->size(anchor); i++) {
    seq->_anchor->push_back(
        seq->_anchor,
        anchor->at(anchor, i));
  }
}

Read *sequencer_process(Sequencer *seq) {
  Strand *s = seq->_strand;
  vector_t *a = seq->_anchor;
  int pa = 0;

  Read *read = new_read();

  for (int i = 0; i < s->size(s); i++) {
    double rate = (double)rand() / RAND_MAX;

    if (pa < a->size(a)) {
      int anchor_pos = *(double *)a->at(a, pa) * s->size(s);

      if (anchor_pos < i) {
        pa++;
      }
      else if (anchor_pos == i)
      {
        read->anchor->push_back(
            read->anchor,
            &(int){read->size(read)});
        read->error->push_back(read->error, &(ReadErrorType){READ_ERROR_NONE});
        read->push_back(read, s->at(s, i));
        continue;
      }
    }

    // substitution
    if (seq->er.sub >= rate) {
      Nucleotide n = rand_nucleotide();
      
      while (n == s->at(s, i)) {
        n = rand_nucleotide();
      }

      read->push_back(read, n);
      read->error->push_back(read->error, &(ReadErrorType){READ_ERROR_SUB});
    }
    // deletion
    else if (seq->er.del >= rate) {
      read->error->push_back(read->error, &(ReadErrorType){READ_ERROR_DEL});
      continue;
    }
    // insertion
    else if (seq->er.ins >= rate) {
      read->error->push_back(read->error, &(ReadErrorType){READ_ERROR_INS});
      read->push_back(read, rand_nucleotide());
      i--;
    }
    // no error
    else {
      read->error->push_back(read->error, &(ReadErrorType){READ_ERROR_NONE});
      read->push_back(read, s->at(s, i));
    }
  }

  return read;
}
