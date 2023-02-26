#ifndef STRAND_H
#define STRAND_H

#define STRAND_INIT_SIZE 3

enum StrandType { STRAND, READ };
typedef enum StrandType StrandType;
enum Nucleotide { A, C, G, T };
typedef enum Nucleotide Nucleotide;

struct Strand {
/* public */
  int len;          // used size
  StrandType type;
  Nucleotide *seq;

/* private */
  int size;         // total size
};
typedef struct Strand Strand;

/* public */
Strand *strand_new(StrandType type);
void strand_free(Strand *s);
Strand *strand_copy(Strand *s);
int strand_resize(Strand *s, int len);
int strand_append(Strand *s, Nucleotide n);
void strand_print(Strand *s);

/* private */
int strand_expand(Strand *s);
int strand_shrink(Strand *s);

#endif
