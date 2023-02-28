#ifndef STRAND_H
#define STRAND_H

#define STRAND_INIT_SIZE 3
#define NUCLEOTIDE_NUM   5

enum StrandType { STRAND, READ, CONSENSUS };
typedef enum StrandType StrandType;
enum Nucleotide { A, C, G, T, N };
typedef enum Nucleotide Nucleotide;
enum StrandState { CREDIBLE, VARIANT, OMITTED };
typedef enum StrandState StrandState;

struct Strand {
/* public */
  int len;            // used size
  StrandType type;
  Nucleotide *seq;
  StrandState state;  // default: CREDIBLE

/* private */
  int size;           // total size

/* TODO 
  double confidence; */
};
typedef struct Strand Strand;

/* public */
Strand *strand_new(StrandType type);
void strand_free(Strand *s);
void strand_clear(Strand *s);
Strand *strand_copy(Strand *s);
int strand_resize(Strand *s, int len);
int strand_append(Strand *s, Nucleotide n);
int strand_cmp_editdistance(Strand *s1, Strand *s2);
void strand_print(Strand *s);

void nucleotide_print(Nucleotide *n, int len);
int nucleotide_cmp(Nucleotide *n1, Nucleotide *n2, int l1, int l2);

/* private */
int strand_expand(Strand *s);
int strand_shrink(Strand *s);

#endif
