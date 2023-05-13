#ifndef VECTOR_H
#define VECTOR_H

typedef struct Vector {
/* public */
  int  (*size)(struct Vector *);
  int  (*at)(struct Vector *, int index);
  void (*push_back)(struct Vector *, int value);
  int  (*pop_back)(struct Vector *);
  void (*clear)(struct Vector *);

/* private */
  int _size;      // actual size
  int _capacity;  // total capacity
  int *_elem;     // elements

  void (*set_capacity)(struct Vector *, int nsize);
  void (*expand)(struct Vector *);
  void (*shrink)(struct Vector *);
} Vector;

Vector *new_vector();
void free_vector(Vector *);
int compare_vector(Vector *, Vector *);

#endif
