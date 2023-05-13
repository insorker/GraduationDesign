#include "vector.h"
#include <malloc.h>
#include <limits.h>
#include <assert.h>

#define VECTOR_MIN_CAPACITY 3
#define VECTOR_MAX_CAPACITY INT_MAX
#define VECTOR_MAX_SIZE (VECTOR_MAX_CAPACITY / 4) // 4 是大于 3 的最大二次幂

int  vector_size(Vector *);
int  vector_at(Vector *, int index);
void vector_push_back(Vector *, int value);
int  vector_pop_back(Vector *);
void vector_clear(Vector *);

void vector_set_capacity(Vector *, int nsize);
void vector_expand(Vector *);
void vector_shrink(Vector *);


Vector *new_vector(int t_size) {
  Vector *vec = (Vector *)malloc(sizeof(Vector));

/* public */
  vec->size = vector_size;
  vec->at = vector_at;
  vec->push_back = vector_push_back;
  vec->pop_back = vector_pop_back;
  vec->clear = vector_clear;

/* private */
  vec->_size = 0;
  vec->_capacity = VECTOR_MIN_CAPACITY;
  vec->_elem = (int *)malloc(vec->_capacity * sizeof(int));

  vec->set_capacity = vector_set_capacity;
  vec->expand = vector_expand;
  vec->shrink = vector_shrink;

  return vec;
}

void free_vector(Vector *vec) {
  free(vec->_elem);
  free(vec);
}

int compare_vector(Vector *vec_a, Vector *vec_b) {
  if (vec_a->size(vec_a) != vec_b->size(vec_b)) {
    return 1;
  }

  for (int i = 0; i < vec_a->size(vec_a); i++) {
    if (vec_a->at(vec_a, i) != vec_b->at(vec_b, i)) {
      return 1;
    }
  }

  return 0;
}


int vector_size(Vector *vec) {
  return vec->_size;
}

int  vector_at(Vector *vec, int index) {
  assert(index >= 0);
  assert(index < vec->_size);

  return vec->_elem[index];
}

void vector_push_back(Vector *vec, int value) {
  vec->set_capacity(vec, vec->_size + 1);
  vec->_elem[vec->_size] = value;
  vec->_size += 1;
}

int vector_pop_back(Vector *vec) {
  int value = vec->_elem[vec->_size - 1];

  vec->set_capacity(vec, vec->_size - 1);
  vec->_size -= 1;

  return value;
}

void vector_clear(Vector *vec) {
  while (vec->_size) {
    vec->pop_back(vec);
  }
}


void vector_set_capacity(Vector *vec, int nsize) {
  assert(nsize >= 0);
  assert(nsize <= VECTOR_MAX_SIZE);

  while (vec->_capacity <= nsize)
    vec->expand(vec);
  while (vec->_capacity > VECTOR_MIN_CAPACITY
      && vec->_capacity > nsize * 2)
    vec->shrink(vec);
}

void vector_expand(Vector *vec) {
  vec->_capacity = vec->_capacity * 2;
  vec->_elem = (int *)realloc(vec->_elem, vec->_capacity * sizeof(int));
}

void vector_shrink(Vector *vec) {
  vec->_capacity = vec->_capacity / 2;
  vec->_elem = (int *)realloc(vec->_elem, vec->_capacity * sizeof(int));
}
