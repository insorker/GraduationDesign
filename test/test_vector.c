#include "test.h"
#include "vector.h"

void test_vector_push_pop_clear() {
  PRINT_TEST_FUNC();

  Vector *vec = new_vector();

  vec->push_back(vec, 1);
  printf("%d\n", vec->at(vec, 0));
  printf("%d\n", vec->pop_back(vec));

  for (int i = 0; i < 100; i++) {
    vec->push_back(vec, i);
  }
  for (int i = 0; i < 100; i++) {
    printf("%d ", vec->at(vec, i));
  }
  printf("\n");
  printf("%d\n", vec->size(vec));
  vec->clear(vec);
  printf("%d\n", vec->size(vec));

  free_vector(vec);

  printf("\n");
}

int main() {
  PRINT_TEST_FILE();

  test_vector_push_pop_clear();
}
