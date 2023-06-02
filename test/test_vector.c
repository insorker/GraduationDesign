#include "test.h"
#include "vector.h"

void test_vector_push_pop_clear() {
  PRINT_TEST_FUNC();

  vector_t *vec = new_vector(sizeof(int));

  vec->push_back(vec, &(int){1});
  printf("%d\n", *(int *)vec->at(vec, 0));

  free_vector(vec);

  PRINT_TEST_FUNC();
}

int main() {
  PRINT_TEST_FILE();

  test_vector_push_pop_clear();
}
