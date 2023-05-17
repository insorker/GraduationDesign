#ifndef Read_TUPLE_H
#define Read_TUPLE_H

#include "read.h"

typedef struct ReadTuple {
  int size;
  Read *elem[5];
} ReadTuple;

#endif
