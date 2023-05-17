#include "scs.h"
#include <stdlib.h>
#include <assert.h>

Read *scs_step(TRTuple tuple, int le, int ri);

Read *scs(TRTuple tuple) {
  return scs_step(tuple);
}

Read *scs_step(TRTuple tuple, int le, int ri) {
}
