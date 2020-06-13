#include "../config.h"

#include <stdlib.h>

struct peano_hilbert_data {
  size_t key;
  int index;
  long long count;
};

size_t peano_hilbert_key(int x, int y, int z, int bits);
int peano_compare(const void* a, const void* b);
