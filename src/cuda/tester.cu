#include "tester.h"
#include <iostream>
#include <vector>
#ifdef __cplusplus
extern "C" {
#endif
void testing_linkage(int a, float *b, float c) {
  std::vector<float> b_value_list;
  b_value_list.reserve(a);
  for (int i = 0; i < a; i++) {
    (*b) = (*b) + c;
    b_value_list.push_back((*b));
    std::cout << "Vector value is " << b_value_list[i] << " b value is " << (*b)
              << std::endl;
  }
  std::cout << "Final value of b is " << (*b) << std::endl;
}
#ifdef __cplusplus
}
#endif
