#include <random>
#include <iostream>

#include "misc.hpp"
#include "minimap_ordering.hpp"

uint64_t hash(uint64_t x, const uint64_t m) {
  x = (~x + (x << 21)) & m;
  x = x ^ x >> 24;
  x = (x + (x << 3) + (x << 8)) & m;
  x = x ^ x >> 14;
  x = (x + (x << 2) + (x << 4)) & m;
  x = x ^ x >> 28;
  x = (x + (x << 31)) & m;
  return x;
}


int main(int argc, char *argv[]) {
  minimap_ordering args(argc, argv);

  const uint64_t mask    = Mask<2>::bases(args.k_arg);
  const uint64_t nb_mers = mask + 1;

  for(uint64_t m = 0; m < nb_mers; ++m)
    std::cout << hash(m, mask) << '\n';

  return 0;
}
