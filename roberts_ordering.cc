#include <iostream>

#include "roberts_ordering.hpp"
#include "misc.hpp"

int main(int argc, char* argv[]) {
  roberts_ordering args(argc, argv);

  static const uint64_t odd_bits   = 0x5555555555555555ULL;
  static const uint64_t even_bits  = 0xaaaaaaaaaaaaaaaaULL;
  static const uint64_t odd_bases  = 0x3333333333333333ULL;
  static const uint64_t even_bases = 0xddddddddddddddddULL;

  uint64_t nb_mers = Mask<2>::bases(args.k_arg) + 1;
  for(uint64_t m = 0; m < nb_mers; ++m) {
    uint64_t nm =
      (even_bases & ((~m & even_bits) | (m & odd_bits)))
      | (odd_bases & ((m & even_bits) | (~m & odd_bits)));
    std::cout << nm;
  }

  return 0;
}
