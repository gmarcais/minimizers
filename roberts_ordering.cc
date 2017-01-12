#include <iostream>

#include "roberts_ordering.hpp"
#include "misc.hpp"

int main(int argc, char* argv[]) {
  roberts_ordering args(argc, argv);

  static const uint64_t odd_bits   = 0x5555555555555555ULL;
  static const uint64_t even_bits  = 0xaaaaaaaaaaaaaaaaULL;
  static const uint64_t odd_bases  = 0x3333333333333333ULL;
  static const uint64_t even_bases = 0xccccccccccccccccULL;

  const uint64_t mask    = Mask<2>::bases(args.k_arg);
  const uint64_t nb_mers = mask + 1;
  const uint64_t obases = odd_bases & mask;
  const uint64_t ebases = even_bases & mask;
  for(uint64_t m = 0; m < nb_mers; ++m) {
    const uint64_t nm =
      (ebases & ((~m & even_bits) | (m & odd_bits)))
      | (obases & ((m & even_bits) | (~m & odd_bits)));
    std::cout << nm << '\n';
  }

  return 0;
}
