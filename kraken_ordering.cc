#include <random>
#include <iostream>

#include "misc.hpp"
#include "kraken_ordering.hpp"

int main(int argc, char *argv[]) {
  kraken_ordering args(argc, argv);

  std::mt19937_64 rng;
  seed_prg(rng);

  const uint64_t mask    = Mask<2>::bases(args.k_arg);
  const uint64_t nb_mers = mask + 1;
  const uint64_t rnd = rng() & mask;

  for(uint64_t m = 0; m < nb_mers; ++m)
    std::cout << (m ^ rnd) << '\n';

  return 0;
}
