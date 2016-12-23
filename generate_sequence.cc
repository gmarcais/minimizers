#include <iostream>
#include <random>

#include "generate_sequence.hpp"
#include "misc.hpp"

int main(int argc, char* argv[]) {
  std::ios::sync_with_stdio(false);
  generate_sequence args(argc, argv);

  const char* bases;
  std::string alphabet;
  if(args.size_given) {
    std::ostringstream os;
    for(int i = 0; i < args.size_arg; ++i)
      os << i;
    alphabet = os.str();
    bases = alphabet.c_str();
  } else {
    bases = args.alphabet_arg;
  }

  std::mt19937_64 rng;
  seed_prg(rng);
  std::uniform_int_distribution<int> base(0, strlen(bases) - 1);

  std::cout << ">sequence" << args.n_arg << '\n';
  for(size_t i = 0; i < args.n_arg; ++i)
    std::cout << bases[base(rng)];
  std::cout << '\n';

  return 0;
}
