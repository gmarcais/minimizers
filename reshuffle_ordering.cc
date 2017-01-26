#include <iostream>


#include "reshuffle_ordering.hpp"
#include "misc.hpp"

template<int BITS>
int shuffle_minimizers(const reshuffle_ordering& args) {
const size_t      k        = args.k_arg;
  const uint64_t  mask    = Mask<BITS>::bases(k);
  const uint64_t  nb_mers = mask + 1;

  std::vector<uint64_t> order(nb_mers, nb_mers);

  std::ifstream is(args.minimizers_arg);
  if(!is.good())
    reshuffle_ordering::error() << "Failed to open minimizers file '" << args.minimizers_arg << '\'';

  // Read minimizers
  std::string line;
  uint64_t o = 0;
  while(std::getline(is, line)) {
    uint64_t m = string_to_mer<BITS>(line, k);
    order[m] = o++;
  }

  // Create random order of thos
  std::mt19937_64 prg;
  seed_prg(prg,
           args.save_seed_given ? args.save_seed_arg : nullptr,
           args.load_seed_given ? args.load_seed_arg : nullptr);
  std::vector<uint64_t> rorder(o, 0);
  for(size_t i = 0; i < rorder.size(); ++i)
    rorder[i] = i;
  std::shuffle(rorder.begin(), rorder.end(), prg);

  for(uint64_t m = 0; m < nb_mers; ++m) {
    if(order[m] == nb_mers) { // Not a minimizer, come after
      std::cout << o++ << '\n';
    } else { // A minimizer, low randomized order
      std::cout << rorder[order[m]] << '\n';
    }
  }

  return 0;
}

int main(int argc, char *argv[]) {
  std::ios::sync_with_stdio(false);
  reshuffle_ordering args(argc, argv);

  const int alphabet_size =
    args.size_given ? initialize_codes(args.size_arg) : initialize_codes(args.alphabet_arg);

  if(alphabet_size == 2)
    return shuffle_minimizers<1>(args);
  else if(alphabet_size == 4)
    return shuffle_minimizers<2>(args);
  else {
    std::cerr << "Alphabet size not supported" << std::endl;
    return 1;
  }

  return 0;
}
