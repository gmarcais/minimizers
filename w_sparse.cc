#include <fstream>
#include <vector>

#include "w_sparse.hpp"
#include "minimizers.hpp"
#include "misc.hpp"

template<int BITS>
std::vector<bool> read_uhmers(const char* path, size_t k, size_t nb_mers) {
  std::vector<bool> uhmers(nb_mers, false);
  std::ifstream is(path);
  std::string line;
  while(std::getline(is, line)) {
    uint64_t m = string_to_mer<BITS>(line, k);
    uhmers[m] = true;
  }
  return uhmers;
}

template<int BITS>
int compute_sparsiness(const w_sparse& args) {
  const size_t      k        = args.k_arg;
  const size_t      w        = args.w_arg;

  const uint64_t  mask    = (~(uint64_t)0) >> (8 * sizeof(uint64_t) - 2 * k);
  const uint64_t  nb_mers = mask + 1;

  std::vector<bool> uhmers = read_uhmers<BITS>(args.universal_arg, k, nb_mers);
  //  const std::string sequence = read_fasta(args.fasta_arg);

  std::vector<bool> window_p(w+1, false);
  fasta_iterator fasta_it(args.fasta_arg);
  auto           it   = begin(fasta_it);
  const auto     last = end(fasta_it);
  slide_mer<BITS>   mer(k);

  // Prime k-1 bases of mer
  for(size_t i = 0; i < k-1 && it != last; ++i, ++it)
    mer.appende(*it);

  // Prime w k-mers of window
  size_t uhmer_in_w = 0;
  for(size_t i = 0; i < w && it != last; ++i, ++it) {
    mer.appende(*it);
    const bool is_uhmer = uhmers[mer.mer];
    window_p[i] = uhmers[mer.mer];
    uhmer_in_w += is_uhmer;
  }

  size_t sparse_window = 0; // Number of sparse windows
  size_t wi = w;
  size_t total = 0;
  for( ; it != last; ++it, ++total) {
    mer.appende(*it);
    const bool is_uhmer = uhmers[mer.mer];
    window_p[wi] = is_uhmer;
    uhmer_in_w += is_uhmer;
    sparse_window += (uhmer_in_w == 1);
    wi = (wi + 1) % (w+1);
    uhmer_in_w -= window_p[wi];
  }

  std::cout << sparse_window << ' '
            << ((double)sparse_window / total)
            << '\n';

  return 0;
}

int main(int argc, char *argv[]) {
  w_sparse args(argc, argv);

  const int alphabet_size =
    args.size_given ? initialize_codes(args.size_arg) : initialize_codes(args.alphabet_arg);

  if(alphabet_size == 2)
    return compute_sparsiness<1>(args);
  else if(alphabet_size == 4)
    return compute_sparsiness<2>(args);
  else {
    std::cerr << "Alphabet size not supported" << std::endl;
    return 1;
  }

}
