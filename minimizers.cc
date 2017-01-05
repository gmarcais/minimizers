#include <fstream>
#include <iostream>

#include "minimizers.hpp"
#include "misc.hpp"
#include "compute_minimizers.hpp"

std::vector<uint64_t> order;
struct order_greater {
  bool operator()(uint64_t x, uint64_t y) const { return order[y] < order[x]; }
};
struct order_lesser {
  bool operator()(uint64_t x, uint64_t y) const { return order[x] < order[y]; }
};

template<int BITS>
int order_minimizer(const minimizers& args) {
  const std::string sequence = read_fasta(args.fasta_arg);
  const size_t      k        = args.k_arg;
  mer_pos::window            = args.w_arg;
  const size_t      w        = mer_pos::window;

  const uint64_t  mask    = Mask<BITS>::bases(k);
  const uint64_t  nb_mers = mask + 1;
  std::mt19937_64 prg;
  seed_prg(prg);

  std::unordered_set<uint64_t> minimizers;

  order.resize(nb_mers, (uint64_t)0);

  // Create counters for each class
  std::vector<size_t> counters;
  counters.resize(1);
  counters[0] = nb_mers;

  // Add non-universal mers as the highest elements in the ordering
  if(args.universal_given) {
    const int nonu_counter = args.conjugacy_flag ? k : 1;
    counters.resize(nonu_counter + 1, (size_t)0);
    std::fill(order.begin(), order.end(), nonu_counter);
    counters[0]		   = 0;
    counters[nonu_counter] = nb_mers;

    std::ifstream is(args.universal_arg);
    if(!is.good())
      minimizers::error() << "Failed to open universal k-mer set file '" << args.universal_arg << '\'';
    std::string line;
    while(std::getline(is, line)) {
      uint64_t m = string_to_mer<BITS>(line, k);
      order[m] = 0;
      ++counters[0];
      --counters[nonu_counter];
    }
    if(!is.eof())
      minimizers::error() << "Error while reading universal k-mer set file '" << args.universal_arg << '\'';
  }

  // Create conjugacy classes for elements that are not yet classified
  // (i.e. order[m] == 0).
  if(args.conjugacy_flag) {
    if(counters.size() < k)
      counters.resize(k, (size_t)0);
    for(uint64_t m = 0; m < nb_mers; ++m) {
      if(order[m] == 0) {
	const auto o = compute_order<BITS>(m, k);
	order[m] = k - o;
	--counters[0];
	++counters[k - o];
      }
    }
  }

  // Compute partial sums
  size_t sum = 0;
  for(size_t i = 0; i < counters.size(); ++i) {
    size_t t = counters[i];
    counters[i] = sum;
    sum += t;
  }

  if(args.debug_flag) {
    std::cerr << "counters:";
    for(const auto s : counters)
      std::cerr << ' ' << s;
    std::cerr << ' ' << sum << '\n';
    std::cerr << "counters:\n";
    for(uint64_t m = 0; m < nb_mers; ++m)
      std::cerr << mer_to_string<BITS>(m, k) << ':' << order[m] << '\n';
  }

  // Assign final order
  if(args.randomize_flag) {
    std::vector<uint64_t> visit(nb_mers);
    for(uint64_t m = 0; m < nb_mers; ++m)
      visit[m] = m;
    std::shuffle(visit.begin(), visit.end(), prg);
    for(uint64_t m : visit)
      order[m] = counters[order[m]]++;
  } else {
    for(uint64_t m = 0; m < nb_mers; ++m)
      order[m] = counters[order[m]]++;
  }

  if(args.debug_flag) {
    std::cerr << "Order: " << (void*)&order << "\n";
    for(uint64_t m = 0; m < nb_mers; ++m)
      std::cerr << mer_to_string<BITS>(m, k) << ':' << order[m] << '\n';
  }

  auto ms = compute_minimizers<order_lesser, BITS>()(sequence, k, minimizers);

  if(args.minimizers_given) {
    std::ofstream os(args.minimizers_arg);
    for(const auto m : minimizers) {
      os << mer_to_string<BITS>(m, k);
      if(args.debug_flag)
        os << ' ' << order[m];
      os << '\n';
    }
  }

  //  const double expected = ((double)(nb_mers - 1) / (double)(nb_mers + w)) * (2.0 / (w + 1));
  const double expected = 2.0 / (w + 1);
  const double actual = (double)ms.nb() / (double)(sequence.size() - k + 1);
  std::cout << "minimizers: " << minimizers.size() << '\n'
            << "mean: " << ms.mean() << ' ' << ms.sum() << '/' << ms.nb() << '\n'
            << "stddev: " << ms.stddev() << '\n'
            << "density: " << actual << (actual < expected ? " < " : " > ")  << expected << '\n'
            << '\t' << ms.nb() << '/' << (sequence.size() - k + 1) << ' ' << ms.total() << '\n';

  return 0;
}

int main(int argc, char *argv[]) {
  const minimizers args(argc, argv);

  const int alphabet_size =
    args.size_given ? initialize_codes(args.size_arg) : initialize_codes(args.alphabet_arg);

  if(alphabet_size == 2)
    return order_minimizer<1>(args);
  else if(alphabet_size == 4)
    return order_minimizer<2>(args);
  else {
    std::cerr << "Alphabet size not supported" << std::endl;
    return 1;
  }

  return 0;
}
