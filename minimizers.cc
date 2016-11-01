#include <fstream>
#include <iostream>

#include "minimizers.hpp"
#include "misc.hpp"
#include "compute_minimizers.hpp"

std::vector<uint64_t> order;
struct order_greater {
  bool operator()(uint64_t x, uint64_t y) const { return order[y] < order[x]; }
};

int main(int argc, char *argv[]) {
  const minimizers args(argc, argv);

  const std::string sequence = read_fasta(args.fasta_arg);
  const size_t      k        = args.k_arg;
  mer_pos::window            = args.w_arg;

  const uint64_t  mask    = (~(uint64_t)0) >> (8 * sizeof(uint64_t) - 2 * k);
  const uint64_t  nb_mers = mask + 1;
  std::mt19937_64 prg;
  seed_prg(prg);

  std::unordered_set<uint64_t> minimizers;

  order.resize(nb_mers, (uint64_t)0);

  // Create counters for each class
  std::vector<size_t> counters;
  if(args.conjugacy_flag) {
    counters.resize(k, (size_t)0);
    for(uint64_t m = 0; m < nb_mers; ++m) {
      const auto o = compute_order(m, k);
      order[m] = k - o;
      ++counters[k - o];
    }
  } else {
    counters.resize(1);
    counters[0] = nb_mers;
    //    std::fill(order.begin(), order.end(), 0);
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
      std::cerr << mer_to_string(m, k) << ':' << order[m] << '\n';
  }

  // Assign final order
  if(args.randomize_flag) {
    std::vector<uint64_t> visit(nb_mers);
    for(uint64_t m = 0; m < nb_mers; ++m)
      visit[m] = m;
    std::shuffle(visit.begin(), visit.end(), prg);
    for(uint64_t m : visit)
      order[m] = ++counters[order[m]];
  } else {
    for(uint64_t m = 0; m < nb_mers; ++m)
      order[m] = ++counters[order[m]];
  }

  if(args.debug_flag) {
    std::cerr << "Order:\n";
    for(uint64_t m = 0; m < nb_mers; ++m)
      std::cerr << mer_to_string(m, k) << ':' << order[m] << '\n';
  }

  auto ms = compute_minimizers<order_greater>()(sequence, k, minimizers);

  if(args.minimizers_given) {
    std::ofstream os(args.minimizers_arg);
    for(const auto m : minimizers)
      os << mer_to_string(m, k) << '\n';
  }

  double expected = ((double)(nb_mers - 1) / (double)(nb_mers + mer_pos::window)) * (2.0 / (mer_pos::window + 1));
  double actual = (double)ms.nb() / (double)(sequence.size() - k + 1);
  std::cout << "minimizers: " << minimizers.size() << '\n'
            << "mean: " << ms.mean() << '\n'
            << "stddev: " << ms.stddev() << '\n'
            << "density: " << actual << (actual < expected ? " < " : " > ")  << expected << '\n';

  return 0;
}
