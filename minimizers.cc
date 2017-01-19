#include <fstream>
#include <iostream>

#include "minimizers.hpp"
#include "misc.hpp"
#include "compute_minimizers.hpp"
#include "mean_stdev.hpp"

std::vector<uint64_t> order;
struct order_greater {
  bool operator()(uint64_t x, uint64_t y) const { return order[y] < order[x]; }
};
struct order_lesser {
  bool operator()(uint64_t x, uint64_t y) const { return order[x] < order[y]; }
};

template<int BITS>
struct order_canonical_lesser {
  bool operator()(uint64_t x, uint64_t y) const {
    return order[canonical<BITS>(x)] < order[canonical<BITS>(y)];
  }
};

template<typename T>
void write_histo(std::ostream& os, const std::vector<T>& h, bool normalize = false) {
  double sum = 1.0;
  if(normalize) {
    for(const auto& x : h)
      sum += x;
  }

  for(size_t i = 0; i < h.size(); ++i) {
    if(h[i] > 0)
      os << i << ' ' << (h[i] / sum) << '\n';
  }
}

template<int BITS>
int order_minimizer(const minimizers& args) {
  //  const std::string sequence = read_fasta(args.fasta_arg);
  const size_t      k        = args.k_arg;
  mer_len                    = k;
  mer_pos::window            = args.w_arg;
  const size_t      w        = mer_pos::window;

  const uint64_t  mask    = Mask<BITS>::bases(k);
  const uint64_t  nb_mers = mask + 1;
  std::mt19937_64 prg;
  seed_prg(prg,
           args.save_seed_given ? args.save_seed_arg : nullptr,
           args.load_seed_given ? args.load_seed_arg : nullptr);

  // Data structure and action to take upon finding a minimizer
  mean_stdev distances;
  std::unordered_set<uint64_t> minimizers;
  std::vector<std::pair<size_t, size_t>> minimizers_counts(nb_mers, std::make_pair(0, 0));

  size_t min_pos    = std::numeric_limits<size_t>::max();
  size_t prev_start = 0;
  auto act = [&](const mer_pos& mp) -> void {
    if(min_pos != std::numeric_limits<size_t>::max()) {
      distances.sample(mp.pos - min_pos);
      uint64_t m = args.canonical_flag ? canonical<BITS>(mp.mer) : mp.mer;
      minimizers.insert(m);
      ++(minimizers_counts[m].first);
      minimizers_counts[m].second += mp.win_start - prev_start;
    }
    min_pos    = mp.pos;
    prev_start = mp.win_start;
  };

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
  } else if(args.order_given) {
    // Read directly an order from the given file
    std::ifstream is(args.order_arg);
    if(!is.good())
      minimizers::error() << "Failed to open order file '" << args.order_arg << '\'';
    size_t i = 0;
    std::vector<bool> used(nb_mers, false);
    for( ; i < nb_mers && is.good(); ++i) {
      is >> order[i];
      if(order[i] >= nb_mers)
        minimizers::error() << "Order " << order[i] << " is greater than the number of mers " << nb_mers;
      if(used[order[i]])
        minimizers::error() << "Order " << order[i] << " used twice";
      used[order[i]] = true;
    }
    if(i < nb_mers)
      minimizers::error() << "Too few elements in order file. Expected " << nb_mers << ", got " << i;
  } else {
    for(uint64_t m = 0; m < nb_mers; ++m)
      order[m] = counters[order[m]]++;
  }

  if(args.debug_flag) {
    std::cerr << "Order: " << (void*)&order << "\n";
    for(uint64_t m = 0; m < nb_mers; ++m)
      std::cerr << mer_to_string<BITS>(m, k) << ':' << order[m] << '\n';
  }

  fasta_iterator it(args.fasta_arg);
  size_t seen_mers = args.canonical_flag
    ? compute_minimizers<order_canonical_lesser<BITS>, BITS>()(begin(it), end(it), k, act)
    : compute_minimizers<order_lesser, BITS>()(begin(it), end(it), k, act);

  if(args.minimizers_given || args.prefix_given) {
    const std::string path = args.minimizers_given ? std::string(args.minimizers_arg) : std::string(args.prefix_arg) + "_minimizers.txt";
    std::ofstream os(path);
    for(const auto m : minimizers) {
      os << mer_to_string<BITS>(m, k);
      if(args.debug_flag)
        os << ' ' << order[m];
      os << '\n';
    }
    if(!os.good())
      minimizers::error() << "Error writing minimizers to file '" << path << '\'';
  }

  if(args.distance_histo_given || args.prefix_given) {
    const std::string path = args.distance_histo_given ? std::string(args.distance_histo_arg) : std::string(args.prefix_arg) + "_distance.histo";
    std::ofstream os(path);
    write_histo(os, distances.histo(), true);
    if(!os.good())
      minimizers::error() << "Error writing distance histo to file '" << path << '\'';
  }

  if(args.count_histo_given || args.data_histo_given || args.prefix_given) {
    std::vector<double> count_histo;
    std::vector<double> data_histo;
    size_t total = 0;
    for(size_t i = 0; i < minimizers_counts.size(); ++i) {
      if(minimizers_counts[i].first > 0) {
        const auto x = (minimizers_counts[i].first - 1) / 50;
        if(x >= count_histo.size())
          count_histo.resize(x + 1, 0);
        ++count_histo[x];
        ++total;
      }
      if(minimizers_counts[i].second > 0) {
        const auto y = (minimizers_counts[i].second - 1) / 100;
        if(y >= data_histo.size())
          data_histo.resize(y + 1, 0);
        ++data_histo[y];
      }
    }
    for(auto& x : count_histo)
      x /= total;
    for(auto& y : data_histo)
      y /= total;

    if(args.count_histo_given || args.prefix_given) {
      const std::string path = args.count_histo_given ? std::string(args.count_histo_arg) : std::string(args.prefix_arg) + "_count.histo";
      std::ofstream os(path);
      write_histo(os, count_histo);
      if(!os.good())
        minimizers::error() << "Error writing count histo to file '" << path << '\'';
    }

    if(args.data_histo_given || args.prefix_given) {
      const std::string path = args.data_histo_given ? std::string(args.data_histo_arg) : std::string(args.prefix_arg) + "_data.histo";
      std::ofstream os(path);
      write_histo(os, data_histo);
      if(!os.good())
        minimizers::error() << "Error writing data histo to file '" << path << '\'';
    }
  }

  //  const double expected = ((double)(nb_mers - 1) / (double)(nb_mers + w)) * (2.0 / (w + 1));
  const double expected = 2.0 / (w + 1);
  const double actual = (double)distances.nb() / (double)seen_mers;
  std::cout << "minimizers: " << minimizers.size() << '\n'
            << "mean: " << distances.mean() << ' ' << distances.sum() << '/' << distances.nb() << '\n'
            << "stddev: " << distances.stddev() << '\n'
            << "density: " << actual << (actual < expected ? " < " : " > ")  << expected << '\n'
            << '\t' << distances.nb() << '/' << seen_mers << ' ' << (seen_mers - w + 1) << '\n';

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
