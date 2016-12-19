#include <iostream>
#include <random>

#include "minimizers_large.hpp"
#include "misc.hpp"
#include "compute_minimizers.hpp"
#include "proba_order.hpp"

const MarsagliaMatrix *order_mat = nullptr;
struct order_mat_mult {
  bool operator()(uint64_t x, uint64_t y) const { return (*order_mat * x) < (*order_mat * y); }
};

const proba_order* order_p = nullptr;
struct order_proba_order {
  bool operator()(uint64_t x, uint64_t y) const {
    switch(order_p->compare_u(x, y)) {
    case -1: return true;
    case 1: return false;
    case 0: return (*order_mat * x) < (*order_mat * y);
    }
    return false; // Should never get there
  }
};

template<typename PRG>
void generate_sequence(std::string& seq, PRG& prg) {
  std::uniform_int_distribution<int> base_dist(0, 3);
  const char* bases = "ACGT";
  for(size_t i = 0; i < seq.size(); ++i)
    seq[i] = bases[base_dist(prg)];
}

template<typename PRG>
MarsagliaMatrix new_matrix(const proba_order& o, PRG& prg) {
  std::uniform_int_distribution<int> param_dist(0, MarsagliaMatrix::nb_param - 1);
  std::uniform_int_distribution<int> order_dist(0, MarsagliaMatrix::nb_order - 1);

  MarsagliaMatrix res(param_dist(prg), order_dist(prg));
  return o.exists(res) ? new_matrix(o, prg) : res;
}

void print_used(const std::vector<bool>& used, uint64_t max, const MarsagliaMatrix& mat, int k) {
  std::cerr << "--------------------\n";
  for(uint64_t m = 0; m < used.size(); ++m) {
    std::cerr << mer_to_string(m, k) << ' ' << used[m] << ' ' << max << ' ' << (mat * m) << '\n';
  }
}

int main(int argc, char* argv[]) {
  const minimizers_large args(argc, argv);
  std::mt19937_64 prg;
  seed_prg(prg);

  const int      k       = args.k_arg;
  const uint64_t nb_mers = (uint64_t)1 << (2 * k);
  mer_pos::window        = args.w_arg;
  const int      w       = mer_pos::window;

  std::uniform_int_distribution<int> param_dist(0, MarsagliaMatrix::nb_param - 1);
  std::uniform_int_distribution<int> order_dist(0, MarsagliaMatrix::nb_order - 1);

  std::string sequence(args.size_arg, '\0');
  //  std::vector<bool> used(nb_mers, false);

  // First iteration
  generate_sequence(sequence, prg);

  const MarsagliaMatrix first_matrix(param_dist(prg), order_dist(prg));
  uint64_t max = 0;
  order_mat = &first_matrix;
  auto ms = compute_minimizers<order_mat_mult>()(sequence.cbegin(), sequence.cend(), args.k_arg,
                                                 [&](mer_pos& mp) { max = std::max(max, first_matrix * mp.mer);
                                                                    //                                                                    used[mp.mer] = true;
                                                 });
  //  print_used(used, max, first_matrix, k);

  const double expected = ((double)(nb_mers - 1) / (double)(nb_mers + w)) * (2.0 / (w + 1));
  double       actual   = (double)ms.nb() / (double)(sequence.size() - k + 1);
  std::cout << "max: " << (100 * (double)max / nb_mers) << max << '\n'
            << "mean: " << ms.mean() << '\n'
            << "stddev: " << ms.stddev() << '\n'
            << "density: " << actual << (actual < expected ? " < " : " > ")  << expected << '\n';

  proba_order order(first_matrix, max);
  order_p = &order;

  for(uint32_t i = 1; i < args.nb_iterations_arg; ++i) {
    generate_sequence(sequence, prg);
    const MarsagliaMatrix matrix = new_matrix(order, prg);
    order_mat                    = &matrix;
    max                          = 0;
    //    std::fill(used.begin(), used.end(), false);
    auto ms = compute_minimizers<order_proba_order>()(sequence.cbegin(), sequence.cend(), args.k_arg,
                                                      [&](mer_pos& mp) { max = std::max(max, matrix * mp.mer);
                                                                         //                                                                         used[mp.mer] = false;
                                                      });
    order.append(matrix, max);

    double       actual   = (double)ms.nb() / (double)(sequence.size() - k + 1);
    std::cout << "order: " << order << '\n'
              << "mean: " << ms.mean() << '\n'
              << "stddev: " << ms.stddev() << '\n'
              << "density: " << actual << (actual < expected ? " < " : " > ")  << expected << '\n';
  }

  return 0;
}
