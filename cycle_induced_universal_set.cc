#include "misc.hpp"
#include "cycle_induced_universal_set.hpp"

inline uint64_t rotate_right(uint64_t m, unsigned int off) { return (m >> 2) | (((m & 0x3)) << off); }

// Return the distance to the minimum k-mer in the conjugacy class
unsigned int distance_to_min(uint64_t m, unsigned int k) {
  const unsigned int off   = 2 * (k - 1);
  uint64_t           min_m = m;
  unsigned int       dist  = 0;

  //  std::cout << "-> " << k << ' ' << mer_to_string(m, k) << '\n';
  unsigned int i = 1;
  for(uint64_t c_m = rotate_right(m, off); c_m != m; ++i, c_m = rotate_right(c_m, off)) {
    if(c_m < min_m) {
      dist  = i;
      min_m = c_m;
    }
    //    std::cout << i << ' ' << mer_to_string(c_m, k) << ' ' << mer_to_string(min_m, k) << ' ' << dist << '\n';
  }

  return dist;
}

int main(int argc, char* argv[]) {
  cycle_induced_universal_set args(argc, argv);

  // Compute number of mers in G_{k-1}
  const uint64_t nb_mers = 1 << (2 * (args.k_arg - 1));
  const uint64_t mask    = nb_mers - 1;

  std::vector<unsigned int> values(nb_mers);

  // Compute value for every mer in G_{k-1}
  for(uint64_t m = 0; m < nb_mers; ++m) {
    values[m] = distance_to_min(m, args.k_arg - 1) % args.w_arg;
    // std::cout << mer_to_string(m, args.k_arg - 1) << ' ' << values[m] << '\n';
  }

  // For every edge of G_{k-1} (i.e. node of G_k), output if edge is
  // select, that is if it goes from a high value to a no-less value.
  for(uint64_t m = 0; m < nb_mers; ++m) {
    for(uint64_t i = 0; i < 4; ++i) {
      const uint64_t nm = (m << 2) | i;
      // std::cout << mer_to_string(m, args.k_arg - 1) << '(' << values[m] << ") -> "
      //           << mer_to_string(nm & mask, args.k_arg - 1) << '(' << values[nm & mask] << ") "
      //           << m << '\n';
      if(values[nm & mask] <= values[m])
        std::cout << mer_to_string(nm, args.k_arg) << '\n';
    }
  }

  return 0;
}
