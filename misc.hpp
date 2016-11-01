#ifndef __MISC_H__
#define __MISC_H__

#include <stdexcept>
#include <cmath>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <functional>

class mean_stdev {
  double        A = 0, Q = 0;
  unsigned long n = 0;
public:

  void sample(const double v) {
    const double y = v - A;
    ++n;
    A += y / n;
    Q += y * (v - A);
  }

  unsigned long nb() const { return n; }
  double mean() const {
    if(n == 0)
      throw std::domain_error("No samples");
    return A;
  }
  double variance() const {
    if(n < 2)
      throw std::domain_error("Not enough samples");
    return Q / (n - 1);
  }
  double stddev() const { return std::sqrt(variance()); }
};
inline mean_stdev& operator<<(mean_stdev& ms, const double v) {
  ms.sample(v);
  return ms;
}

std::string read_fasta(const char* path);

inline int base_to_code(char c) {
  switch(c) {
  case 'A': case 'a': return 0;
  case 'C': case 'c': return 1;
  case 'G': case 'g': return 2;
  case 'T': case 't': return 3;
  default: return -1;
  }
}

extern const char* conv; // == "ACGT"
std::string mer_to_string(uint64_t x, size_t k);


struct slide_mer {
  const uint64_t mask;
  const size_t   k;
  size_t         len = 0;
  uint64_t       mer = 0;

  slide_mer(size_t k_)
    : mask(~(uint64_t)0 >> (sizeof(uint64_t) * 8 - 2 * k_))
    , k(k_)
  { }

  void append(const char c) {
    int code = base_to_code(c);
    if(code >= 0) {
      mer = ((mer << 2) | code) & mask;
      len = std::min(len + 1, k);
    } else {
      len = 0;
    }
  }

  bool full() const { return len == k; }
};

// Read mers from a file and add to a set
std::unordered_set<uint64_t> read_mers(const char* path, const size_t k);

int compute_order(uint64_t mer, size_t k);

// Seed a random generator
template <typename EngineT, std::size_t StateSize = EngineT::state_size>
void
seed_prg(EngineT& engine)
{
  using          engine_type    = typename EngineT::result_type;
  using          device_type    = std::random_device::result_type;
  using          seedseq_type   = std::seed_seq::result_type;
  constexpr auto bytes_needed   = StateSize * sizeof(engine_type);
  constexpr auto numbers_needed = (sizeof(device_type) < sizeof(seedseq_type))
    ? (bytes_needed / sizeof(device_type))
    : (bytes_needed / sizeof(seedseq_type));
  std::array<device_type, numbers_needed> numbers {};
  std::random_device rnddev {};
  std::generate(numbers.begin(), numbers.end(), std::ref(rnddev));
  std::seed_seq seedseq(numbers.cbegin(), numbers.cend());
  engine.seed(seedseq);
}

#endif /* __MISC_H__ */
