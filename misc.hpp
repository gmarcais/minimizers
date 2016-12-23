#ifndef __MISC_H__
#define __MISC_H__

#include <stdexcept>
#include <cmath>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

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

// Mask
template<int BITS>
struct Mask {
  static const uint64_t base = (~(uint64_t)0) >> (sizeof(uint64_t) * 8 - BITS);
  static constexpr uint64_t bases(size_t s) { return (~(uint64_t)0) >> (sizeof(uint64_t) * 8 - BITS * s); }
};

std::string read_fasta(const char* path);

inline int base_to_code(char c) {
  switch(c) {
  case 'A': case 'a': case '0': return 0;
  case 'C': case 'c': case '1': return 1;
  case 'G': case 'g': case '2': return 2;
  case 'T': case 't': case '3': return 3;
  case '4': return 4;
  case '6': return 5;
  case '7': return 6;
  case '8': return 7;
  default: return -1;
  }
}

extern const char* conv; // == "ACGT"
template<int BITS>
std::string mer_to_string(uint64_t x, size_t k) {
  std::string res;
  for(int i = k - 1; i >= 0; i--)
    res += conv[(x >> (i * BITS)) & Mask<BITS>::base];
  return res;

}

template<int BITS>
uint64_t string_to_mer(const std::string& line, size_t k) {
  uint64_t m = 0;
  for(size_t i = 0; i < k; ++i) {
    int code = base_to_code(line[i]);
    if(code < 0) {
      std::cerr << "Invalid characters in mer " << line << '\n';
      exit(1);
    }
    m = (m << BITS) | (code & Mask<BITS>::base);
  }
  return m;
}

template<int BITS>
struct slide_mer {
  const uint64_t mask;
  const size_t   k;
  size_t         len = 0;
  uint64_t       mer = 0;

  slide_mer(size_t k_)
    : mask(Mask<BITS>::bases(k_))
    , k(k_)
  { }

  void append(const char c) {
    int code = base_to_code(c);
    if(code >= 0) {
      mer = ((mer << BITS) | code) & mask;
      len = std::min(len + 1, k);
    } else {
      len = 0;
    }
  }

  bool full() const { return len == k; }
};

extern const char* bconf; // == "01"
std::string bmer_to_string(uint64_t x, size_t k);
uint64_t string_to_bmer(const std::string& line, size_t k);

// Read mers from a file and add to a set
template<int BITS>
std::unordered_set<uint64_t> read_mers(const char* path, const size_t k) {
  std::unordered_set<uint64_t> res;
  std::string                  line;
  std::ifstream                is(path);
  size_t                       nb_mers = 0;

  while(std::getline(is, line)) {
    if(line.size() != k) {
      std::cerr << "All mers must have length " << k << ", got '" << line << "'\n";
      exit(1);
    }
    uint64_t m = string_to_mer<BITS>(line, k);
    auto tmp = res.insert(m);
    if(!tmp.second) {
      std::cerr << "Duplicated mer in set " << mer_to_string<BITS>(m, k) << ':' << line << '\n';
      exit(1);
    }
    ++nb_mers;
  }

  if(nb_mers != res.size()) {
    std::cerr << "Error, size of set != nb mers inserted\n";
    exit(1);
  }

  return res;
}

template<int BITS>
int compute_order(uint64_t mer, size_t k) {
  const int off  = BITS * (k - 1);
  uint64_t  rmer = mer;
  int       o    = 0;

  do {
    rmer = (rmer >> BITS) | ((rmer & Mask<BITS>::base) << off);
    ++o;
  } while(rmer ^ mer);
  return o;
}

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
