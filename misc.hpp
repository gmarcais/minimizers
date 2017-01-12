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
#include <iterator>

// Mask
template<int BITS>
struct Mask {
  static const uint64_t base = (~(uint64_t)0) >> (sizeof(uint64_t) * 8 - BITS);
  static constexpr uint64_t bases(size_t s) { return (~(uint64_t)0) >> (sizeof(uint64_t) * 8 - BITS * s); }
};

std::string read_fasta(const char* path);

extern int ascii_to_code[256];
extern size_t mer_len;

// Only works if BITS is a power of 2
template<int BITS>
uint64_t reverse_bits(uint64_t x) {
  switch(BITS) {
  case  1: x = ((0x5555555555555555ULL & x) <<  1) | ((0xaaaaaaaaaaaaaaaaULL & x) >>  1);
  case  2: x = ((0x3333333333333333ULL & x) <<  2) | ((0xccccccccccccccccULL & x) >>  2);
  case  4: x = (( 0xf0f0f0f0f0f0f0fULL & x) <<  4) | ((0xf0f0f0f0f0f0f0f0ULL & x) >>  4);
  case  8: x = ((  0xff00ff00ff00ffULL & x) <<  8) | ((0xff00ff00ff00ff00ULL & x) >>  8);
  case 16: x = ((    0xffff0000ffffULL & x) << 16) | ((0xffff0000ffff0000ULL & x) >> 16);
  case 32: x = ((                        x) << 32) | ((                        x) >> 32);
  }
  return x;
}
template<int BITS>
inline uint64_t rev(uint64_t x) {
  return reverse_bits<BITS>(x) >> (8 * sizeof(uint64_t) - BITS * mer_len);
}
template<int BITS>
inline uint64_t rev_comp(uint64_t x) {
  return rev<BITS>(~(uint64_t)0 - x);
}
template<int BITS>
uint64_t canonical(uint64_t x) {
  uint64_t rcx = rev_comp<BITS>(x);
  return rcx < x ? rcx : x;
}


// Initialize ascii_to_code and return the size of the alphabet
int initialize_codes(const char* s);
int initialize_codes(int s);

inline int base_to_code(char c) {
  return ascii_to_code[(unsigned char)c];
}

extern const char* conv; // == "ACGT"

template<int BITS>
struct mer_to_string {
  uint64_t x;
  uint64_t k;
  explicit mer_to_string(uint64_t x_, uint64_t k_ = mer_len) : x(x_), k(k_) { }

  operator std::string() const {
    std::string res;
    for(int i = k - 1; i >= 0; i--)
      res += conv[(x >> (i * BITS)) & Mask<BITS>::base];
    return res;
  }
};

template<int BITS>
std::ostream& operator<<(std::ostream& os, const mer_to_string<BITS>& m) {
  for(int i = m.k - 1; i >= 0; i--)
    os<< conv[(m.x >> (i * BITS)) & Mask<BITS>::base];
  return os;
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

  void truncate() { len = 0; }

  void append(const char c) {
    int code = base_to_code(c);
    if(code >= 0) {
      mer = ((mer << BITS) | code) & mask;
      len = std::min(len + 1, k);
    } else {
      len = 0;
    }
  }

  void appende(const char c) {
    int code = base_to_code(c);
    if(code < 0)
      throw std::runtime_error(std::string("Unexpected character '") + c + '\'');
    mer = ((mer << BITS) | code) & mask;
    len = std::min(len + 1, k);
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

  if(!is.good())
    throw std::runtime_error(std::string("Failed to open file '") + path + "\'");

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
seed_prg(EngineT& engine, const char* save = nullptr, const char* load = nullptr)
{
  using          engine_type    = typename EngineT::result_type;
  using          device_type    = std::random_device::result_type;
  using          seedseq_type   = std::seed_seq::result_type;
  constexpr auto bytes_needed   = StateSize * sizeof(engine_type);
  constexpr auto numbers_needed = (sizeof(device_type) < sizeof(seedseq_type))
    ? (bytes_needed / sizeof(device_type))
    : (bytes_needed / sizeof(seedseq_type));
  std::array<device_type, numbers_needed> numbers {};

  if(load) {
    std::ifstream is(load);
    size_t i = 0;
    for( ; is && i < numbers_needed; ++i)
      is >> numbers[i];
    if(i != numbers_needed)
      std::runtime_error(std::string("Failed loading seed from '") + load + "'");
  } else {
    std::random_device rnddev {};
    std::generate(numbers.begin(), numbers.end(), std::ref(rnddev));
  }
  std::seed_seq seedseq(numbers.cbegin(), numbers.cend());
  engine.seed(seedseq);

  if(save) {
    std::ofstream os(save);
    for(size_t i = 0; i < numbers_needed; ++i)
      os << numbers[i] << '\n';
    if(!os.good())
      throw std::runtime_error(std::string("Failed writing seed to '") + save + "'");
  }
}

// Fasta sequence input iterator
class fasta_iterator : public std::iterator<std::input_iterator_tag, char> {
  std::filebuf* m_fb;
  std::istream* m_is;
  char          m_val;
  char          m_separator;

public:
  fasta_iterator(char sep = 'N')
    : m_fb(nullptr)
    , m_is(nullptr)
    , m_separator(sep)
  { }
  explicit fasta_iterator(std::istream& is, char sep = 'N')
    : m_fb(nullptr)
    , m_is(new std::istream(is.rdbuf()))
    , m_separator(sep)
  {
    skip_headers();
    ++*this;
  }
  explicit fasta_iterator(const char* path, char sep = 'N')
    : m_fb(new std::filebuf)
    , m_is(new std::istream(m_fb))
    , m_separator(sep)
  {
    if(!m_fb->open(path, std::ios::in))
      throw std::runtime_error(std::string("Failed to open fasta file '") + path + "'");
    skip_headers();
    ++*this;
  }
  explicit fasta_iterator(const std::string& path) : fasta_iterator(path.c_str()) { }
  fasta_iterator(const fasta_iterator& rhs)
    : m_fb(nullptr)
    , m_is(rhs.m_is ? new std::istream(rhs.m_is->rdbuf()) : nullptr)
    , m_val(rhs.m_val)
    , m_separator(rhs.m_separator)
  { }

  ~fasta_iterator() {
    delete m_is;
    delete m_fb;
  }

  bool operator==(const fasta_iterator& rhs) const { return m_is == rhs.m_is; }
  bool operator!=(const fasta_iterator& rhs) const { return m_is != rhs.m_is; }

  fasta_iterator& operator++() {
    bool m_nl = false;
    for(m_val = m_is->get(); m_val == '\n'; m_val = m_is->get())
      m_nl = true;
    if(m_nl && m_val == '>') {
      skip_headers();
      m_val = m_separator;
    }
    if(m_is->eof()) {
      delete m_is;
      m_is = nullptr;
    }
    return *this;
  }
  fasta_iterator operator++(int) {
    fasta_iterator res(*this);
    ++*this;
    return res;
  }
  char operator*() const { return m_val; }

protected:
  void skip_headers() {
    while(m_is->peek() == '>')
      m_is->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
};

inline fasta_iterator begin(fasta_iterator& it) { return it; }
inline fasta_iterator end(fasta_iterator& it) { return fasta_iterator(); }

#endif /* __MISC_H__ */
