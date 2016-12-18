#include "misc.hpp"

#include <string>
#include <fstream>
#include <iostream>


std::string read_fasta(const char* path) {
  std::string   res, line;
  std::ifstream is(path);

  int c = is.peek();
  if(c != '>') {
    std::cerr << "Invalid fasta file" << std::endl;
    exit(1);
  }

  for( ; c != EOF; ) {
    std::getline(is, line);
    if(!res.empty())
      res += 'N';
    for(c = is.peek(); c != '>' && c != EOF; c = is.peek()) {
      std::getline(is, line);
      res += line;
    }
  }

  return res;
}

const char* conv = "ACGT";
std::string mer_to_string(uint64_t x, size_t k) {

  std::string res;
  for(int i = k - 1; i >= 0; i--)
    res += conv[(x >> (2 * i)) & 0x3];
  return res;
}

uint64_t string_to_mer(const std::string& line, size_t k) {
  uint64_t m = 0;
  for(size_t i = 0; i < k; ++i) {
    int code = base_to_code(line[i]);
    if(code < 0) {
      std::cerr << "Invalid characters in mer " << line << '\n';
      exit(1);
    }
    m = (m << 2) | code;
  }
  return m;
}

const char* bconv = "01";
std::string bmer_to_string(uint64_t x, size_t k) {
  std::string res;
  for(int i = k - 1; i >= 0; i--)
    res += bconv[(x >> i) & 0x1];
  return res;
}

int num_to_code(char c) {
  switch(c) {
  case '0': return 0;
  case '1': return 1;
  default: return -1;
  }
}
uint64_t string_to_bmer(const std::string& line, size_t k) {
  uint64_t m = 0;
  for(size_t i = 0; i < k; ++i) {
    int code = num_to_code(line[i]);
    if(code < 0) {
      std::cerr << "Invalid characters in mer " << line << '\n';
      exit(1);
    }
    m = (m << 1) | code;
  }
  return m;
}

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
    uint64_t m = string_to_mer(line, k);
    auto tmp = res.insert(m);
    if(!tmp.second) {
      std::cerr << "Duplicated mer in set " << mer_to_string(m, k) << ':' << line << '\n';
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

int compute_order(uint64_t mer, size_t k) {
  const int off  = 2 * (k - 1);
  uint64_t  rmer = mer;
  int       o    = 0;

  do {
    rmer = (rmer >> 2) | ((rmer & 0x3) << off);
    ++o;
  } while(rmer ^ mer);
  return o;
}
