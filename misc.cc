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

const char* conv = "012345678";

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
