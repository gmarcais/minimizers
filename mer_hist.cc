#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <limits>

#include "misc.hpp"
#include "mer_hist.hpp"

int main(int argc, char* argv[]) {
  mer_hist args(argc, argv);

  initialize_codes(2);

  std::ifstream is(args.fasta_arg);
  if(!is.good())
    mer_hist::error() << "Error opening fasta file '" << args.fasta_arg << '\'';

  slide_mer<1>        mer(args.k_arg);
  std::vector<size_t> counts((size_t)1 << args.k_arg, 0);
  std::vector<size_t> histo;

  char c = is.peek();
  for( ; c != EOF; c = is.peek()) {
    if(c == '>') { // Skip header
      mer.truncate();
      is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      continue;
    }

    for(c = is.get(); c != EOF && c != '\n'; c = is.get()) {
      mer.appende(c);
      if(mer.full()) {
        auto nc = ++counts[mer.mer];
        if(nc > 1)
          --histo[nc - 1];
        if(nc >= histo.size())
          histo.resize(nc + 1, 0);
        ++histo[nc];
      }
    }
  }

  for(size_t i = 1; i < histo.size(); ++i)
    if(histo[i] > 0)
      std::cout << i << ' ' << histo[i] << '\n';

  return 0;
}
