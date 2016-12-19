#include <cstdlib>
#include <unordered_set>
#include <vector>
#include <fstream>

#include "misc.hpp"
#include "binary_longest_path.hpp"

binary_longest_path args;

std::unordered_set<uint64_t> read_mers(const char* path, unsigned int k) {
  std::unordered_set<uint64_t> res;
  std::ifstream is(path);
  std::string mer;

  if(!is.good()) {
    std::cerr << "Failed to open " << path << std::endl;
    exit(1);
  }
  while(is >> mer) {
    uint64_t m = string_to_bmer(mer, k);
    res.insert(m);
  }

  return res;
}

bool dfs_visit(uint64_t m, std::vector<size_t>& paths, std::vector<bool>& marks,
               const std::unordered_set<uint64_t>& mers, const uint64_t mask) {
  if(paths[m] > 0) return true;
  if(marks[m]) return false;

  marks[m] = true;

  size_t max_depth = 1;
  for(int i = 0; i < 2; ++i) {
    const uint64_t nm = (mask & (m << 1)) | i;
    if(mers.find(nm) != mers.end()) continue;
    if(!dfs_visit(nm, paths, marks, mers, mask)) {
      std::cout << bmer_to_string(nm, args.k_arg) << '|';
      marks[m] = false;
      return false;
    }
    max_depth = std::max(max_depth, paths[nm] + 1);
  }

  paths[m] = max_depth;
  marks[m] = false;
  return true;
}

int main(int argc, char* argv[]) {
  args.parse(argc, argv);

  const auto mers = read_mers(args.vertices_arg, args.k_arg);
  const uint64_t nb_mers = 1 << args.k_arg;
  const uint64_t mask    = nb_mers - 1;
  std::vector<size_t> paths(nb_mers, 0);
  std::vector<bool> marks(nb_mers, false);

  size_t max_path = 0;
  for(uint64_t m = 0; m < nb_mers; ++m) {
    if(!dfs_visit(m, paths, marks, mers, mask)) {
      std::cerr << "\nCycle!" << std::endl;
      exit(1);
    }
    max_path = std::max(max_path, paths[m]);
  }

  std::cout << max_path << '\n';

  return 0;
}
