#include <cstdlib>
#include <unordered_set>
#include <vector>
#include <fstream>

#include "misc.hpp"
#include "longest_path.hpp"

longest_path args;

template<int BITS>
bool dfs_visit(uint64_t m, std::vector<size_t>& paths, std::vector<bool>& marks,
               const std::unordered_set<uint64_t>& mers, const uint64_t mask) {
  if(paths[m] > 0) return true;
  if(marks[m]) return false;

  marks[m] = true;

  size_t max_depth = 1;
  for(int i = 0; i < (1 << BITS); ++i) {
    const uint64_t nm = (mask & (m << BITS)) | i;
    if(mers.find(nm) != mers.end()) continue;
    if(!dfs_visit<BITS>(nm, paths, marks, mers, mask)) {
      std::cout << mer_to_string<BITS>(nm, args.k_arg) << '|';
      marks[m] = false;
      return false;
    }
    max_depth = std::max(max_depth, paths[nm] + 1);
  }

  paths[m] = max_depth;
  marks[m] = false;
  return true;
}

template<int BITS>
int compute_max_path() {
  const auto mers = read_mers<BITS>(args.vertices_arg, args.k_arg);
  const uint64_t mask    = Mask<BITS>::bases(args.k_arg);
  const uint64_t nb_mers = mask + 1;
  std::vector<size_t> paths(nb_mers, 0);
  std::vector<bool> marks(nb_mers, false);

  size_t max_path = 0;
  for(uint64_t m = 0; m < nb_mers; ++m) {
    if(!dfs_visit<BITS>(m, paths, marks, mers, mask)) {
      std::cerr << "\nCycle!" << std::endl;
      return 1;
    }
    max_path = std::max(max_path, paths[m]);
  }

  std::cout << max_path << '\n';
  return 0;
}

int main(int argc, char* argv[]) {
  args.parse(argc, argv);

    const int alphabet_size =
    args.size_given ? initialize_codes(args.size_arg) : initialize_codes(args.alphabet_arg);

  if(alphabet_size == 2)
    return compute_max_path<1>();
  else if(alphabet_size == 4)
    return compute_max_path<2>();
  else {
    std::cerr << "Alphabet size not supported" << std::endl;
    return 1;
  }

  return 0;
}
