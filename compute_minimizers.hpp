#ifndef __COMPUTE_MINIMIZERS_H__
#define __COMPUTE_MINIMIZERS_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "misc.hpp"


struct mer_pos {
  uint64_t mer;
  size_t   pos;
  static size_t window;
  // bool operator<(const mer_pos& rhs) const {
  //   return pos + window < rhs.pos || mer < rhs.mer;
  // }
};

template<typename MerComp>
bool greater(const mer_pos& x, const mer_pos& y) {
  return x.pos + mer_pos::window < y.pos || MerComp()(x.mer, y.mer);
}

template<typename MerComp>
struct compute_minimizers {
  mean_stdev operator()(const std::string& seq, const size_t k,
                        std::unordered_set<uint64_t>& minimizers) {
    std::vector<mer_pos> heap;
    mean_stdev           ms;
    slide_mer            mer(k);
    size_t               pos     = 0;
    size_t               min_pos = 0; // Position of last minizer
    const size_t         base_w  = k + mer_pos::window - 1;
    static const size_t  updates = 1024 * 1024;

    for(auto it = seq.cbegin(); it != seq.cend(); ++it, ++pos) {
      if(pos % (updates) == 0)
        std::cerr << '\r' << (seq.size() / updates) << ' ' << (pos / updates) << std::flush;
      if(pos >= base_w) {
        while(!heap.empty() && heap.front().pos + mer_pos::window < pos) {
          std::pop_heap(heap.begin(), heap.end(), greater<MerComp>);
          heap.pop_back();
        }
        if(!heap.empty() && heap.front().pos != min_pos) {
          //        std::cout << heap.front().pos << ' ' << pos << ' ' << mer_to_string(heap.front().mer, k) << ' ' << seq.substr(pos - base_w, base_w) <<  '\n';
          minimizers.insert(heap.front().mer);
          ms.sample(heap.front().pos - min_pos);
          min_pos = heap.front().pos;
        }
      }

      mer.append(*it);
      if(mer.full()) {
        heap.push_back({ mer.mer, pos });
        std::push_heap(heap.begin(), heap.end(), greater<MerComp>);
      }
    }
    std::cerr << '\n';

    return ms;
  }
};


#endif /* __COMPUTE_MINIMIZERS_H__ */
