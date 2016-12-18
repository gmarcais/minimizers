#ifndef __COMPUTE_MINIMIZERS_H__
#define __COMPUTE_MINIMIZERS_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "misc.hpp"


struct mer_pos {
  uint64_t      mer;
  size_t        pos;
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
  template<typename String, typename Container>
  inline mean_stdev operator()(const String& seq, const size_t k,
                               Container& minimizers) {
    return this->operator()(seq.cbegin(), seq.cend(), k,
                            [&](const mer_pos& mp) { minimizers.insert(mp.mer); });
  }

  template<typename Iterator, typename Action>
  mean_stdev operator()(Iterator first, Iterator last, const size_t k,
                        Action act) {
    std::vector<mer_pos> heap;
    mean_stdev           ms;
    slide_mer            mer(k);
    size_t               pos     = 0;
    size_t               min_pos = 0; // Position of last minizer
    const size_t         base_w  = k + mer_pos::window - 1;

    for( ; first != last; ++first, ++pos) {
      if(pos >= base_w) {
        while(!heap.empty() && heap.front().pos + mer_pos::window < pos) {
          std::pop_heap(heap.begin(), heap.end(), greater<MerComp>);
          heap.pop_back();
        }
        if(!heap.empty() && heap.front().pos != min_pos) {
          act(heap.front());
          ms.sample(heap.front().pos - min_pos);
          min_pos = heap.front().pos;
        }
      }

      mer.append(*first);
      if(mer.full()) {
        heap.push_back({ mer.mer, pos });
        std::push_heap(heap.begin(), heap.end(), greater<MerComp>);
      }
    }
    //    std::cerr << '\n';

    return ms;
  }
};


#endif /* __COMPUTE_MINIMIZERS_H__ */
