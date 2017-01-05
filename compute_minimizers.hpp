#ifndef __COMPUTE_MINIMIZERS_H__
#define __COMPUTE_MINIMIZERS_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "misc.hpp"

extern std::vector<uint64_t> order;

struct mer_pos {
  uint64_t      mer;
  size_t        pos;
  static size_t window;
  mer_pos() : mer(0), pos(0) { }
  mer_pos(uint64_t m, size_t p) : mer(m), pos(p) { }
  // bool operator<(const mer_pos& rhs) const {
  //   return pos + window < rhs.pos || mer < rhs.mer;
  // }
};

template<typename MerComp>
struct mer_pos_comp {
  MerComp comp;
  inline bool operator()(const mer_pos& x, const mer_pos& y) const { return comp(x.mer, y.mer); }
};

template<typename MerComp>
bool greater(const mer_pos& x, const mer_pos& y) {
  return x.pos + mer_pos::window < y.pos || MerComp()(x.mer, y.mer);
}

template<int BITS>
void print(std::ostream& os, const std::vector<mer_pos>& mers, int k) {
  os << '{';
  for(const auto& it : mers)
    os << ' ' << mer_to_string<BITS>(it.mer, k) << ':' << it.pos << ':' << order[it.mer];
  os << " }";
}

template<typename MerComp, int BITS>
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
    mean_stdev           ms;
    slide_mer<BITS>      mer(k);
    size_t pos_i             = 0; // Index in circular buffer;
    size_t               pos = 0;
    std::vector<mer_pos> mers(mer_pos::window); // Circular buffer
    mer_pos_comp<MerComp> comp;

    // Fill up first mer
    for(size_t i = 0; i < k && first != last; ++first, ++i, ++pos)
      mer.appende(*first);
    mers[0] = mer_pos(mer.mer, pos - k);

    // Fill up first window
    size_t min_pos_i = 0;
    for(pos_i = 1; first != last && pos_i < mer_pos::window; ++first, ++pos, ++pos_i) {
      mer.appende(*first);
      mers[pos_i] = mer_pos(mer.mer, pos - k + 1);
      if(comp(mers[pos_i], mers[min_pos_i]))
        min_pos_i = pos_i;
    }
    if(pos_i < mer_pos::window)
      return ms;

    size_t min_pos = min_pos_i;
    act(mer_pos(mers[min_pos_i]));
    ms.count();

    for(pos_i = 0; first != last; ++first, ++pos, pos_i = (pos_i + 1) % mer_pos::window) {
      ms.count();
      mer.appende(*first);
      mer_pos mp(mer.mer, pos - k + 1);
      if(comp(mp, mers[min_pos_i])) { // a new minimimum arrived
        mers[pos_i] = mp;
        act(mp);
        ms.sample(pos - min_pos);
        min_pos_i = pos_i;
        min_pos = pos;
      } else if(min_pos_i == pos_i) { // the current min fell out of the window
        mers[pos_i] = mp;
        min_pos_i = std::min_element(mers.cbegin(), mers.cend(), comp) - mers.cbegin();
        act(mers[min_pos_i]);
        ms.sample(pos - min_pos);
        min_pos = pos;
      } else {
        mers[pos_i] = mp;
      }
    }

    return ms;
  }
};


#endif /* __COMPUTE_MINIMIZERS_H__ */
