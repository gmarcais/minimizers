#ifndef __COMPUTE_MINIMIZERS_H__
#define __COMPUTE_MINIMIZERS_H__

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include "misc.hpp"

extern std::vector<uint64_t> order;

enum Reason { NotMin, DropOff, NewMin };
struct mer_pos {
  uint64_t      mer;
  size_t        pos, win_start, win_pos;
  Reason        reason;
  static size_t window;
  mer_pos() : mer(0), pos(0) { }
  mer_pos(uint64_t m, size_t p, size_t wp) : mer(m), pos(p), win_start(wp), win_pos(0), reason(NotMin) { }
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

template<typename T>
void shift(std::vector<T>& a) {
  auto it = a.begin();
  auto nit = it + 1;
  const auto end = a.end();

  for( ; nit != end; ++it, ++nit)
    *it = *nit;
}

template<typename MerComp, int BITS>
struct compute_minimizers {
  template<typename String, typename Container>
  inline size_t operator()(const String& seq, const size_t k,
                               Container& minimizers) {
    return this->operator()(seq.cbegin(), seq.cend(), k,
                            [&](const mer_pos& mp) { minimizers.insert(mp.mer); });
  }

  template<typename Iterator, typename Action>
  size_t operator()(Iterator first, Iterator last, const size_t k,
                        Action act) {
    slide_mer<BITS>       mer(k);
    size_t                seen_mers = 0;
    size_t pos_i                  = 0; // Index in circular buffer;
    size_t                pos     = 0;
    std::vector<mer_pos>  mers(mer_pos::window); // Circular buffer
    mer_pos_comp<MerComp> comp;
    auto                  left_comp = [&comp](const mer_pos& a, const mer_pos& b) -> bool {
      return comp(a, b) || (!comp(b, a) && a.pos < b.pos);
    };
    std::vector<int> tie_vec(mer_pos::window, 0);

    // Fill up first mer
    for(size_t i = 0; i < k - 1 && first != last; ++first, ++i, ++pos)
      mer.appende(*first);

    // Fill up first window
    size_t min_pos_i = 0;
    for(pos_i = 0; first != last && pos_i < mer_pos::window - 1; ++first, ++pos, ++pos_i) {
      mer.appende(*first);
      ++seen_mers;
      mers[pos_i] = mer_pos(mer.mer, pos - k + 1, 0);
      if(comp(mers[pos_i], mers[min_pos_i]))
        min_pos_i = pos_i;
    }

    act(mers[min_pos_i]);

    for( ; first != last; ++first, ++pos, pos_i = (pos_i + 1) % mer_pos::window) {
      mer.appende(*first);
      ++seen_mers;
      mer_pos mp(mer.mer, pos - k + 1, pos + 2 - k - mer_pos::window);
      if(comp(mp, mers[min_pos_i])) { // a new minimimum arrived
        mp.reason = NewMin;
        mp.win_pos = mer_pos::window - 1;
        mers[pos_i] = mp;
        act(mp);
        min_pos_i = pos_i;
      } else if(min_pos_i == pos_i) { // || !comp(mers[min_pos_i], mp)) {
        // the current min fell out of the windowk
        // Don't care about new tie as we are using left most policy
        mp.reason = DropOff;
        mers[pos_i] = mp;
        const auto min_it = std::min_element(mers.cbegin(), mers.cend(), left_comp);

        min_pos_i = min_it - mers.cbegin();

        mers[pos_i].win_pos = (min_pos_i + (min_pos_i > pos_i ? 0 : mer_pos::window)) - pos_i - 1;
        act(mers[min_pos_i]);
      } else {
        mers[pos_i] = mp;
      }

      // Debug...
      // for(int ii = 0; ii < mer_pos::window; ++ii) {
      //   const int pi = (pos_i + 1 + ii) % mer_pos::window;
      //   const auto& mp = mers[pi];
      //   std::cerr << mer_to_string<BITS>(mp.mer, k) << (pi == min_pos_i ? '*' : ':') << mp.pos << ' ';
      // }
      // std::cerr << '\n';
    }
    return seen_mers;
  }

  template<typename Iterator, typename Action>
  void calc2(Iterator first, Iterator last, const size_t k, Action act,
                   bool no_first = false) {
    slide_mer<BITS>       mer(k);
    std::vector<mer_pos>  mers(mer_pos::window);
    size_t                min_pos = mer_pos::window + 1;
    mer_pos_comp<MerComp> comp;
    size_t                pos     = 0;

     // Fill up first mer
    for(size_t i = 0; i < k - 1 && first != last; ++first, ++i, ++pos)
      mer.appende(*first);

    // Fill up first window
    for(size_t pos_i = 0; first != last && pos_i < mer_pos::window - 1; ++first, ++pos, ++pos_i) {
      mer.appende(*first);
      mers[pos_i] = mer_pos(mer.mer, pos - k + 1);
    }

    if(no_first)
      min_pos = std::min_element(mers.begin(), mers.end() - 1, comp)->pos;

    for( ; first != last; ++first, ++pos, shift(mers)) {
      mer.appende(*first);
      mers[mer_pos::window - 1] = mer_pos(mer.mer, pos - k + 1);
      auto new_min = std::min_element(mers.begin(), mers.end(), comp);
      if(new_min->pos != min_pos) {
        new_min->win_pos = new_min - mers.begin();
        act(*new_min);
        min_pos = new_min->pos;
      }
      // for(const auto& mp : mers)
      //   std::cerr << mer_to_string<BITS>(mp.mer, k) << (mp.pos == min_pos ? '*' : ':') << mp.pos << ' ';
      // std::cerr << '\n';
    }
  }
};


#endif /* __COMPUTE_MINIMIZERS_H__ */
