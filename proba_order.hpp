#ifndef __PROBA_ORDER_H__
#define __PROBA_ORDER_H__

#include "marsaglia_binary_matrices.hpp"

// Ordering of k-mers based on repeated randomized ordering
class proba_order {
protected:
  struct mat_thresh {
    const MarsagliaMatrix matrix;
    const uint64_t        threshold;
  };
  std::vector<mat_thresh> m_matrices;

public:
  proba_order(const MarsagliaMatrix& m, uint64_t threshold)
  { append(m, threshold); };

  bool exists(const MarsagliaMatrix& m) const { // Does this matrix m exists in m_matrices?
    return std::any_of(m_matrices.cbegin(), m_matrices.cend(),
                       [&](const mat_thresh& mt) { return mt.matrix == m; });
  }
  void append(const MarsagliaMatrix& m, uint64_t threshold) {
    m_matrices.push_back({m, threshold});
  }

  // Return -1 if x is in a hitting set and not y, 1 if the other way
  // around, 0 if both x and y are in the same hitting or not hitting
  // set category.
  int compare_u(uint64_t x, uint64_t y) const {
    bool x_u = false, y_u = false; // Are x and y probably universal mers?
    for(const auto mt : m_matrices) {
      x_u = x_u || (mt.matrix * x < mt.threshold);
      y_u = y_u || (mt.matrix * y < mt.threshold);
    }
    if( x_u && !y_u) return -1;
    if(!x_u &&  y_u) return 1;
    return 0;
  }

  bool is_less(uint64_t x, uint64_t y) const {
    switch(compare_u(x, y)) {
    case -1: return true;
    case  1: return false;
    case  0: return (m_matrices.back().matrix) * x < (m_matrices.back().matrix * y);
    }
    return false; // Should never get there
  }

  friend std::ostream& operator<<(std::ostream& os, const proba_order& order);
};

inline std::ostream& operator<<(std::ostream& os, const proba_order& order) {
  os << '{';
  for(const auto m : order.m_matrices)
    os << " [" << m.matrix.param() << ',' << m.matrix.order() << ',' << m.threshold << ']';
  return os << " }";
}

#endif // __PROBA_ORDER_H__
