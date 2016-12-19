#ifndef __MARSAGLIA_BINARY_MATRICES_H__
#define __MARSAGLIA_BINARY_MATRICES_H__

#include <cassert>
#include <stdexcept>

// From Xorshift RNGs paper by George MARSAGLIA, the following
// functions/matrices are invertible matrices on [0,2^64). (They also
// generate large period RNG, but that is not the direct focus here).

class MarsagliaMatrix {
public:
  static const uint16_t nb_param = 275;
  static const uint16_t nb_order = 8;
  static const int parameters[nb_param][3];

  MarsagliaMatrix(int parameter, int order)
    : m_param(parameter % nb_param)
    , m_order(order % nb_order)
  { }

  bool operator==(const MarsagliaMatrix& rhs) const { return m_param == rhs.m_param && m_order == rhs.m_order; }
  bool operator!=(const MarsagliaMatrix& rhs) const { return !(*this == rhs); }

  uint64_t operator*(uint64_t y) const {
    const int* const p = parameters[m_param];

    for(int i = 0; i < 20; ++i) {
      switch(m_order) {
      case 0: y^=y<<p[0]; y^=y>>p[1]; y^=y<<p[2]; break;
      case 1: y^=y<<p[2]; y^=y>>p[1]; y^=y<<p[0]; break;
      case 2: y^=y>>p[0]; y^=y<<p[1]; y^=y>>p[2]; break;
      case 3: y^=y>>p[2]; y^=y<<p[1]; y^=y>>p[0]; break;
      case 4: y^=y<<p[0]; y^=y<<p[2]; y^=y>>p[1]; break;
      case 5: y^=y<<p[2]; y^=y<<p[0]; y^=y>>p[1]; break;
      case 6: y^=y>>p[0]; y^=y>>p[2]; y^=y<<p[1]; break;
      case 7: y^=y>>p[2]; y^=y>>p[0]; y^=y<<p[1]; break;
      }
    }
    return y;
  }

  int16_t param() const { return m_param; }
  int16_t order() const { return m_order; }

protected:
  const int16_t m_param;
  const int16_t m_order;
};

#endif // __MARSAGLIA_BINARY_MATRICES_H__
