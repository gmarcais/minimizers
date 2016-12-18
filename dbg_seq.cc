#include <iostream>
#include "dbg_seq.hpp"


const int params[][2] = {
  {0, 0}, {0, 0},  {0, 0},
  {1, 0}, {1, 0},  {2, 0},  {1, 0}, {1, 0},
  {1, 5}, {4, 0},  {3, 0},  {2, 0}, {3, 4},
  {1, 3}, {1, 11}, {1, 0},  {2, 3}, {3, 0},
  {7, 0}, {1, 5},  {3, 0},  {2, 0}, {1, 0},
  {5, 0}, {1, 3},  {3, 0},  {1, 7}, {1, 7},
  {3, 0}, {2, 0},  {3, 15}, {3, 0}, {1, 27},
};

struct binary_out {
  int        m_wi;
  const int  m_n;
  const bool m_do_wrap;
  int*       m_wrap;


  binary_out(int n, bool do_wrap) : m_wi(0), m_n(n), m_do_wrap(do_wrap)
  {
    m_wrap = new int[n-1];
  }
  ~binary_out() {
    if(m_do_wrap) {
      for(int i = 0; i < m_n - 1; ++i)
        std::cout << m_wrap[i];
    }
    delete [] m_wrap;
  }

  void operator()(int b) {
    std::cout << b;

    if(m_wi < m_n - 1) {
      m_wrap[m_wi] = b;
      ++m_wi;
    }
  }
};

template<typename Output>
void seq_1(int n, int s, bool do_wrap) {
  Output out(n, do_wrap);

  int x[n];
  int k = 0;
  int j = s;
  int r = 0;

  x[0] = 1;
  for(int i = 1; i < n; ++i)
    x[i] = 0;

  while(true) {
    out(x[k]);
    if(x[k])
      r = 0;
    else {
      ++r;
      if(r == n - 1) {
        out(0);
        break;
      }
    }
    k = ((k == 0) ? n : k) - 1;
    j = ((j == 0) ? n : j) - 1;
    x[k] ^= x[j];
  }
}

template<typename Output>
void seq_2(int n, int s, int t, bool do_wrap) {
  Output out(n, do_wrap);

  int x[n];
  int k = 0;
  int j = s;
  int r = 0;
  int i = t;
  int h = s + t;

  x[0] = 1;
  for(int l = 1; l < n; ++l)
    x[l] = 0;

  while(true) {
    out(x[k]);

    if(x[k])
      r = 0;
    else {
      ++r;
      if(r == n - 1) {
        std::cout << "0";
        break;
      }
    }

    k = ((k == 0) ? n : k) - 1;
    j = ((j == 0) ? n : j) - 1;
    i = ((i == 0) ? n : i) - 1;
    h = ((h == 0) ? n : h) - 1;

    x[k] ^= x[j] ^ x[i] ^ x[h];
  }
}

int main(int argc, char *argv[]) {
  std::ios::sync_with_stdio(false);
  dbg_seq args(argc, argv);

  if(args.order_arg >= sizeof(params) / sizeof(int[2]) || args.order_arg < 1) {
    dbg_seq::error() << "order must be in [1, " << (sizeof(params) / sizeof(int[2])) << "]";
  }

  if(args.order_arg == 1) {
    binary_out out(1, args.wrap_flag);
    out(0), out(1);
  } else if(args.order_arg == 2) {
    binary_out out(2, args.wrap_flag);
    out(0), out(0), out(1), out(1);
  } else {
    if(params[args.order_arg][1] == 0) {
      seq_1<binary_out>(args.order_arg, params[args.order_arg][0], args.wrap_flag);
    } else {
      seq_2<binary_out>(args.order_arg, params[args.order_arg][0], params[args.order_arg][1], args.wrap_flag);
    }
  }
  std::cout << '\n';

  return 0;
}
