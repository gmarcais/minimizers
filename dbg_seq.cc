#include <unistd.h>
#include <iostream>
#include "dbg_seq.hpp"


const int params[][2] = {
  {0, 0}, {0, 0},  {0, 0},
  {1, 0}, {1, 0},  {2, 0},  {1, 0}, {1, 0},
  {1, 5}, {4, 0},  {3, 0},  {2, 0}, {3, 4},
  {1, 3}, {1, 11}, {1, 0},  {2, 3}, {3, 0},
  {7, 0}, {1, 5},  {3, 0},  {2, 0}, {1, 0},
  {5, 0}, {1, 3},  {3, 0},  {1, 7}, {1, 7},
  {3, 0}, {2, 0},  {1, 15}, {3, 0}, {1, 27},
};

struct binary_out {
  int        m_wi;
  const int  m_n;
  const bool m_do_wrap;
  int*       m_wrap;


  binary_out(int n, bool do_wrap)
    : m_wi(0), m_n(n), m_do_wrap(do_wrap)
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

void out(int x, ssize_t& length) {
  if(length != 0) {
    std::cout << x;
    if(length > 0)
      --length;
  }
}

void seq_1(int n, int s, ssize_t length) {
  int x[n];
  int k = 0;
  int j = s;
  int r = 0;

  x[0] = 1;
  for(int i = 1; i < n; ++i)
    x[i] = 0;

  while(length != 0) {
    out(x[k], length);
    if(x[k])
      r = 0;
    else {
      ++r;
      if(r == n - 1) {
        out(0, length);
        break;
      }
    }
    k = ((k == 0) ? n : k) - 1;
    j = ((j == 0) ? n : j) - 1;
    x[k] ^= x[j];
  }

  if(length != 0)
    seq_1(n, s, length);
}

void seq_2(int n, int s, int t, ssize_t length) {

  int x[n];
  int k = 0;
  int j = s;
  int r = 0;
  int i = t;
  int h = s + t;

  x[0] = 1;
  for(int l = 1; l < n; ++l)
    x[l] = 0;

  while(length != 0) {
    out(x[k], length);

    if(x[k])
      r = 0;
    else {
      ++r;
      if(r == n - 1) {
        out(0, length);
        break;
      }
    }

    k = ((k == 0) ? n : k) - 1;
    j = ((j == 0) ? n : j) - 1;
    i = ((i == 0) ? n : i) - 1;
    h = ((h == 0) ? n : h) - 1;

    x[k] ^= x[j] ^ x[i] ^ x[h];
  }

  if(length != 0) // wrap around
    seq_2(n, s, t, length);
}

void seq(int n, ssize_t length) {
  if(n == 0) {
    // NOOP
  } else if(n == 1) {
    while(length != 0) {
      out(0, length);
      out(1, length);
    }
  } else if(n == 2) {
    while(length != 0) {
      out(0, length), out(0, length);
      out(1, length), out(1, length);
    }
  } else if(params[n][1] == 0) {
    seq_1(n, params[n][0], length);
  } else {
    seq_2(n, params[n][0], params[n][1], length);
  }
}

int main(int argc, char *argv[]) {
  std::ios::sync_with_stdio(false);
  dbg_seq args(argc, argv);

  if(args.order_arg >= (int)(sizeof(params) / sizeof(int[2])) || args.order_arg < 0) {
    dbg_seq::error() << "order must be in [0, " << (sizeof(params) / sizeof(int[2])) << "]";
  }

  ssize_t length = (ssize_t)1 << args.order_arg;
  if(args.wrap_flag)
    length += args.order_arg - 1;
  else if(args.infinity_flag)
    length = -1;
  else if(args.length_given)
    length = args.length_arg;

  if(args.fasta_flag)
    std::cout << ">dbg_binary order:" << args.order_arg << " length:" << length << '\n';
  seq(args.order_arg, length);
  if(isatty(1) || args.fasta_flag)
    std::cout << '\n';

  return 0;
}
