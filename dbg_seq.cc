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

void seq_1(int n, int s, bool do_wrap) {
  int x[n];
  int wrap[n-1];
  int k = 0;
  int j = s;
  int r = 0;

  x[0] = 1;
  for(int i = 1; i < n; ++i)
    x[i] = 0;

  size_t wi = 0;
  while(true) {
    std::cout << x[k];
    if(wi < n - 1) {
      wrap[wi] = x[k];
      ++wi;
    }

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
    x[k] ^= x[j];
  }

  if(do_wrap) {
    for(int i = 0; i < n - 1; ++i)
      std::cout << wrap[i];
  }
  std::cout << '\n';
}

void seq_2(int n, int s, int t, bool do_wrap) {
  int x[n];
  int wrap[n-1];
  int k = 0;
  int j = s;
  int r = 0;
  int i = t;
  int h = s + t;

  x[0] = 1;
  for(int l = 1; l < n; ++l)
    x[l] = 0;

  size_t wi = 0;
  while(true) {
    std::cout << x[k];
    if(wi < n - 1) {
      wrap[wi] = x[k];
      ++wi;
    }

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
  if(do_wrap) {
    for(int i = 0; i < n - 1; ++i)
      std::cout << wrap[i];
  }
  std::cout << '\n';

}

int main(int argc, char *argv[]) {
  dbg_seq args(argc, argv);

  if(args.order_arg >= sizeof(params) / sizeof(int[2]) || args.order_arg < 1) {
    dbg_seq::error() << "order must be in [1, " << (sizeof(params) / sizeof(int[2])) << "]";
  }

  if(args.order_arg == 1) {
    std::cout << "01\n";
  } else if(args.order_arg == 2) {
    std::cout << "0011";
    if(args.wrap_flag)
      std::cout << '0';
    std::cout << '\n';
  } else {
    if(params[args.order_arg][1] == 0) {
      seq_1(args.order_arg, params[args.order_arg][0], args.wrap_flag);
    } else {
      seq_2(args.order_arg, params[args.order_arg][0], params[args.order_arg][1], args.wrap_flag);
    }
  }

  return 0;
}
