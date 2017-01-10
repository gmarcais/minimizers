#include <stdexcept>
#include "dbg_gen.hpp"
#include "misc.hpp"

void out(int x, ssize_t& length) {
  if(length != 0) {
    std::cout << conv[x];
    if(length > 0)
      --length;
  }
}

void gen(int t, int p, std::vector<int>& a, const int k, const int s, ssize_t& length) {
  if(t > k) {
    if(k % p == 0) {
      for(int i = 1; i <= p && length != 0; ++i)
        out(a[i], length);
    }
  } else {
    a[t] = a[t-p];
    gen(t + 1, p, a, k, s, length);
    if(length == 0) return;
    for(int i = a[t-p] + 1; i < s && length != 0; ++i) {
      a[t] = i;
      gen(t + 1, t, a, k, s, length);
    }
  }
}

void debruijn(const int k, const int s, ssize_t length) {
  std::vector<int> a(k + 1, 0);

  gen(1, 1, a, k, s, length);
  if(length != 0)
    debruijn(k, s, length);
}

ssize_t power(ssize_t x, ssize_t y) {
  ssize_t res = 1;
  for( ; y; y >>= 1, x *= x)
    if(y & 1)
      res *= x;
  return res;
}

int main(int argc, char* argv[]) {
  std::ios::sync_with_stdio(false);
  dbg_gen args(argc, argv);

  const int alphabet_size =
    args.size_given ? initialize_codes(args.size_arg) : initialize_codes(args.alphabet_arg);

  ssize_t length = power(alphabet_size, args.order_arg);
  if(args.wrap_flag)
    length += args.order_arg - 1;
  else if(args.wrapbases_given)
    length += args.wrapbases_arg;
  else if(args.infinity_flag)
    length = -1;
  else if(args.length_given)
    length = args.length_arg;

  if(args.fasta_flag)
    std::cout << ">dbg_binary order:" << args.order_arg << " length:" << length << '\n';

  debruijn(args.order_arg, alphabet_size, length);

  if(isatty(1) || args.fasta_flag)
    std::cout << '\n';

  return 0;
}
