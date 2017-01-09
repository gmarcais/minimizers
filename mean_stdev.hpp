#ifndef __MEAN_STDEV_H__
#define __MEAN_STDEV_H__

class mean_stdev {
  double              A = 0, Q = 0, S = 0;
  unsigned long       n = 0, c = 0;
  std::vector<size_t> h;

public:

  void sample(const unsigned int v) {
    const double y = v - A;
    S += v;
    ++n;
    A += y / n;
    Q += y * (v - A);
    if(v >= h.size())
      h.resize(v + 1, 0);
    ++h[v];
  }
  void count() { ++c; }

  unsigned long nb() const { return n; }
  unsigned long total() const { return c; }
  double mean() const {
    if(n == 0)
      throw std::domain_error("No samples");
    return A;
  }
  double variance() const {
    if(n < 2)
      throw std::domain_error("Not enough samples");
    return Q / (n - 1);
  }
  double sum() const { return S; }
  double stddev() const { return std::sqrt(variance()); }
  const std::vector<size_t>& histo() const { return h; }
};
inline mean_stdev& operator<<(mean_stdev& ms, const double v) {
  ms.sample(v);
  return ms;
}

#endif // __MEAN_STDEV_H__
