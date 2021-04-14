#ifndef DIPLOMA_POLAR_ENCODER_H
#define DIPLOMA_POLAR_ENCODER_H

#include <common.h>

struct PolarEncoder {
    int n, k;
    double sigma;
    std::vector<size_t> rev;
    std::vector<bool> frozen;

    PolarEncoder(int _n, int _k, double _prob) {
        sigma = _prob;
        n = _n;
        k = _k;
        frozen.resize(n, 0);
        std::vector<std::pair<double, size_t>> rates;
        for (size_t i = 0; i < n; i++)
            rates.emplace_back(-get_rate(n, i + 1), i);
        std::sort(rates.begin(), rates.end());
        for (size_t i = 0; i < n - k; i++)
            frozen[rates[i].second] = true;
        rev = genReversedIndex(n);

    }

    double get_rate(size_t cur_n, size_t cur_i) {
        if (cur_n == 1 && cur_i == 1)
            return 2 / (sigma * sigma);
        if (cur_i % 2 == 1) {
            double rate = get_rate(cur_n / 2, (cur_i + 1) / 2);
            return 2 * rate;
        }
        double rate = get_rate(cur_n / 2, cur_i / 2);
        return bigXi(rate);
    }


    std::vector<bool> encode(const std::vector<bool> &word) {
        std::vector<bool> before_shuffle(n);
        int j = 0;
        for (size_t i = 0; i < n; i++)
            if (frozen[i])
                before_shuffle[i] = 0;
            else
                before_shuffle[i] = word[j++];
        for (unsigned long i = 2; i <= n; i *= 2)
            for (unsigned long q = 0; q < n; q += i) {
                unsigned long start = q + i - 1;
                unsigned long end = q + i / 2 - 1;
                for (unsigned long cur_k = start; cur_k > end; --cur_k) {
                    before_shuffle[cur_k - i / 2] = (before_shuffle[cur_k] ^ before_shuffle[cur_k - i / 2]);
                }
            }
        std::vector<bool> codeword(n, 0);
        for (size_t i = 0; i < rev.size(); i++)
            codeword[i] = before_shuffle[rev[i]];
        return before_shuffle;
    }

    double bigXi(double x) {
        if (x > 12)
            return 0.9861 * x - 2.3152;
        if (x > 3.5)
            return x * (9.005 * 0.001 * x + 0.7694) - 0.9507;
        if (x > 1)
            return x * (0.062883 * x + 0.3678) - 0.1627;
        return x * (0.2202 * x + 0.06448);
    }

    std::vector<size_t> genReversedIndex(size_t i) {
        std::vector<size_t> res;
        res.resize(i);
        res[0] = 0;
        for (size_t i1 = 1; i1 < i; i1 <<= 1u) {
            for (size_t j = 0; j < i1; ++j) {
                res[j] <<= 1u;
                res[j + i1] = res[j] + 1;
            }
        }
        return res;
    }


};

#endif //DIPLOMA_POLAR_ENCODER_H
