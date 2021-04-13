#ifndef DIPLOMA_POLAR_ENCODER_H
#define DIPLOMA_POLAR_ENCODER_H

#include "common.h"

struct PolarEncoder {
    std::vector<std::vector<bool>> g;
    std::vector<std::vector<bool>> t;
    int n, k;
    double prob;
    std::vector<bool> frozen;

    PolarEncoder(int _n, int _k, double _prob) {
        prob = _prob;
        n = _n;
        k = _k;
        frozen.resize(n, 0);
        std::vector<std::pair<double, size_t>> rates;
        for (size_t i = 0; i < n; i++)
            rates.emplace_back(-get_rate(n, i + 1), i);
        std::sort(rates.begin(), rates.end());
        for (size_t i = 0; i < n - k; i++)
            frozen[rates[i].second] = true;
        create_A();
        print(t);

    }

    double get_rate(size_t cur_n, size_t cur_i) {
        if (cur_n == 1 && cur_i == 1)
            return prob;
        if (cur_i % 2 == 1) {
            double rate = get_rate(cur_n / 2, (cur_i + 1) / 2);
            return rate * 2 - rate * rate;
        }
        double rate = get_rate(cur_n / 2, cur_i / 2);
        return rate * rate;
    }

    void create_A() {
        std::vector<std::vector<bool>> x;
        x.resize(2, std::vector<bool>(2, true));
        x[0][1] = 0;
        int n_cur = n;
        t.resize(1, std::vector<bool>(1, 1));
        while (n_cur > 1) {
            auto ans = kron_mul(t, x);
            t = ans;
            n_cur /= 2;
        }

    }
};

#endif //DIPLOMA_POLAR_ENCODER_H
