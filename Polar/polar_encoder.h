#ifndef DIPLOMA_POLAR_ENCODER_H
#define DIPLOMA_POLAR_ENCODER_H

#include <common.h>

struct PolarEncoder {
    std::vector<std::vector<bool>> g;
    std::vector<std::vector<bool>> A;
    int n, k;
    double sigma;
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
        for (size_t i = 0; i < frozen.size(); i++) {
            if (frozen[i])
                std::cout << i << ' ';
        }
        create_A();
        auto mat = createB(n);
        g = mul_matrixes(mat, A);
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

    void create_A() {
        std::vector<std::vector<bool>> x;
        x.resize(2, std::vector<bool>(2, true));
        x[0][1] = 0;
        int n_cur = n;
        A.resize(1, std::vector<bool>(1, 1));
        while (n_cur > 1) {
            auto ans = kron_mul(A, x);
            A = ans;
            n_cur /= 2;
        }
    }

    std::vector<std::vector<bool>> calc_r(std::vector<std::vector<bool>> x) {
        std::vector<std::vector<bool>> ans(x.size());
        for (size_t i = 0; i < ans.size(); i++) {
            ans[i].resize(x[i].size());
            for (size_t j = 0; j < x[i].size(); j++) {
                if (j % 2 == 0)
                    ans[i][j / 2] = x[i][j];
                else
                    ans[i][x[i].size() / 2 + j / 2] = x[i][j];
            }
        }
        return ans;
    }

    std::vector<std::vector<bool>> createB(int size) {
        std::vector<std::vector<bool>> I;
        I.resize(2, std::vector<bool>(2, false));
        I[0][0] = true;
        I[1][1] = true;
        std::vector<std::vector<bool>> ans = I;
        for (size_t q = 4; q <= size; q *= 2) {
            auto x = kron_mul(I, ans);
            ans = calc_r(x);
        }
        return ans;
    }

    std::vector<bool> encode(const std::vector<bool> &word) {
        std::vector<bool> ans;
        int j = 0;
        for (size_t i = 0; i < n; i++)
            if (frozen[i])
                ans.push_back(0);
            else
                ans.push_back(word[j++]);
        std::cout << "PASSED WORD:\n";
        for (size_t i = 0; i < ans.size(); i++)
            std::cout << ans[i] << ' ';
        std::cout << "\n";
        return mulVectorMatrix(ans, g);
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


};

#endif //DIPLOMA_POLAR_ENCODER_H
