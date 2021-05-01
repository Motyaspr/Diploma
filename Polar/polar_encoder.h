#ifndef DIPLOMA_POLAR_ENCODER_H
#define DIPLOMA_POLAR_ENCODER_H

std::vector<std::vector<bool>>
mul_matrixes(const std::vector<std::vector<bool>> &a, const std::vector<std::vector<bool>> &b) {
    std::vector<std::vector<bool>> ans(a.size(), std::vector<bool>(b[0].size(), 0));
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j < b[0].size(); j++) {
            for (int q = 0; q < a[0].size(); q++) {
                ans[i][j] = ((a[i][q] & b[q][j]) ^ ans[i][j]);
            }
        }
    }
    return ans;
}

std::vector<std::vector<bool>>
kron_mul(const std::vector<std::vector<bool>> &a, const std::vector<std::vector<bool>> &b) {
    std::vector<std::vector<bool>> ans(a.size() * b.size(), std::vector<bool>(a[0].size() * b[0].size(), 0));
    for (size_t i = 0; i < a.size(); i++) {
        for (size_t j = 0; j < a[0].size(); j++) {
            for (size_t x = 0; x < b.size(); x++)
                for (size_t y = 0; y < b[0].size(); y++)
                    ans[x * a.size() + i][y * a[0].size() + j] = (a[i][j] & b[x][y]);
        }
    }
    return ans;
}

struct PolarEncoder {
    int n, k;
    double sigma;
    std::vector<std::vector<bool>> g;
    std::vector<std::vector<bool>> A;
    std::vector<size_t> rev;
    std::vector<bool> frozen;

    PolarEncoder(int _n, int _k) {
        sigma = 0.0;
        n = _n;
        k = _k;
        frozen.resize(n, 0);
        rev = genReversedIndex(n);
        create_A();
        auto mat = createB(n);
        g = mul_matrixes(mat, A);

    }

    void reuse_frozen(double prob) {
        sigma = prob;
        std::fill(frozen.begin(), frozen.end(), false);
        std::vector<std::pair<double, size_t>> rates;
        for (size_t i = 0; i < n; i++)
            rates.emplace_back(-get_rate(n, i + 1), i);
        std::sort(rates.begin(), rates.end());
        for (size_t i = 0; i < n - k; i++)
            frozen[rates[i].second] = true;
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

    std::vector<std::vector<bool>> getRealGenMatrix() {
        std::vector<std::vector<bool>> ans;
        for (size_t i = 0; i < g.size(); i++) {
            if (frozen[i])
                ans.push_back(g[i]);
        }
        return ans;
    }

    std::vector<bool> encode(const std::vector<bool> &word, bool f) {
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
        if (f)
            return before_shuffle;
        return codeword;
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

};

#endif //DIPLOMA_POLAR_ENCODER_H
