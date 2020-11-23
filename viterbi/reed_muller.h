//
// Created by Matvey Sprikut on 23.11.2020.
//

#ifndef DIPLOMA_REED_MULLER_H
#define DIPLOMA_REED_MULLER_H

#include <vector>
#include <iostream>

bool get_bit(size_t x, size_t i, size_t n) {
    return ((x >> (n - 1 - i)) & 1);
}

struct ReedMuller {
    int r, m, n;
    std::vector<std::vector<bool>> generated;

    int calc_C_m_k(int cur) const {
        int mul = 1;
        for (int i = cur + 1; i <= m; i++)
            mul *= i;
        for (size_t i = 2; i <= m - cur; i++)
            mul /= i;
        return mul;
    }

    int calc_N() {
        int ans = 1;
        for (int i = 1; i <= r; i++)
            ans += calc_C_m_k(i);
        return ans;
    }



    std::vector<bool> get_ones() const {
        std::vector<bool> a = std::vector<bool>(1 << m, true);
        return a;
    }

    std::vector<std::vector<bool>> get_G1() const {
        std::vector<std::vector<bool>> t(m);
        for (size_t i = 0; i < m; ++i)
            t[i].resize(1 << m);
        for (size_t j = 0; j < (1 << m); j++)
            for (size_t i = 0; i < m; i++) {
                t[i][j] = get_bit(j, i, m);
            }
        return t;
    }

    std::vector<bool> mul(std::vector<bool> &a, std::vector<bool> &b) {
        std::vector<bool> ans;
        ans.reserve(a.size());
        for (size_t i = 0; i < a.size(); i++)
            ans.push_back(a[i] & b[i]);
        return ans;
    }

    void gen_rest(std::vector<std::vector<bool>> &t) {
        for (size_t i = 2; i <= r; i++)
            brute(t, i, std::vector<size_t>());
    }

    void brute(std::vector<std::vector<bool>> &t, size_t k, std::vector<size_t> cur) {
        if (cur.size() == k) {
            generated.push_back(mul(t[cur[0]], t[cur[1]]));
            return;
        }
        size_t start = (cur.size()) ? cur.back() : -1;
        for (size_t i = start + 1; i < m; i++) {
            cur.push_back(i);
            brute(t, k, cur);
            cur.pop_back();
        }
    }

    ReedMuller(int new_r, int new_m) {
        r = new_r;
        m = new_m;
        n = calc_N();
        std::cout << n << "\n";
        generated.push_back(get_ones());
        auto t = get_G1();
        for (const auto & i : t)
            generated.push_back(i);
        gen_rest(t);
    }

    void print() {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < (1 << m); j++)
                std::cout << generated[i][j] << ' ';
            std::cout << "\n";
        }
    }

};


#endif //DIPLOMA_REED_MULLER_H
