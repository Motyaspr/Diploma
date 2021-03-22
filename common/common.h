//
// Created by Matvey Sprikut on 20.03.2021.
//

#ifndef DIPLOMA_COMMON_H
#define DIPLOMA_COMMON_H

#include <utility>
#include <vector>
#include <iostream>
#include <set>
#include <algorithm>
#include <iomanip>
#include <map>
#include <random>
#include <cmath>
#include "reed_muller.h"

const int INF = 10000;

const int ITER = 100;


size_t get_first_one(std::vector<bool> t) {
    for (size_t i = 0; i < t.size(); i++)
        if (t[i])
            return i;
    return INF;
}

struct matrix {
    size_t n{}, m{};
    std::vector<std::vector<bool>> mt;
    std::vector<size_t> ones;

    matrix() {
        n = 0, m = 0;
        mt.resize(n);
        ones.resize(n);
    }

    matrix(size_t nn, size_t mm) {
        n = nn;
        m = mm;
        mt.resize(n);
        for (size_t i = 0; i < n; i++) {
            mt[i].resize(m, false);
        }
    }

    matrix(std::vector<std::vector<bool>> &mat) {
        n = mat.size();
        ones.resize(n);
        for (size_t i = 0; i < n; i++) {
            m = mat[0].size();
            mt.push_back(mat[i]);
        }
    }

    static bool span_cmp(std::vector<bool> a, std::vector<bool> b) {
        return get_first_one(a) < get_first_one(b);
    }

    void change_first() {
        for (size_t iter = 0; iter < m; iter++) {
            sort(mt.begin(), mt.end(), span_cmp);
            int cur_ind = -1;
            int cur_row = -1;
            for (size_t i = 0; i < n; i++) {
                size_t v = get_first_one(mt[i]);
                if (v == cur_ind)
                    for (size_t j = 0; j < m; j++)
                        mt[i][j] = (mt[cur_row][j] ^ mt[i][j]);
                else {
                    cur_ind = get_first_one(mt[i]);
                    cur_row = i;
                }
            }

        }
        sort(mt.begin(), mt.end(), span_cmp);
    }

    void change_last() {
        std::vector<std::vector<size_t>> ends(m);
        for (size_t i = 0; i < n; i++)
            for (size_t j = m - 1; j >= 0; j--)
                if (mt[i][j]) {
                    ends[j].push_back(i);
                    break;
                }
        for (size_t i = 0; i < m; i++) {
            if (ends[i].size() < 2)
                continue;
            for (size_t j = 0; j < ends[i].size() - 1; j++)
                for (size_t k = 0; k < m; k++)
                    mt[ends[i][j]][k] = (mt[ends[i][j]][k] ^ mt[ends[i].back()][k]);
        }
    }

    void to_span() {
        change_first();
        for (size_t i = 0; i < m; i++)
            change_last();
        for (size_t i = 0; i < n; i++)
            ones[i] = get_first_one(mt[i]);
    }

    void print() {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++)
                std::cout << mt[i][j];
            std::cout << "\n";
        }
    }

    std::vector<std::pair<size_t, size_t>> get_start_end() {
        std::vector<std::pair<size_t, size_t>> ans;
        for (size_t i = 0; i < n; i++) {
            size_t first_ind = -1;
            size_t second_ind = -1;
            for (size_t j = 0; j < m; j++) {
                if (mt[i][j]) {
                    if (first_ind == -1)
                        first_ind = j;
                    second_ind = j;
                }
            }
            ans.emplace_back(first_ind, second_ind);
        }
        return ans;
    }
};


std::vector<bool>
get_message(matrix &t, std::vector<bool> &word) {
    std::vector<bool> ans(t.n, 0);
    for (size_t i = 0; i < t.n; i++) {
        ans[i] = word[t.ones[i]];
        if (word[t.ones[i]])
            for (size_t j = 0; j < word.size(); j++)
                word[j] = (word[j] ^ t.mt[i][j]);
    }
    return ans;
}

std::vector<std::vector<size_t>> profile(matrix a) {
    auto indexes = a.get_start_end();
    std::vector<std::vector<size_t>> lengths;
    lengths.resize(a.m + 1);
    for (size_t i = 0; i < indexes.size(); i++)
        for (size_t j = indexes[i].first; j < indexes[i].second; j++)
            lengths[j + 1].push_back(i);
    return lengths;
}


void add_bit(size_t &x, size_t i, size_t n) {
    x = (x | (1 << (n - 1 - i)));
}

struct edge {
    size_t nxt;
    int w;

    edge() {
        w = -1;
        nxt = 0;
    };

    edge(size_t nxt1, int w1) : nxt(nxt1), w(w1) {};
};

struct vertex {
    size_t ind;
    edge zer, one;

    vertex() = default;
};

int calc_w(std::vector<std::pair<size_t, bool>> a, std::vector<std::pair<size_t, bool>> b) {
    int res = 0;
    for (size_t i = 0; i < a.size(); i++)
        res = (res ^ (a[i].second & b[i].second));
    return res;
}

void to_vector(unsigned long long &x, std::vector<bool> &ans) {
    for (auto &&an : ans) {
        an = (x & 1);
        x >>= 1;
    }

}

std::vector<std::vector<vertex>> create_graph(matrix t) {
    std::vector<std::vector<vertex>> trellis;
    auto indexes = profile(t);
    trellis.resize(indexes.size());
    for (size_t i = 0; i < trellis.size(); i++) {
        trellis[i].resize(1 << indexes[i].size());
        for (size_t j = 0; j < trellis[i].size(); j++)
            trellis[i][j].ind = j;
    }

    for (size_t i = 0; i < trellis.size() - 1; i++) {
        std::set<size_t> cur_nxt, cur;
        for (size_t j = 0; j < indexes[i + 1].size(); j++)
            cur_nxt.insert(j);
        for (size_t j = 0; j < indexes[i].size(); j++)
            cur.insert(j);
        std::vector<std::pair<size_t, size_t>> same;
        for (size_t fir = 0; fir < indexes[i].size(); ++fir)
            for (size_t sec = 0; sec < indexes[i + 1].size(); sec++)
                if (indexes[i][fir] == indexes[i + 1][sec]) {
                    same.emplace_back(fir, sec);
                    cur_nxt.erase(cur_nxt.find(sec));
                    cur.erase(cur.find(fir));
                    break;
                }
        int diff_nxt = (cur_nxt.size() == 0 ? -1 : *cur_nxt.begin());
        int diff = (cur.size() == 0 ? -1 : *cur.begin());
        std::vector<size_t> need;
        for (size_t j = 0; j < same.size(); j++) {
            need.push_back(indexes[i][same[j].first]);
        }
        if (diff_nxt != -1)
            need.push_back(indexes[i + 1][diff_nxt]);
        if (diff != -1)
            need.push_back(indexes[i][diff]);
        std::vector<std::pair<size_t, bool>> column;
        for (size_t j = 0; j < need.size(); j++)
            column.push_back({need[j], t.mt[need[j]][i]});
        sort(column.begin(), column.end());
        for (size_t j = 0; j < trellis[i].size(); j++) {
            size_t to_gen = 0;
            std::vector<std::pair<size_t, bool>> my_column;
            for (size_t v = 0; v < same.size(); v++)
                if (get_bit(j, same[v].first, indexes[i].size())) {
                    add_bit(to_gen, same[v].second, indexes[i + 1].size());
                    my_column.emplace_back(indexes[i + 1][same[v].second], 1);
                } else
                    my_column.emplace_back(indexes[i + 1][same[v].second], 0);
            if (diff != -1) {
                my_column.emplace_back(indexes[i][diff], get_bit(j, diff, indexes[i].size()));
            }
            if (diff_nxt != -1)
                my_column.emplace_back(indexes[i + 1][diff_nxt], 0);
            sort(my_column.begin(), my_column.end());
            int weight = calc_w(column, my_column);
            trellis[i][j].zer = {to_gen, weight};
            if (diff_nxt != -1) {
                add_bit(to_gen, diff_nxt, indexes[i + 1].size());
                trellis[i][j].one = {to_gen, weight ^ 1};
            }
        }
    }
    return trellis;
}

std::vector<bool> gen_rand_vect(size_t n) {
    std::vector<bool> ans;
    for (size_t i = 0; i < n; i++)
        if (rand() % 2)
            ans.push_back(true);
        else
            ans.push_back(false);
    return ans;
}

int cmp(std::vector<bool> &a, std::vector<bool> &b) {
    int cnt = 0;
    for (size_t i = 0; i < a.size(); i++)
        if (a[i] != b[i])
            return 1;
    return 0;
}

std::vector<bool> code_matrix(matrix &a, const std::vector<bool> &code) {
    std::vector<bool> res;
    for (size_t i = 0; i < a.m; i++) {
        bool ans = 0;
        for (size_t j = 0; j < code.size(); j++)
            ans ^= (a.mt[j][i] & code[j]);
        res.push_back(ans);
    }
    return res;
}

std::vector<double> add_noise(std::vector<bool> &message, std::vector<double> &noise) {
    std::vector<double> res;
    for (size_t i = 0; i < message.size(); i++) {
        if (!message[i])
            res.push_back(-1 + noise[i]);
        else
            res.push_back(1 + noise[i]);
    }
    return res;
}

std::vector<std::vector<bool>> transpose(const std::vector<std::vector<bool>> &x) {
    std::vector<std::vector<bool>> ans;
    ans.resize(x[0].size(), std::vector<bool>(x.size(), 0));
    for (size_t i = 0; i < ans.size(); i++) {
        for (size_t j = 0; j < ans[0].size(); j++)
            ans[i][j] = x[j][i];
    }
    return ans;
}

std::pair<int, int> get_ones(const std::vector<bool> &x) {
    int fir_ind = -1;
    int last_ind = -1;
    for (size_t i = 0; i < x.size(); i++)
        if (x[i]) {
            fir_ind = i;
            break;
        }
    if (fir_ind < 0)
        return {fir_ind, last_ind};
    for (size_t j = x.size() - 1; j >= 0; j--)
        if (x[j]) {
            last_ind = j;
            break;
        }
    return {fir_ind, last_ind};
}

int type_row(const std::vector<bool> &x, const int &l, const int &m, const int &r) {
    auto p = get_ones(x);
    if (p.first == -1)
        return -1;
    if (p.first >= l && p.second <= m)
        return 0;
    if (p.first >= m && p.second <= r)
        return 1;
    if (p.first >= l && p.second <= r)
        return 2;
    return 3;
}

std::vector<bool> copy(const std::vector<bool> &x, const int &l, const int &r) {
    return std::vector<bool>(x.begin() + l, x.begin() + r);
}

std::vector<std::vector<bool>> prepare_matrix(const std::vector<std::vector<bool>> &x, int &cnt) {
    auto t = x;
    int last = 0;
    for (int i = 0; i < t.size(); i++) {
        auto p = get_ones({t[i].begin() + last, t[i].end()});
        if (p.first == -1) {
            continue;
        }
        for (int j = 0; j < t.size(); j++)
            std::swap(t[j][last], t[j][p.first + last]);
        for (size_t j = 0; j < t.size(); j++) {
            if (i == j)
                continue;
            if (t[j][last])
                add(t[j], t[i]);
        }
        last++;
        if (last == t[0].size())
            break;
    }
    std::vector<std::vector<bool>> ans;
    for (size_t i = 0; i < t.size(); i++) {
        if (get_ones(t[i]).first == -1)
            cnt--;
        else
            ans.emplace_back(std::move(x[i]));
    }
    return ans;
}

void gauss(std::vector<std::vector<bool>> a, std::vector<bool> &ans) {
    int n = (int) a.size();
    int m = (int) a[0].size() - 1;

    std::vector<int> where(m, -1);
    for (int col = 0, row = 0; col < m && row < n; ++col) {
        for (int i = row; i < n; ++i)
            if (a[i][col]) {
                swap(a[i], a[row]);
                break;
            }
        if (!a[row][col])
            continue;
        where[col] = row;

        for (int i = 0; i < n; ++i)
            if (i != row && a[i][col])
                add(a[i], a[row]);
        ++row;
    }
    for (size_t i = 0; i < m; i++)
        ans[i] = a[where[i]].back();
}

size_t get_ind(const std::vector<std::vector<bool>> &system_sols, size_t x, size_t sz) {
    std::vector<bool> num;
    int n = system_sols.size();
    for (size_t i = 0; i < n; i++) {
        num.emplace_back((x >> (n - i - 1)) & 1);
    }
    std::vector<bool> answer(system_sols[0].size(), false);
    for (size_t i = 0; i < n; i++)
        if (num[i])
            add(answer, system_sols[i]);
    size_t ans = 0;
    for (size_t i = answer.size() - sz; i < answer.size(); i++) {
        ans *= 2;
        ans += answer[i];
    }
    return ans;
}


#endif //DIPLOMA_COMMON_H
