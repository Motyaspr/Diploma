//
// Created by Matvey Sprikut on 08.11.2020.
//

#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <algorithm>
#include "reed_muller.h"
#include <iomanip>
#include <map>
#include <random>
#include <cmath>


const int INF = 10000;

size_t get_first_one(std::vector<bool> t) {
    for (size_t i = 0; i < t.size(); i++)
        if (t[i])
            return i;
    return INF;
}

struct matrix {
    size_t n{}, m{};
    std::vector<std::vector<bool>> mt;

    matrix() {
        n = 0, m = 0;
        mt.resize(n);
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

    edge(size_t  nxt1, int w1) : nxt(nxt1), w(w1) {};
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

double get_prob(double x, bool y, double sigma_2) {
    double t = 1 / sqrt(2 * M_PI * sigma_2);
    double v = (y == 1) ? exp(-((x - 1) * (x - 1)) / (2 * sigma_2)) : exp(-((x + 1) * (x + 1)) / (2 * sigma_2));
    return t * v;
}

std::vector<size_t> decode_bool(matrix &t, std::vector<bool> &code) {
    auto trellis = create_graph(t);
    std::vector<std::vector<int>> dp, p;
    dp.resize(trellis.size());
    p.resize(trellis.size());
    for (size_t i = 0; i < trellis.size(); i++) {
        dp[i].resize(trellis[i].size(), 0.0);
        p[i].resize(trellis[i].size(), -1);
    }
    dp[0][0] = 1.0;
    for (size_t i = 0; i < dp.size() - 1; i++)
        for (size_t j = 0; j < dp[i].size(); j++) {

            if (dp[i + 1][trellis[i][j].zer.nxt] > dp[i][j] + (trellis[i][j].zer.w ^ code[i])) {
                dp[i + 1][trellis[i][j].zer.nxt] = dp[i][j] + (trellis[i][j].zer.w ^ code[i]);
                p[i + 1][trellis[i][j].zer.nxt] = j;
            }
            if (trellis[i][j].one.w != -1) {
                if (dp[i + 1][trellis[i][j].one.nxt] > dp[i][j] + (trellis[i][j].one.w ^ code[i])) {
                    dp[i + 1][trellis[i][j].one.nxt] = dp[i][j] + (trellis[i][j].one.w ^ code[i]);
                    p[i + 1][trellis[i][j].one.nxt] = j;
                }
            }
        }
    std::vector<size_t> errs;
    size_t curj = 0;
    for (size_t i = dp.size() - 1; i >= 1; i--) {
        if (dp[i][curj] != dp[i - 1][p[i][curj]])
            errs.push_back(i - 1);
        curj = p[i][curj];
    }
    reverse(errs.begin(), errs.end());
    return errs;
}

std::vector<bool>
decode(std::vector<double> &code, double sigma_2, std::vector<std::vector<vertex>> trellis) {
    std::vector<std::vector<int>> p;
    std::vector<std::vector<double>> dp;
    dp.resize(trellis.size());
    p.resize(trellis.size());
    for (size_t i = 0; i < trellis.size(); i++) {
        dp[i].resize(trellis[i].size(), INF);
        p[i].resize(trellis[i].size(), -1);
    }
    dp[0][0] = 0.0;
    for (size_t i = 0; i < dp.size() - 1; i++) {
        std::pair<double, double> ws = {get_prob(code[i], false, sigma_2), get_prob(code[i], true, sigma_2)};
        for (size_t j = 0; j < dp[i].size(); j++) {
            double w = (trellis[i][j].zer.w == 0) ? ws.first : ws.second;
            //std::cout << code[i] << ' ' << trellis[i][j].zer.w << ' ' << w << "\n";
            //std::cout << 0 << ' ' << trellis[i][j].zer.w << "\n";
            if (dp[i + 1][trellis[i][j].zer.nxt] > dp[i][j] + w) {
                dp[i + 1][trellis[i][j].zer.nxt] = dp[i][j] + w;
                p[i + 1][trellis[i][j].zer.nxt] = j;
            }
            if (trellis[i][j].one.w != -1) {
                w = (trellis[i][j].one.w == 0) ? ws.first : ws.second;;
                if (dp[i + 1][trellis[i][j].one.nxt] > dp[i][j] + w) {
                    dp[i + 1][trellis[i][j].one.nxt] = dp[i][j] + w;
                    p[i + 1][trellis[i][j].one.nxt] = j;
                }
            }
        }
    }
    std::vector<bool> ans;
    size_t curj = 0;
    for (size_t i = dp.size() - 1; i >= 1; i--) {
        size_t new_curj = p[i][curj];
        if (trellis[i - 1][new_curj].zer.nxt == curj)
            ans.push_back(1 - trellis[i - 1][new_curj].zer.w);
        else
            ans.push_back(1 - trellis[i - 1][new_curj].one.w);
        curj = new_curj;
    }
    reverse(ans.begin(), ans.end());
    return ans;
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

bool cmp(std::vector<bool> &a, std::vector<bool> &b) {
    for (size_t i = 0; i < a.size(); i++)
        if (a[i] != b[i])
            return false;
    return true;
}


int main() {
    srand(time(NULL));
    std::random_device rd{};
    std::mt19937 gen{rd()};

    ReedMuller reedMuller(2, 5);
    matrix t(reedMuller.generated);
    t.print();
    t.to_span();
    t.print();

    std::vector<std::vector<vertex>> trellis = create_graph(t);

    for (double Eb_N0_dB = 0; Eb_N0_dB <= 6.0; Eb_N0_dB += 1.0) {
        int cnt = 0;
        for (size_t i = 0; i < 10000000; i++) {
            double sigma_square = 1.0 / (pow(10.0, Eb_N0_dB / 10));
            std::normal_distribution<> d{0, sigma_square};
            std::vector<bool> coded = code_matrix(t, gen_rand_vect(t.n));
            std::vector<double> noise;
            noise.reserve(coded.size());
            for (size_t j = 0; j < coded.size(); j++)
                noise.push_back(d(gen));
            auto x = add_noise(coded, noise);
            auto decoded = decode(x, sigma_square, trellis);
            if (!cmp(coded, decoded))
                cnt++;
        }
        std::cout << "\n";
        std::cout << Eb_N0_dB << ' ' << (double) cnt / 1000000 << "\n";
    }

}


