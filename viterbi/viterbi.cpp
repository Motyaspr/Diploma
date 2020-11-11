//
// Created by Matvey Sprikut on 08.11.2020.
//

#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <algorithm>

using namespace std;

const int INF = 10000;

struct matrix {
    size_t n{}, m{};
    vector<vector<bool>> mt;

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

    void to_span() {

    }

    void print() {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < m; j++)
                cout << mt[i][j];
            cout << "\n";
        }
    }

    vector<pair<size_t, size_t>> get_start_end() {
        vector<pair<size_t, size_t>> ans;
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

vector<vector<size_t>> profile(matrix a) {
    auto indexes = a.get_start_end();
    vector<vector<size_t>> lengths;
    lengths.resize(a.m + 1);
    for (size_t i = 0; i < indexes.size(); i++)
        for (size_t j = indexes[i].first; j < indexes[i].second; j++)
            lengths[j + 1].push_back(i);
    return lengths;
}

bool get_bit(size_t x, size_t i, size_t n) {
    return ((x >> (n - 1 - i)) & 1);
}

void add_bit(size_t &x, size_t i, size_t n) {
    //cout << x << ' ' << i << ' ' << n << "\n";
    x = (x | (1 << (n - 1 - i)));
    //cout << x << ' ' << i << ' ' << n << "\n";
}

struct edge {
    int w;
    size_t nxt;

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

int calc_w(vector<pair<size_t, bool>> a, vector<pair<size_t, bool>> b) {
    int res = 0;
    for (size_t i = 0; i < a.size(); i++)
        res = (res ^ (a[i].second & b[i].second));
    return res;
}

vector<vector<vertex>> create_graph(matrix t) {
    vector<vector<vertex>> trellis;
    auto indexes = profile(t);
    trellis.resize(indexes.size());
    for (size_t i = 0; i < trellis.size(); i++) {
        trellis[i].resize(1 << indexes[i].size());
        for (size_t j = 0; j < trellis[i].size(); j++)
            trellis[i][j].ind = j;
    }

    for (size_t i = 0; i < trellis.size() - 1; i++) {
        set<size_t> cur_nxt, cur;
        for (size_t j = 0; j < indexes[i + 1].size(); j++)
            cur_nxt.insert(j);
        for (size_t j = 0; j < indexes[i].size(); j++)
            cur.insert(j);
        vector<pair<size_t, size_t>> same;
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
        vector<size_t> need;
        for (size_t j = 0; j < same.size(); j++) {
            need.push_back(indexes[i][same[j].first]);
        }
        if (diff_nxt != -1)
            need.push_back(indexes[i + 1][diff_nxt]);
        if (diff != -1)
            need.push_back(indexes[i][diff]);
        vector<pair<size_t, bool>> column;
        for (size_t j = 0; j < need.size(); j++)
            column.push_back({need[j], t.mt[need[j]][i]});
        sort(column.begin(), column.end());
        for (size_t j = 0; j < trellis[i].size(); j++) {
            size_t to_gen = 0;
            vector<pair<size_t, bool>> my_column;
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

vector<size_t> decode(matrix t, vector<bool> code) {
    auto trellis = create_graph(t);
    vector<vector<int>> dp, p;
    dp.resize(trellis.size());
    p.resize(trellis.size());
    for (size_t i = 0; i < trellis.size(); i++) {
        dp[i].resize(trellis[i].size(), INF);
        p[i].resize(trellis[i].size(), -1);
    }
    dp[0][0] = 0;
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
    vector<size_t> errs;
    size_t curj = 0;
    for (size_t i = dp.size() - 1; i >= 1; i--) {
        if (dp[i][curj] != dp[i - 1][p[i][curj]])
            errs.push_back(i - 1);
        curj = p[i][curj];
    }
    reverse(errs.begin(), errs.end());
    return errs;
}


int main() {
    string s;
    s = "1110010100000000"
        "0110011110000000"
        "0011100000000000"
        "0000101100100000"
        "0000011101011000"
        "0000000101110100"
        "0000000011110110"
        "0000000000001101";
    matrix t(8, 16);
    vector<bool> code{0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0};
    for (size_t i = 0; i < 8; i++)
        for (size_t j = 0; j < 16; j++)
            t.mt[i][j] = (s[16 * i + j] == '1');
    t.print();
    auto ans = decode(t, code);
    cout << ans.size() << "\n";
    for (size_t i = 0; i < ans.size(); i++)
        cout << ans[i] << ' ';
}


