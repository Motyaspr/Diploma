#include "common.h"
#include "reed_muller.h"

struct branch {
    long long ind_l, ind_r, val;
};

struct pMatrix {
    bool is_leaf;
    pMatrix *fir, *sec;
    int l, r, third_sz, fourth_sz;
    std::vector<std::vector<bool>> mt;
    std::vector<branch> rules;
    unsigned long long difficult;
    std::vector<std::pair<double, unsigned long long>> CBT;
    std::map<long long, std::vector<std::pair<long long, long long>>> rules_l, rules_r, rules_cos;
    std::vector<bool> used;
    std::vector<size_t> best;

    pMatrix(int _l, int _r) {
        l = _l;
        r = _r;
        is_leaf = false;
        fir = nullptr;
        sec = nullptr;
        third_sz = 0;
        fourth_sz = 0;
    }

    void combCBT(const matrix &generator, const int &mid) {
        third_sz = 0;
        fourth_sz = 0;
        mt.clear();
        std::vector<std::pair<int, std::vector<bool>>> t;
        for (const auto &i : generator.mt) {
            int type = type_row(i, l, mid, r);
            assert(type != -1);
            if (type == 2)
                third_sz++;
            if (type == 3)
                fourth_sz++;
            t.emplace_back(type, copy(i, l, r));
        }
        std::sort(t.begin(), t.end());
        std::vector<std::vector<bool>> x;
        x.reserve(t.size());
        for (auto &i : t)
            x.emplace_back(std::move(i.second));
        mt = std::move(prepare_matrix(x, fourth_sz));
    }

    void makeCBT(const matrix &generator) {
        combCBT(generator, r);
    }

    void prepareLengthOneOrMore() {
        if (r - l == 1) {
            rules.resize(2, {-1, -1, 0});
            CBT.resize(2);
            if (mt[0][0] == 1)
                rules[1] = {-1, -1, 1};
            return;
        }
        long long cosets_count = (1ll << (fourth_sz));
        long long masks_count = (1ll << (mt.size()));
        CBT.resize(cosets_count);
        used.resize(CBT.size(), false);
        rules.resize(1ll << (r - l), {-1, -1, -1});
        for (size_t i = 0; i < masks_count; i++) {
            auto x = get_ind(mt, i, r, false);
            rules[x] = {-1, -1, (long long) i % cosets_count};
        }
    }

    void printMt() {
        std::cout << l << ' ' << r << "\n";
        for (size_t i = 0; i < mt.size(); i++) {
            for (size_t j = 0; j < mt[i].size(); j++)
                std::cout << (int) mt[i][j];
            std::cout << "\n";
        }
        std::cout << "SIZES[2]=" << third_sz << ' ' << "SIZES[3]=" << fourth_sz << "\n";
    }

    void merge() {
        auto left = transpose(fir->mt);
        auto right = transpose(sec->mt);
        for (auto &i : left)
            i.push_back(false);
        for (auto &i : right)
            i.push_back(false);
        std::vector<bool> ans_l, ans_r;
        ans_l.resize(left[0].size() - 1);
        ans_r.resize(right[0].size() - 1);
        std::vector<std::vector<bool>> left_system_solutions, right_system_solutions;
        for (size_t i = mt.size() - third_sz - fourth_sz; i < mt.size(); i++) {
            int ind = 0;
            for (auto &j : left)
                j.back() = mt[i][ind++];
            for (auto &j : right)
                j.back() = mt[i][ind++];
            gauss(left, ans_l);
            left_system_solutions.push_back(ans_l);
            gauss(right, ans_r);
            right_system_solutions.push_back(ans_r);
        }
        long long cosets_count = (1ll << (fourth_sz));
        long long masks_count = (1ll << (third_sz + fourth_sz));
        CBT.resize(cosets_count, {-INF, 0});
        used.resize(cosets_count, false);
        for (size_t i = 0; i < masks_count; i++) {
            long long fir_ind = get_ind(left_system_solutions, i, fir->fourth_sz, true);
            long long sec_ind = get_ind(right_system_solutions, i, sec->fourth_sz, true);
            long long rule = i % cosets_count;
            rules.push_back({fir_ind, sec_ind, rule});
            rules_l[fir_ind].push_back({sec_ind, rule});
            rules_r[sec_ind].push_back({fir_ind, rule});
            rules_cos[rule].push_back({fir_ind, sec_ind});
        }
    }

};

void prepare(pMatrix *x, const matrix &generator, bool f) {
    x->difficult = std::numeric_limits<unsigned long long>::max();
    x->combCBT(generator, -1);
    int len = x->r - x->l;
    if (len < 50) {
        long long adds = (len - 1) + (1ll << (len - 1)) - 1;
        long long comps = (1ll << (x->mt.size() - 1)) - (1ll << x->fourth_sz);
        x->difficult = adds + comps;
    }

    int mid = (x->l + x->r) / 2;
    if (len == 1) {
        x->is_leaf = true;
        return;
    }
    x->sec = new pMatrix(mid, x->r);
    x->fir = new pMatrix(x->l, mid);
    prepare(x->sec, generator, f);
    prepare(x->fir, generator, f);
    x->combCBT(generator, mid);
    unsigned long long total =
            x->fir->difficult + x->sec->difficult + ((2ull << (x->third_sz + x->fourth_sz)) - (1ull << (x->fourth_sz)));
    if (f) {
        if (len <= 4) {
            x->is_leaf = true;
            return;
        }
    } else {
        if (total < x->difficult) {
            x->difficult = total;
        } else {
            x->is_leaf = true;
            delete x->fir;
            delete x->sec;
        }
    }
}

void run(pMatrix *x, const matrix &generator) {
    if (x->is_leaf) {
        x->makeCBT(generator);
        x->prepareLengthOneOrMore();
//        x->printMt();
        return;
    }
    int mid = (x->l + x->r) / 2;
    run(x->sec, generator);
    run(x->fir, generator);
    x->combCBT(generator, mid);
//    x->printMt();
    x->merge();
}


void Gray(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds, bool f) {
    double total_sum = -data[x->l];
    for (size_t i = x->l + 1; i < x->r; i++) {
        total_sum -= data[i];
        adds++;
    }
    long long inverse = x->rules.size() - 1;
    long long cur_code = 0;
    for (size_t i = 0; i < x->rules.size() / 2; i++) {
        long long cur = cur_code;
        double cur_sum = total_sum;
        if (cur_sum < 0) {
            cur ^= inverse;
            cur_sum = -total_sum;
        }
        if (x->rules[cur].val != -1) {
            if (x->CBT[x->rules[cur].val].first != -INF)
                comps++;
            if (x->rules[cur].val != -1 && x->CBT[x->rules[cur].val].first < total_sum) {
                x->CBT[x->rules[cur].val] = {cur_sum, cur};
            }
            if (x->rules[cur].val != x->rules[cur ^ inverse].val) {
                if (x->CBT[x->rules[cur ^ inverse].val].first != -INF)
                    comps++;
                if (x->CBT[x->rules[cur ^ inverse].val].first < -cur_sum) {
                    x->CBT[x->rules[cur ^ inverse].val] = {-cur_sum, cur ^ inverse};
                }
            }
        }
        long long pred = cur_code;
        cur_code = ((i + 1) ^ ((i + 1) >> 1));
        auto diff = (cur_code ^ pred);
        for (size_t j = x->l; j < x->r; j++)
            if ((diff >> (j - x->l)) & 1) {
                if ((cur_code >> (j - x->l)) & 1)
                    total_sum += 2 * data[j];
                else
                    total_sum -= 2 * data[j];
                break;
            }
        adds++;
    }
    adds--;
    if (f) {
        std::set<std::pair<double, size_t>> s;
        for (size_t i = 0; i < x->CBT.size(); i++)
            s.insert({-x->CBT[i].first, i});
        for (auto it : s) {
            x->best.push_back(it.second);
            x->used[it.second] = true;
        }
    }

}

unsigned long long decode(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds) {
    if (x->is_leaf) {
        x->CBT.assign(x->CBT.size(), {-INF, 0});
        if (x->r - x->l == 1) {
            int curs = (data[x->l] > 0) ? 1 : 0;
            double value = (data[x->l] > 0) ? data[x->l] : -data[x->l];
            x->CBT.assign(x->CBT.size(), {-100, 0});
            x->CBT[x->rules[curs].val] = {value, curs};
            if (x->mt[0][0]) {
                x->CBT[x->rules[1 ^ curs].val] = {-value, 1 ^ curs};
            }
            return 0;
        }
        Gray(x, data, comps, adds, false);
    } else {
        int mid = (x->l + x->r) / 2;
        decode(x->fir, data, comps, adds);
        decode(x->sec, data, comps, adds);
        for (size_t i = 0; i < x->rules.size(); i++) {
            auto &ind = x->rules[i];
            double sum = x->fir->CBT[ind.ind_l].first + x->sec->CBT[ind.ind_r].first;
            adds++;
            if (i < x->CBT.size()) {
                x->CBT[ind.val] = {sum,
                                   x->fir->CBT[ind.ind_l].second +
                                   ((x->sec->CBT[ind.ind_r].second) << (mid - x->l))};

            } else {
                comps++;
                if (x->CBT[ind.val].first < sum) {
                    x->CBT[ind.val] = {sum,
                                       x->fir->CBT[ind.ind_l].second +
                                       ((x->sec->CBT[ind.ind_r].second) << (mid - x->l))};
                }
            }
        }
    }
    return x->CBT[0].second;
}

size_t main_decode2(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds, size_t b) {
    if (x->best.size() > b)
        return x->best[b];
    if (x->is_leaf) {
        if (x->best.size() == 0) {
            x->CBT.assign(x->CBT.size(), {-INF, 0});
            if (x->r - x->l == 1) {
                int curs = (data[x->l] > 0) ? 1 : 0;
                double value = (data[x->l] > 0) ? data[x->l] : -data[x->l];
                x->CBT.assign(x->CBT.size(), {-100, 0});
                x->CBT[x->rules[curs].val] = {value, curs};
                if (x->mt[0][0]) {
                    x->CBT[x->rules[1 ^ curs].val] = {-value, 1 ^ curs};
                }
                return 0;
            }
            Gray(x, data, comps, adds, true);
        }
        return x->best[b];
    } else {
        int mid = (x->l + x->r) / 2;
        size_t max_ind_r = x->sec->CBT.size();
        size_t ind_l = 0;
        std::set<std::pair<size_t, size_t>> interested_pair;
        std::vector<bool> unused_right(x->sec->CBT.size(), false);
        while (max_ind_r != 0) {
            size_t cnt = 0;
            std::fill(unused_right.begin(), unused_right.end(), false);
            size_t cos_l = main_decode2(x->fir, data, comps, adds, ind_l);
            auto rules_l = x->rules_l[cos_l];
            for (auto &i : rules_l) {
                if (!x->used[i.second]) {
                    unused_right[i.first] = true;
                    cnt++;
                }
            }
            if (cnt == 0) {
                ind_l++;
                break;
                if (ind_l == x->fir->CBT.size())
                    break;
                continue;
            }
            bool f = false;
            for (size_t j = 0; j < max_ind_r; j++) {
                auto t = main_decode2(x->sec, data, comps, adds, j);
                if (unused_right[t]) {
                    interested_pair.insert({cos_l, t});
                    max_ind_r = j;
                    f = true;
                    break;
                }
            }
            ind_l++;
            if (ind_l == x->fir->CBT.size())
                break;
            if (!f)
                break;
        }
        size_t max_ind_l = x->fir->CBT.size();
        size_t ind_r = 0;
        std::vector<bool> unused_left(x->fir->CBT.size(), false);
        while (max_ind_l != 0) {
            int cnt = 0;
            size_t cos_r = main_decode2(x->sec, data, comps, adds, ind_r);
            std::fill(unused_left.begin(), unused_left.end(), false);
            auto rules_r = x->rules_r[cos_r];
            for (auto &i : rules_r) {
                if (!x->used[i.second]) {
                    unused_left[i.first] = true;
                    cnt++;
                }
            }
            if (cnt == 0) {
                ind_r++;
                break;
                if (ind_r == x->sec->CBT.size())
                    break;
                continue;
            }
            bool f = false;
            for (size_t j = 0; j < max_ind_l; j++) {
                auto t = main_decode2(x->fir, data, comps, adds, j);
                if (unused_left[t]) {
                    interested_pair.insert({t, cos_r});
                    max_ind_l = j;
                    f = true;
                    break;
                }
            }
            ind_r++;
            if (ind_r == x->sec->CBT.size())
                break;
            if (!f)
                break;
        }
        double mx = 0;
        std::pair<size_t, size_t> best_pair;
        assert(interested_pair.size() != 0);
        for (auto it : interested_pair) {
            adds++;
            comps++;
            if (x->fir->CBT[it.first].first + x->sec->CBT[it.second].first > mx) {
                mx = x->fir->CBT[it.first].first + x->sec->CBT[it.second].first;
                best_pair = it;
            }
        }
        comps--;
        size_t ans = 0;
        for (size_t i = 0; i < x->rules_l[best_pair.first].size(); i++)
            if (x->rules_l[best_pair.first][i].first == best_pair.second) {
                ans = x->rules_l[best_pair.first][i].second;
                x->used[ans] = true;
                x->best.push_back(ans);
                x->CBT[ans] = {mx, x->fir->CBT[best_pair.first].second +
                                   (x->sec->CBT[best_pair.second].second << (mid - x->l))};
                break;
            }
        return ans;
    }
}

void clear_structures(pMatrix *x) {
    x->best.clear();
    std::fill(x->used.begin(), x->used.end(), false);
    std::fill(x->CBT.begin(), x->CBT.end(), std::make_pair(-INF, 0));
    if (x->is_leaf)
        return;
    clear_structures(x->fir);
    clear_structures(x->sec);
}

unsigned long long decode2(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds) {
    clear_structures(x);
    return x->CBT[main_decode2(x, data, comps, adds, 0)].second;
}


void check(int r, int m, bool f) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    ReedMuller reedMuller(r, m);
    matrix t(reedMuller.generated);
    auto *ptr = new pMatrix(0, t.m);
    t.to_span();
    prepare(ptr, t, f);
    run(ptr, t);
    std::cout << "RM(" << r << ", " << m << ") created\n";
    std::cout << ptr->difficult << "\n";
    long long comps = 0, adds = 0;
    for (double Eb_N0_dB = 2.0; Eb_N0_dB <= 6.0; Eb_N0_dB += 1.) {
        int cnt = 0;
        comps = 0;
        adds = 0;
        double sigma_square = 0.5 * ((double) t.m / t.n) * ((double) pow(10.0, -Eb_N0_dB / 10));
        std::normal_distribution<> d{0, sqrt(sigma_square)};
        for (size_t i = 0; i < ITER; i++) {
            std::vector<bool> word = gen_rand_vect(t.n);
            std::vector<bool> coded = code_matrix(t, word);
            std::vector<double> noise;
            noise.reserve(coded.size());
            for (size_t j = 0; j < coded.size(); j++)
                noise.push_back(d(gen));
            auto x = add_noise(coded, noise);
            std::vector<bool> recieve(coded.size(), false);

            auto ans = (f) ? decode2(ptr, x, comps, adds) : decode(ptr, x, comps, adds);
            to_vector(ans, recieve);
            auto decoded = get_message(t, recieve);
            cnt += cmp(decoded, word);
        }
        std::cout.precision(7);
        std::cout << std::fixed << (int) Eb_N0_dB << ' ' << (double) cnt / ITER << " " << (comps + adds) / ITER
                  << "\n";
    }
//    std::cout << "Count of adds and cmps:" << (comps + adds) / (ITER * 7) << "\n";
//    std::cout << "Count of adds:" << (adds) / (ITER * 7) << "\n";
//    std::cout << "Count of cmps:" << (comps) / (ITER * 7) << "\n";
    delete ptr;
}

int main() {
    srand(time(NULL));
//    check(1, 3);
//    check(2, 5);
//    check(2, 6);
    check(3, 6, true);
//    check(4, 6);
    return 0;
}

