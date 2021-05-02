#include "common.h"
#include "reed_muller.h"
#include <unordered_map>
#include "polar_encoder.h"

struct branch {
    long long ind_l, ind_r, val;

    bool operator==(branch const &a) const {
        return a.ind_l == ind_l && a.ind_r == ind_r && a.val == val;
    }

    branch(long long _ind_l, long long _ind_r, long long _val) {
        ind_l = _ind_l;
        ind_r = _ind_r;
        val = _val;
    }
};

struct pMatrix {
    bool is_leaf;
    pMatrix *fir, *sec;
    int l, r, third_sz, fourth_sz;
    std::vector<std::vector<bool>> mt;
    std::vector<branch> rules;
    unsigned long long difficult;
    std::vector<std::pair<double, std::vector<bool>>> CBT;
    std::vector<std::vector<std::pair<long long, long long>>> rules_l, rules_r, rules_cos;
    std::vector<bool> used;
    std::vector<size_t> best;
    std::map<std::pair<size_t, size_t>, double> mp;
    std::vector<branch> interested_pairs;

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

        mt = prepare_matrix(x, fourth_sz);
    }

    void makeCBT(const matrix &generator) {
        combCBT(generator, r);
    }

    void prepareLengthOneOrMore() {
        if (r - l == 1) {
            rules.resize(2, {-1, -1, 0});
            CBT.resize(2);
            if (mt.size() == 1)
                rules[1] = {-1, -1, 1};
            return;
        }
        if (mt.size() == 0) {
            CBT.resize(1);
            used.resize(1);
            rules.resize(1ll << (r - l), {-1, -1, -1});
            for (size_t i = 0; i < rules.size(); i++) {
                rules[i] = {-1, -1, 0};
            }
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
        if (!left.empty())
            ans_l.resize(left[0].size() - 1);
        if (right.size() != 0)
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
        CBT.resize(cosets_count, {-INF, std::vector<bool>()});
        used.resize(cosets_count, false);
        rules_l.resize(fir->CBT.size());
        rules_r.resize(sec->CBT.size());
        rules_cos.resize(cosets_count);
//        std::cout << "l=" << l << ' ' << "r=" << r << "\n";
//        std::cout << "masks_count: " << masks_count << "cosets_count: " << cosets_count << "\n";
        for (size_t i = 0; i < masks_count; i++) {

            long long fir_ind = get_ind(left_system_solutions, i, fir->fourth_sz, true);
            long long sec_ind = get_ind(right_system_solutions, i, sec->fourth_sz, true);
            long long rule = i % cosets_count;
            rules.emplace_back(fir_ind, sec_ind, rule);
            rules_l[fir_ind].push_back({sec_ind, rule});
            rules_r[sec_ind].push_back({fir_ind, rule});
            rules_cos[rule].push_back({fir_ind, sec_ind});
        }
    };

    ~pMatrix() {
        if (fir != nullptr)
            delete fir;
        if (sec != nullptr)
            delete sec;

    }


};

void deletePtr(pMatrix *x) {
    if (x->is_leaf)
        delete x;
    else {
        deletePtr(x->fir);
        deletePtr(x->sec);
        delete x;
    }
}

void prepare(pMatrix *x, const matrix &generator, bool f) {
    x->difficult = std::numeric_limits<unsigned long long>::max();
    x->combCBT(generator, -1);
    if (x->mt.size() == 0) {
        x->is_leaf = true;
        return;
    }
    int len = x->r - x->l;
    if (len < 50) {
        long long adds = (len - 1) + (1ll << (len - 1)) - 1;
        long long comps = (1ll << (x->mt.size() - 1)) - (1ll << x->fourth_sz);
        x->difficult = adds + comps;
    }

    int mid = (x->l + x->r) / 2;
    if (len == 2) {
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
        if (len <= 8) {
            x->is_leaf = true;
            return;
        }
    } else {
        if (total < x->difficult) {
            x->difficult = total;
        } else {
            x->is_leaf = true;
        }
    }
}

void run(pMatrix *x, const matrix &generator) {
//    std::cout << x->l << ' ' << x->r << "\n";
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
//    std::cout << x->l << ' ' << x->r << "\n";
}


void Gray(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds, bool f) {
    double total_sum = -data[x->l];
    for (size_t i = x->l + 1; i < x->r; i++) {
        total_sum -= data[i];
        adds++;
    }
    unsigned long long inverse = x->rules.size() - 1;
    unsigned long long cur_code = 0;
    std::vector<bool> to_push(x->r - x->l, false);
    for (size_t i = 0; i < x->rules.size() / 2; i++) {
        unsigned long long cur = cur_code;
        double cur_sum = total_sum;
        if (cur_sum < 0) {
            cur ^= inverse;
            cur_sum = -total_sum;
        }
        if (x->rules[cur].val != -1) {
            if (x->CBT[x->rules[cur].val].first != -INF)
                comps++;
            if (x->rules[cur].val != -1 && x->CBT[x->rules[cur].val].first < total_sum) {
                to_vector(cur, to_push);
                x->CBT[x->rules[cur].val] = {cur_sum, to_push};
                std::fill(to_push.begin(), to_push.end(), false);
            }
            if (x->rules[cur].val != x->rules[cur ^ inverse].val) {
                if (x->CBT[x->rules[cur ^ inverse].val].first != -INF)
                    comps++;
                if (x->CBT[x->rules[cur ^ inverse].val].first < -cur_sum) {
                    to_vector((cur ^ inverse), to_push);
                    x->CBT[x->rules[cur ^ inverse].val] = {-cur_sum, to_push};
                    std::fill(to_push.begin(), to_push.end(), false);
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

std::vector<bool> decode(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds) {
    x->CBT.assign(x->CBT.size(), {-INF, std::vector<bool>()});
    if (x->is_leaf) {
        x->CBT.assign(x->CBT.size(), {-INF, std::vector<bool>()});
        std::vector<bool> to_push(x->r - x->l, false);
        if (x->r - x->l == 1) {
            int curs = (data[x->l] > 0) ? 1 : 0;
            double value = (data[x->l] > 0) ? data[x->l] : -data[x->l];
            to_vector(curs, to_push);
            x->CBT[x->rules[curs].val] = {value, to_push};
            std::fill(to_push.begin(), to_push.end(), false);
            if (x->mt.size() == 1) {
                to_vector(curs ^ 1, to_push);
                x->CBT[x->rules[1 ^ curs].val] = {-value, to_push};
            }
            return {};
        }
        Gray(x, data, comps, adds, false);
    } else {
        decode(x->fir, data, comps, adds);
        decode(x->sec, data, comps, adds);
        for (size_t i = 0; i < x->rules.size(); i++) {
            auto &ind = x->rules[i];
            double sum = x->fir->CBT[ind.ind_l].first + x->sec->CBT[ind.ind_r].first;
            adds++;
            if (i >= x->CBT.size())
                comps++;
            if (x->CBT[ind.val].first < sum) {
                x->CBT[ind.val] = {sum,
                                   concat(x->fir->CBT[ind.ind_l].second,
                                          x->sec->CBT[ind.ind_r].second)};
            }
        }
    }
    return x->CBT[0].second;
}

size_t log2(size_t x) {
    if (x <= 2)
        return 1;
    size_t cnt = 1;
    while (x >= 2) {
        x /= 2;
        cnt++;
    }
    return cnt;
}

void insert_to_vector(pMatrix *x, branch t, long long &comps, long long &adds) {
    for (size_t i = 0; i < x->interested_pairs.size(); i++)
        if (x->interested_pairs[i] == t)
            return;
    if (x->mp[{t.ind_l, t.ind_r}] == 0) {
        adds++;
        x->mp[{t.ind_l, t.ind_r}] = x->fir->CBT[t.ind_l].first + x->sec->CBT[t.ind_r].first;
    }
    double value = x->mp[{t.ind_l, t.ind_r}];
    if (x->interested_pairs.size() == 0) {
        x->interested_pairs.push_back(t);
        return;
    }
    comps += log2(x->interested_pairs.size());
    size_t ind = x->interested_pairs.size();
    for (size_t i = 0; i < x->interested_pairs.size(); i++) {
        if (value > x->mp[{x->interested_pairs[i].ind_l, x->interested_pairs[i].ind_r}]) {
            ind = i;
            break;
        }
    }
    x->interested_pairs.insert(x->interested_pairs.begin() + ind, t);
}

size_t main_decode2(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds, size_t b) {
    if (x->best.size() > b)
        return x->best[b];
    if (x->is_leaf) {
        if (x->best.size() == 0) {
            std::vector<bool> to_push(x->r - x->l, false);
            x->CBT.assign(x->CBT.size(), {-INF, {}});
            if (x->r - x->l == 1) {
                int curs = (data[x->l] > 0) ? 1 : 0;
                double value = (data[x->l] > 0) ? data[x->l] : -data[x->l];
                x->CBT.assign(x->CBT.size(), {-INF, {}});
                to_vector(curs, to_push);
                x->CBT[x->rules[curs].val] = {value, to_push};
                to_push.assign(to_push.size(), false);
                if (x->mt[0][0]) {
                    to_vector(1 ^ curs, to_push);
                    x->CBT[x->rules[1 ^ curs].val] = {-value, to_push};
                }
                return 0;
            }
            Gray(x, data, comps, adds, true);
        }
        return x->best[b];
    } else {

        size_t max_ind_r = x->sec->CBT.size();
        size_t ind_l = 0;
        std::vector<long long> unused_right(x->sec->CBT.size());
        while (max_ind_r != 0) {
            size_t cnt = 0;
            std::fill(unused_right.begin(), unused_right.end(), -1);
            size_t cos_l = main_decode2(x->fir, data, comps, adds, ind_l);
            auto rules_l = x->rules_l[cos_l];
            for (auto &i : rules_l) {
                if (!x->used[i.second]) {
                    unused_right[i.first] = i.second;
                    cnt++;
                }
            }
            if (cnt == 0) {
                ind_l++;
                if (ind_l == x->fir->CBT.size())
                    break;
                continue;
            }
            bool f = false;
            for (size_t j = 0; j < max_ind_r; j++) {
                auto t = main_decode2(x->sec, data, comps, adds, j);
                if (unused_right[t] != -1) {
                    insert_to_vector(x, branch(cos_l, t, unused_right[t]), comps, adds);
                    max_ind_r = j;
                    f = true;
//                    std::cout << x->l << ' ' << x->r << ' ' << ind_l << ' ' << max_ind_r << "\n";
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
        std::vector<long long> unused_left(x->fir->CBT.size(), -1);
        while (max_ind_l != 0) {
            int cnt = 0;
            size_t cos_r = main_decode2(x->sec, data, comps, adds, ind_r);
            std::fill(unused_left.begin(), unused_left.end(), -1);
            auto rules_r = x->rules_r[cos_r];
            for (auto &i : rules_r) {
                if (!x->used[i.second]) {
                    unused_left[i.first] = i.second;
                    cnt++;
                }
            }
            if (cnt == 0) {
                ind_r++;
                if (ind_r == x->sec->CBT.size())
                    break;
                continue;
            }
            bool f = false;
            for (size_t j = 0; j < max_ind_l; j++) {
                auto t = main_decode2(x->fir, data, comps, adds, j);
                if (unused_left[t] != -1) {
                    insert_to_vector(x, branch(t, cos_r, unused_left[t]), comps, adds);
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
        //assert(x->interested_pairs.size() != 0);
        double mx = x->mp[{x->interested_pairs[0].ind_l, x->interested_pairs[0].ind_r}];
        std::pair<size_t, size_t> best_pair = {x->interested_pairs[0].ind_l, x->interested_pairs[0].ind_r};
        size_t ans = x->interested_pairs[0].val;
        x->used[ans] = true;
        x->best.push_back(ans);
        x->CBT[ans] = {mx, concat(x->fir->CBT[best_pair.first].second,
                                  x->sec->CBT[best_pair.second].second)};
        x->interested_pairs.erase(x->interested_pairs.begin());
        std::vector<branch> tmp;
        size_t need = 0;
        for (auto interested_pair : x->interested_pairs)
            if (interested_pair.val != ans) {
                tmp.push_back(interested_pair);
            } else
                need = 1;
        if (need)
            x->interested_pairs = tmp;
        return ans;
    }
}

void clear_structures(pMatrix *x) {
    x->best.clear();
    std::fill(x->used.begin(), x->used.end(), false);
    std::fill(x->CBT.begin(), x->CBT.end(), std::make_pair(-INF, std::vector<bool>()));
    x->mp.clear();
    x->interested_pairs.clear();
    if (x->is_leaf)
        return;
    clear_structures(x->fir);
    clear_structures(x->sec);
}

std::vector<bool> decode2(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds) {
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
    for (double Eb_N0_dB = 0.0; Eb_N0_dB <= 6.0; Eb_N0_dB += 1.) {
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
            auto recieved = (f) ? decode2(ptr, x, comps, adds) : decode(ptr, x, comps, adds);
            auto decoded = get_message(t, recieved);
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

void check_polar(size_t n, size_t k, bool f) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    long long comps = 0, adds = 0;
    PolarEncoder q = PolarEncoder(n, k);
    for (double Eb_N0_dB = 0.0; Eb_N0_dB <= 5.0; Eb_N0_dB += 1.) {
        double sigma_square = 0.5 * ((double) n / k) * ((double) pow(10.0, -Eb_N0_dB / 10));
        std::normal_distribution<> d{0, sqrt(sigma_square)};
        q.reuse_frozen(sqrt(sigma_square));
        auto v = q.getRealGenMatrix();
        matrix t(v);
        auto *ptr = new pMatrix(0, t.m);
        t.to_span();
        prepare(ptr, t, f);
        std::cout << ptr->difficult << "\n";
        run(ptr, t);

        std::cout << "Created \n";
        int cnt = 0;
        comps = 0;
        adds = 0;
        for (size_t i = 0; i < ITER; i++) {
            std::vector<bool> word = gen_rand_vect(t.n);
            std::vector<bool> coded = mulVectorMatrix(word, v);
//            std::vector<bool> new_coded = mulVectorMatrix(word, v);

            std::vector<double> noise;
            noise.reserve(coded.size());
            for (size_t j = 0; j < coded.size(); j++)
                noise.push_back(d(gen));
            auto x = add_noise(coded, noise);
            auto recieve = (f) ? decode2(ptr, x, comps, adds) : decode(ptr, x, comps, adds);
            cnt += cmp(recieve, coded);
        }
        std::cout.precision(7);
        std::cout << std::fixed << (int) Eb_N0_dB << ' ' << (double) cnt / ITER << " " << (comps + adds) / ITER
                  << "\n";
//        deletePtr(ptr);
    }
//    std::cout << "Count of adds and cmps:" << (comps + adds) / (ITER * 7) << "\n";
//    std::cout << "Count of adds:" << (adds) / (ITER * 7) << "\n";
//    std::cout << "Count of cmps:" << (comps) / (ITER * 7) << "\n";

}

int main() {
    srand(time(NULL));
//    check(1, 3, true);
//    check(2, 5, true);
//    check(2, 6, false);
//    check(3, 6, true);
//    check(4, 6, false);
    check_polar(256, 128, false);
    return 0;
}

