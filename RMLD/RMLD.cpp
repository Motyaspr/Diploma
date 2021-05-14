#include "common.h"
#include "reed_muller.h"
#include <unordered_map>
#include "polar_encoder.h"
#include "thread"

const int threads_count = 32;

const int MIN_LEN_LEAF = 8;

std::mutex myMutex;

struct branch {
    long long ind_l, ind_r, val;

    branch() = default;

    bool operator==(branch const &a) const {
        return a.ind_l == ind_l && a.ind_r == ind_r && a.val == val;
    }

    bool operator<(branch const &a) const {
        return (a.ind_l > ind_l || ((a.ind_l == ind_l) && a.ind_r > ind_r));
    }

    branch(long long _ind_l, long long _ind_r, long long _val) {
        ind_l = _ind_l;
        ind_r = _ind_r;
        val = _val;
    }
};


struct pair_hash {
    template<class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // Mainly for demonstration purposes, i.e. works but is overly simple
        // In the real world, use sth. like boost.hash_combine
        return h1 ^ h2;
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
    std::vector<std::vector<std::pair<long long, long long>>> rules_l, rules_r, rules_l2, rules_r2;
    std::vector<std::vector<bool>> solutions_for_left_part, oneD_solutions_for_left_part;
    std::vector<std::vector<bool>> solutions_for_right_part, oneD_solutions_for_right_part;
    std::vector<size_t> best;
    std::unordered_map<std::pair<size_t, size_t>, double, pair_hash> mp;
    std::set<branch> interested_pairs;
    std::unordered_map<long long, long long> unused_left, unused_right;
    std::vector<std::vector<bool>> left_system_solutions, right_system_solutions;

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
            rules.resize(1ll << (r - l), {-1, -1, -1});
            for (size_t i = 0; i < rules.size(); i++) {
                rules[i] = {-1, -1, 0};
            }
            return;
        }
        long long cosets_count = (1ll << (fourth_sz));
        long long masks_count = (1ll << (mt.size()));
        CBT.resize(cosets_count);
        rules.resize(1ll << (r - l), {-1, -1, -1});
        for (size_t i = 0; i < masks_count; i++) {
            auto x = get_ind(mt, i, r, false, false);
            rules[x] = {-1, -1, (long long) i % cosets_count};
        }
    }

    void printMt() {
        //std::cout << l << ' ' << r << "\n";
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
        std::vector<std::vector<bool>> mt1, mt2;
//        std::cout << "MT: " << mt.size() << ' ' << mt[0].size() << "\n";
//        printMt();
//        std::cout << "fir " << fir->mt.size() << ' ' << fir->mt[0].size() << "\n";
//        fir->printMt();
//        std::cout << "fir cbt: " << fir->CBT.size() << "\n";
//        std::cout << "sec " << sec->mt.size() << ' ' << sec->mt[0].size() << "\n";
//        sec->printMt();
//        std::cout << "a'size =" << fir->mt.size() - fir->fourth_sz << "\n";
        std::vector<std::vector<std::vector<bool>>> all_left_creation;
        for (size_t i = mt.size() - fourth_sz - third_sz; i < mt.size(); i++) {
            std::vector<bool> half1;
            std::vector<bool> half2;
            for (size_t j = 0; j < mt[0].size(); j++)
                if (j < fir->mt[0].size())
                    half1.push_back(mt[i][j]);
                else
                    half2.push_back(mt[i][j]);
            mt1.push_back(half1);
            mt2.push_back(half2);
        }
        auto mt1T = transpose(mt1);
        auto mt2T = transpose(mt2);
        for (size_t i = 0; i < fir->mt.size() - fir->fourth_sz; i++)
            for (size_t j = 0; j < fir->mt[i].size(); j++)
                mt1T[j].push_back(fir->mt[i][j]);
        for (size_t i = 0; i < sec->mt.size() - sec->fourth_sz; i++)
            for (size_t j = 0; j < sec->mt[i].size(); j++)
                mt2T[j].push_back(sec->mt[i][j]);
        auto q = get_free(mt1T);
        for (size_t i = 0; i < q.size(); i++) {
            std::vector<bool> tmp(mt1T[0].size(), false);
            tmp[q[i]] = true;
            mt1T.push_back(tmp);
        }
        q = get_free(mt2T);
        for (size_t i = 0; i < q.size(); i++) {
            std::vector<bool> tmp(mt2T[0].size(), false);
            tmp[q[i]] = true;
            mt2T.push_back(tmp);
        }
        for (size_t i = 0; i < mt1T.size(); i++)
            mt1T[i].push_back(false);
        for (size_t i = 0; i < mt2T.size(); i++)
            mt2T[i].push_back(false);

        for (auto &i : left)
            i.push_back(false);
        for (auto &i : right)
            i.push_back(false);
        std::vector<bool> ans_l, ans_r;
        if (!left.empty())
            ans_l.resize(left[0].size() - 1);
        if (right.size() != 0)
            ans_r.resize(right[0].size() - 1);
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
//        CBT.resize(cosets_count, {-INF, std::vector<bool>()});
        rules_l.resize(1 << fir->fourth_sz);
        rules_r.resize(1 << sec->fourth_sz);
        rules.resize(masks_count);
        update(0, masks_count, cosets_count);
        std::vector<bool> ans_(mt1T[0].size() - 1, false);
        std::vector<bool> ans_2(mt2T[0].size() - 1, false);

        for (size_t j = sec->mt.size() - sec->fourth_sz; j < sec->mt.size(); j++) {
            auto answer = sec->mt[j];
            for (size_t k = 0; k < sec->mt[j].size(); k++)
                mt2T[k].back() = answer[k];
            assert(gauss(mt2T, ans_2));
            solutions_for_right_part.push_back(ans_2);
        }

        for (size_t j = fir->mt.size() - fir->fourth_sz; j < fir->mt.size(); j++) {
            auto answer = fir->mt[j];
            for (size_t k = 0; k < fir->mt[j].size(); k++)
                mt1T[k].back() = answer[k];
            assert(gauss(mt1T, ans_));
            solutions_for_left_part.push_back(ans_);
        }
        for (size_t i = 0; i < mt1T.size(); i++)
            mt1T[i].back() = false;
        for (size_t i = 0; i < mt2T.size(); i++)
            mt2T[i].back() = false;

        for (size_t i = fir->mt[0].size(); i < mt1T.size(); i++) {
            mt1T[i].back() = true;
            gauss(mt1T, ans_);
            oneD_solutions_for_left_part.push_back(ans_);
            mt1T[i].back() = false;
        }
        for (size_t i = sec->mt[0].size(); i < mt2T.size(); i++) {
            mt2T[i].back() = true;
            gauss(mt2T, ans_2);
            oneD_solutions_for_right_part.push_back(ans_2);
            mt2T[i].back() = false;
        }
//        for (size_t i = 0; i < (1ll << fir->fourth_sz); i++) {
//            rules_l2.push_back(to_mul_left(i));
//        }
//        for (size_t i = 0; i < (1ll << sec->fourth_sz); i++) {
//            rules_r2.push_back(to_mul_right(i));
//        }
//        for (size_t i = 0; i < rules_l.size(); i++) {
//            std::sort(rules_l[i].begin(), rules_l[i].end());
//            std::sort(rules_l2[i].begin(), rules_l2[i].end());
//            assert(rules_l[i] == rules_l2[i]);
//        }
//        for (size_t i = 0; i < rules_r.size(); i++) {
//            std::sort(rules_r[i].begin(), rules_r[i].end());
//            std::sort(rules_r2[i].begin(), rules_r2[i].end());
//            assert(rules_r[i] == rules_r2[i]);
//        }
//        auto action = [this](const std::vector<std::vector<bool>> &left_system_solutions,
//                             const std::vector<std::vector<bool>> &right_system_solutions, size_t mini, size_t maxi,
//                             size_t cosets_count) {
//            update(left_system_solutions, right_system_solutions, mini, maxi, cosets_count);
//        };
//
//        std::vector<std::thread> threads;
//        for (size_t i = 0; i < threads_count; i++) {
//            size_t mini = (masks_count / threads_count) * i;
//            size_t maxi = (masks_count / threads_count) * (i + 1);
//            if (i == threads_count - 1)
//                maxi = masks_count;
//            threads.push_back(
//                    std::thread(action, std::ref(left_system_solutions), std::ref(right_system_solutions), mini, maxi,
//                                cosets_count));
//        }
//        for (size_t i = 0; i < threads.size(); i++)
//            threads[i].join();
    }

    std::vector<std::pair<long long, long long>> to_mul_left(long long x) {
        std::vector<std::pair<long long, long long>> all_pairs;
        std::vector<bool> part_l(solutions_for_left_part.size());
        to_vector2(x, part_l);
        auto answer = mulVectorMatrix(part_l, solutions_for_left_part);
        for (size_t i = 0; i < (1ll << oneD_solutions_for_left_part.size()); i++) {
            std::vector<bool> summer(answer.size(), false);
            for (size_t j = 0; j < oneD_solutions_for_left_part.size(); j++)
                if (get_bit(i, j, oneD_solutions_for_left_part.size()))
                    add(summer, oneD_solutions_for_left_part[j]);
            add(summer, answer);
            long long coset = 0;
            for (size_t i = 0; i < summer.size(); i++)
                coset = coset * 2 + summer[i];
            all_pairs.push_back(
                    {get_ind(right_system_solutions, reinterpret_cast<size_t &>(coset), sec->fourth_sz, true, false),
                     coset % (1 << fourth_sz)});
        }
        return all_pairs;
    }

    std::vector<std::pair<long long, long long>> to_mul_right(long long x) {
        std::vector<std::pair<long long, long long>> all_pairs;
        std::vector<bool> part_r(solutions_for_right_part.size());
        to_vector2(x, part_r);
        auto answer = mulVectorMatrix(part_r, solutions_for_right_part);
        for (size_t i = 0; i < (1ll << oneD_solutions_for_right_part.size()); i++) {
            std::vector<bool> summer(answer.size(), false);
            for (size_t j = 0; j < oneD_solutions_for_right_part.size(); j++)
                if (get_bit(i, j, oneD_solutions_for_right_part.size()))
                    add(summer, oneD_solutions_for_right_part[j]);
            add(summer, answer);
            long long coset = 0;
            for (size_t i = 0; i < summer.size(); i++)
                coset = coset * 2 + summer[i];
            all_pairs.push_back(
                    {get_ind(left_system_solutions, reinterpret_cast<size_t &>(coset), fir->fourth_sz, true, false),
                     coset % (1 << fourth_sz)});
        }
        return all_pairs;
    }

    void update(size_t mini, size_t maxi,
                size_t cosets_count) {
        for (size_t i = mini; i < maxi; i++) {
            long long fir_ind = get_ind(left_system_solutions, i, fir->fourth_sz, true, true);
            long long sec_ind = get_ind(right_system_solutions, i, sec->fourth_sz, true, false);
            long long rule = i % cosets_count;
            rules[i] = branch(fir_ind, sec_ind, rule);
            std::lock_guard<std::mutex> myLock(myMutex);
            rules_l[fir_ind].push_back({sec_ind, rule});
            rules_r[sec_ind].push_back({fir_ind, rule});
        }
    }

    ~pMatrix() {
        if (fir != nullptr)
            delete fir;
        if (sec != nullptr)
            delete sec;

    }
};

struct CBT {
    pMatrix *x;
    std::vector<size_t> best;
    std::unordered_map<std::pair<size_t, size_t>, double> mp;
    std::vector<branch> interested_pairs;

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
        if (len <= MIN_LEN_LEAF) {
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
    std::cout << x->l << ' ' << x->r << "\n";
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

//void insert_to_vector(pMatrix *x, branch t, long long &comps, long long &adds) {
//    for (size_t i = 0; i < x->interested_pairs.size(); i++)
//        if (x->interested_pairs[i] == t)
//            return;
//    if (x->mp[{t.ind_l, t.ind_r}] == 0) {
//        adds++;
//        x->mp[{t.ind_l, t.ind_r}] = x->fir->CBT[t.ind_l].first + x->sec->CBT[t.ind_r].first;
//    }
//    double value = x->mp[{t.ind_l, t.ind_r}];
//    if (x->interested_pairs.size() == 0) {
//        x->interested_pairs.push_back(t);
//        return;
//    }
//    comps += log2(x->interested_pairs.size());
//    size_t ind = x->interested_pairs.size();
//    for (size_t i = 0; i < x->interested_pairs.size(); i++) {
//        if (value > x->mp[{x->interested_pairs[i].ind_l, x->interested_pairs[i].ind_r}]) {
//            ind = i;
//            break;
//        }
//    }
//    x->interested_pairs.insert(x->interested_pairs.begin() + ind, t);
//}

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
        while (max_ind_r != 0) {
            size_t cnt = 0;
            x->unused_right.clear();
            size_t cos_l = main_decode2(x->fir, data, comps, adds, ind_l);
            auto rules_l = x->rules_l2[cos_l];
            for (auto &i : rules_l) {
                if (x->CBT[i.second].first == -INF) {
                    x->unused_right[i.first] = i.second;
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
                if (x->unused_right.find(t) != x->unused_right.end()) {
                    x->interested_pairs.insert(branch(cos_l, t, x->unused_right[t]));
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
        while (max_ind_l != 0) {
            int cnt = 0;
            size_t cos_r = main_decode2(x->sec, data, comps, adds, ind_r);
            x->unused_left.clear();
            auto rules_r = x->rules_r[cos_r];
            for (auto &i : rules_r) {
                if (x->CBT[i.second].first == -INF) {
                    x->unused_left[i.first] = i.second;
                    cnt++;
                }
            }
            if (cnt == 0) {
                ind_r++;
                if (x->interested_pairs.size() != 0)
                    break;
                if (ind_r == x->sec->CBT.size())
                    break;
                continue;
            }
            bool f = false;
            for (size_t j = 0; j < max_ind_l; j++) {
                auto t = main_decode2(x->fir, data, comps, adds, j);
                if (x->unused_left.find(t) != x->unused_left.end()) {
                    x->interested_pairs.insert(branch(t, cos_r, x->unused_left[t]));
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
//        if (x->interested_pairs.size() == 0)
//            return -1;
        branch best_pair;
        double mx = -1000000;
        for (auto it : x->interested_pairs) {
            double val = x->mp[{it.ind_l, it.ind_r}];
            if (val == 0.0) {
                x->mp[{it.ind_l, it.ind_r}] = x->fir->CBT[it.ind_l].first + x->sec->CBT[it.ind_r].first;
                val = x->mp[{it.ind_l, it.ind_r}];
                adds++;
            }
            comps++;
            if (val > mx) {
                mx = val;
                best_pair = it;
            }
        }
        x->best.push_back(best_pair.val);
        x->CBT[best_pair.val] = {mx, concat(x->fir->CBT[best_pair.ind_l].second, x->sec->CBT[best_pair.ind_r].second)};
        x->interested_pairs.clear();
        return best_pair.val;
    }
}

void clear_structures(pMatrix *x) {
    x->best.clear();
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

void check_polar(size_t n, size_t k, bool f, int cnt_iter) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    long long comps = 0, adds = 0;
    PolarEncoder q = PolarEncoder(n, k);
    for (double Eb_N0_dB = 0.0; Eb_N0_dB <= 6.0; Eb_N0_dB += 1) {

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
        std::cout << std::fixed << (double) Eb_N0_dB << ' ';
        for (size_t i = 0; i < cnt_iter; i++) {
            std::vector<bool> word = gen_rand_vect(t.n);
            std::vector<bool> coded = mulVectorMatrix(word, v);
//            std::vector<bool> new_coded = mulVectorMatrix(word, v);

            std::vector<double> noise;
            noise.reserve(coded.size());
            for (size_t j = 0; j < coded.size(); j++)
                noise.push_back(d(gen));
            auto x = add_noise(coded, noise);
            auto recieve = (f) ? decode2(ptr, x, comps, adds) : decode(ptr, x, comps, adds);
//            if (f)
//                std::cout << i << " " << (adds + comps) / (i + 1) << "\n";
            cnt += cmp(recieve, coded);
        }
        std::cout << "\n";
        std::cout.precision(7);
        std::cout << (double) cnt / (cnt_iter + 1) << " " << (comps + adds) / (cnt_iter + 1)
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
//   check(3, 6, true);
//    check(3, 6, true);
    check_polar(128, 64, false, 0);
//    check_polar(1024, 512, true, 0);

    return 0;
}
