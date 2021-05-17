#include "common.h"
#include "reed_muller.h"
#include <unordered_map>
#include "polar_encoder.h"
#include <thread>
#include <mutex>
#include <cassert>

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
    std::vector<std::pair<double, std::vector<bool>>> CBT_2;
    std::vector<std::vector<std::pair<long long, long long>>> rules_l, rules_r, rules_l2, rules_r2;
    std::vector<std::vector<bool>> solutions_for_left_part, oneD_solutions_for_left_part;
    std::vector<std::vector<bool>> solutions_for_right_part, oneD_solutions_for_right_part;
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
            CBT_2.resize(2);
            if (mt.size() == 1)
                rules[1] = {-1, -1, 1};
            return;
        }
        if (mt.size() == 0) {
            CBT_2.resize(1);
            rules.resize(1ll << (r - l), {-1, -1, -1});
            for (size_t i = 0; i < rules.size(); i++) {
                rules[i] = {-1, -1, 0};
            }
            return;
        }
        long long cosets_count = (1ll << (fourth_sz));
        long long masks_count = (1ll << (mt.size()));
        CBT_2.resize(cosets_count);
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
//        std::cout << "fir cbt: " << fir->CBT_2.size() << "\n";
//        std::cout << "sec " << sec->mt.size() << ' ' << sec->mt[0].size() << "\n";
//        sec->printMt();
//        std::cout << "a'size =" << fir->mt.size() - fir->fourth_sz << "\n";

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
//        long long cosets_count = (1ll << (fourth_sz));
//        long long masks_count = (1ll << (third_sz + fourth_sz));
//        CBT_2.resize(cosets_count, {-INF, std::vector<bool>()});
//        rules_l.resize(1 << fir->fourth_sz);
//        rules_r.resize(1 << sec->fourth_sz);
//        rules.resize(masks_count);
        std::vector<bool> ans_(third_sz + fourth_sz, false);
        std::vector<bool> ans_2(third_sz + fourth_sz, false);

        for (size_t j = sec->mt.size() - sec->fourth_sz; j < sec->mt.size(); j++) {
            auto answer = sec->mt[j];
            for (size_t k = 0; k < sec->mt[j].size(); k++)
                mt2T[k].back() = answer[k];
            gauss(mt2T, ans_2);
            solutions_for_right_part.push_back(ans_2);
        }

        for (size_t j = fir->mt.size() - fir->fourth_sz; j < fir->mt.size(); j++) {
            auto answer = fir->mt[j];
            for (size_t k = 0; k < fir->mt[j].size(); k++)
                mt1T[k].back() = answer[k];
            gauss(mt1T, ans_);
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

    std::vector<std::pair<long long, long long>> to_mul_left(long long x) const {
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
            for (size_t i = 0; i < third_sz + fourth_sz; i++)
                coset = coset * 2 + summer[i];
            all_pairs.push_back(
                    {get_ind(right_system_solutions, reinterpret_cast<size_t &>(coset), sec->fourth_sz, true, false),
                     coset % (1 << fourth_sz)});
        }
        return all_pairs;
    }

    std::vector<std::pair<long long, long long>> to_mul_right(long long x) const {
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
            for (size_t i = 0; i < third_sz + fourth_sz; i++)
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
    const pMatrix *x;
    std::vector<size_t> best;
    std::set<branch> interested_pairs;
    std::map<long long, long long> unused_left, unused_right;
    std::vector<std::pair<double, std::vector<bool>>> final_CBT;
    std::map<long long, long long> index_mapping;
    std::unordered_map<long long, long long> left_index_mapping;
    std::vector<std::vector<std::pair<long long, long long>>> left_neighbors;

    CBT() = default;

    CBT(const pMatrix *x1) {
        x = x1;
    }

    CBT(const CBT& another) {
        x = another.x;
    }

    CBT(CBT&& another) : x(another.x) {
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
//    std::cout << x->l << ' ' << x->r << "\n";
}


void Gray(CBT &t, const std::vector<double> &data, long long &comps, long long &adds, bool f) {
    auto x = t.x;
    t.final_CBT.resize(x->CBT_2.size(), {-INF, std::vector<bool>(x->r - x->l, false)});
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
            if (t.final_CBT[x->rules[cur].val].first != -INF)
                comps++;
            if (x->rules[cur].val != -1 && t.final_CBT[x->rules[cur].val].first < total_sum) {
                to_vector(cur, to_push);
                t.final_CBT[x->rules[cur].val] = {cur_sum, to_push};
                std::fill(to_push.begin(), to_push.end(), false);
            }
            if (x->rules[cur].val != x->rules[cur ^ inverse].val) {
                if (t.final_CBT[x->rules[cur ^ inverse].val].first != -INF)
                    comps++;
                if (t.final_CBT[x->rules[cur ^ inverse].val].first < -cur_sum) {
                    to_vector((cur ^ inverse), to_push);
                    t.final_CBT[x->rules[cur ^ inverse].val] = {-cur_sum, to_push};
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
        for (size_t i = 0; i < t.final_CBT.size(); i++)
            s.insert({-t.final_CBT[i].first, i});
        for (auto it : s) {
            t.best.push_back(it.second);
        }
    }
    for (size_t i = 0; i < t.final_CBT.size(); i++)
        t.index_mapping[i] = i;

}

std::vector<bool>
decode(std::vector<CBT> &r, size_t index, const std::vector<double> &data, long long &comps, long long &adds) {
    CBT &t = r[index];
    t.final_CBT.assign(t.final_CBT.size(), {-INF, std::vector<bool>()});
    auto x = t.x;
    if (x->is_leaf) {
        t.final_CBT.assign(t.final_CBT.size(), {-INF, std::vector<bool>()});
        std::vector<bool> to_push(x->r - x->l, false);
        if (x->r - x->l == 1) {
            int curs = (data[x->l] > 0) ? 1 : 0;
            double value = (data[x->l] > 0) ? data[x->l] : -data[x->l];
            to_vector(curs, to_push);
            t.final_CBT[x->rules[curs].val] = {value, to_push};
            std::fill(to_push.begin(), to_push.end(), false);
            if (x->mt.size() == 1) {
                to_vector(curs ^ 1, to_push);
                t.final_CBT[x->rules[1 ^ curs].val] = {-value, to_push};
            }
            return {};
        }
        Gray(t, data, comps, adds, false);
    } else {
        decode(r, index * 2, data, comps, adds);
        decode(r, index * 2 + 1, data, comps, adds);
        for (size_t i = 0; i < (1 << x->fir->fourth_sz); i++) {
            if (t.left_index_mapping.find(i) == t.left_index_mapping.end()) {
                t.left_index_mapping[i] = t.left_neighbors.size();
                t.left_neighbors.push_back(t.x->to_mul_left(i));
            }
            auto neighbors = t.left_neighbors[t.left_index_mapping[i]];
            for (auto it : neighbors) {
                if (t.index_mapping.find(it.second) == t.index_mapping.end()) {
                    t.index_mapping[it.second] = t.final_CBT.size();
                    t.final_CBT.push_back({-INF, std::vector<bool>()});
                }
                size_t CBT_ind = t.index_mapping[it.second];
                size_t real_first_ind = r[index * 2].index_mapping[i];
                size_t real_second_ind = r[index * 2 + 1].index_mapping[it.first];
//                std::cout << i << ' ' << it.first << ' ' << it.second << "\n";
//                std::cout << real_first_ind << ' ' << real_second_ind << ' ' << CBT_ind << "\n";
                adds++;
                if (t.final_CBT[CBT_ind].first != -INF)
                    comps++;
                double sum = r[index * 2].final_CBT[real_first_ind].first +
                             r[index * 2 + 1].final_CBT[real_second_ind].first;
                if (t.final_CBT[CBT_ind].first < sum) {
                    t.final_CBT[CBT_ind].first = sum;
                    t.final_CBT[CBT_ind].second = concat(r[index * 2].final_CBT[real_first_ind].second,
                                                         r[index * 2 + 1].final_CBT[real_second_ind].second);
                }
            }
        }
//        for (size_t i = 0; i < x->rules.size(); i++) {
//            auto &ind = x->rules[i];
//            double sum = x->fir->CBT_2[ind.ind_l].first + x->sec->CBT_2[ind.ind_r].first;
//            adds++;
//            if (i >= x->CBT_2.size())
//                comps++;
//            if (x->CBT_2[ind.val].first < sum) {
//                x->CBT_2[ind.val] = {sum,
//                                     concat(x->fir->CBT_2[ind.ind_l].second,
//                                            x->sec->CBT_2[ind.ind_r].second)};
//            }
//        }
    }
    return t.final_CBT[0].second;
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
//        x->mp[{t.ind_l, t.ind_r}] = x->fir->CBT_2[t.ind_l].first + x->sec->CBT_2[t.ind_r].first;
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

long long
main_decode2(std::vector<CBT> &r, size_t index, const std::vector<double> &data, long long &comps, long long &adds,
             size_t b, size_t depth, size_t cnt1, size_t cnt2) {
    CBT &t = r[index];
    if (t.best.size() > b)
        return t.best[b];
    auto x = t.x;
    if (x->is_leaf) {
        if (t.best.size() == 0) {
            std::vector<bool> to_push(x->r - x->l, false);
            t.final_CBT.assign(t.final_CBT.size(), {-INF, {}});
            if (x->r - x->l == 1) {
                int curs = (data[x->l] > 0) ? 1 : 0;
                double value = (data[x->l] > 0) ? data[x->l] : -data[x->l];
                t.final_CBT.assign(t.final_CBT.size(), {-INF, {}});
                to_vector(curs, to_push);
                t.final_CBT[x->rules[curs].val] = {value, to_push};
                to_push.assign(to_push.size(), false);
                if (x->mt[0][0]) {
                    to_vector(1 ^ curs, to_push);
                    t.final_CBT[x->rules[1 ^ curs].val] = {-value, to_push};
                }
                return 0;
            }
            Gray(t, data, comps, adds, true);
        }
        return t.best[b];
    } else {
        size_t lCBT_size = (1ll << x->fir->fourth_sz);
        size_t rCBT_size = (1ll << x->sec->fourth_sz);
        size_t max_ind_r = std::min(rCBT_size, cnt1 + depth * cnt2);
        size_t to_count_l = std::min(lCBT_size, cnt1 + depth * cnt2);
        size_t ind_l = 0;
        while (max_ind_r != 0) {
            size_t cnt = 0;
            t.unused_right.clear();
            size_t cos_l = main_decode2(r, index * 2, data, comps, adds, ind_l, depth + 1, cnt1, cnt2);
            if (cos_l == -1)
                break;
            if (t.left_index_mapping.find(cos_l) == t.left_index_mapping.end()) {
                t.left_index_mapping[cos_l] = t.left_neighbors.size();
                t.left_neighbors.push_back(x->to_mul_left(cos_l));
            }
            auto rules_l = t.left_neighbors[t.left_index_mapping[cos_l]];
            for (auto &i : rules_l) {
                if (t.index_mapping.find(i.second) == t.index_mapping.end()) {
                    t.unused_right[i.first] = i.second;
                    cnt++;
                }
            }
            if (cnt == 0) {
                ind_l++;
                if (ind_l == to_count_l)
                    break;
                continue;
            }
            bool f = false;
            for (size_t j = 0; j < max_ind_r; j++) {
                auto t_neighbor = main_decode2(r, index * 2 + 1, data, comps, adds, j, depth + 1, cnt1, cnt2);
                if (t_neighbor == -1)
                    break;
                if (t.unused_right.find(t_neighbor) != t.unused_right.end()) {
                    t.interested_pairs.insert(branch(cos_l, t_neighbor, t.unused_right[t_neighbor]));
                    max_ind_r = j;
                    f = true;
                    break;
                }
            }
            ind_l++;
            if (ind_l == to_count_l)
                break;
        }
//        size_t max_ind_l = lCBT_size;
//        size_t ind_r = 0;
//        while (max_ind_l != 0) {
//            int cnt = 0;
//            size_t cos_r = main_decode2(r, index * 2 + 1, data, comps, adds, ind_r, depth + 1);
//            t.unused_left.clear();
//            if (t.right_index_mapping.find(cos_r) == t.right_index_mapping.end()) {
//                t.right_index_mapping[cos_r] = t.right_neighbors.size();
//                t.right_neighbors.push_back(x->to_mul_right(cos_r));
//            }
//            auto rules_r = t.right_neighbors[t.right_index_mapping[cos_r]];
//            for (auto &i : rules_r) {
//                if (t.index_mapping.find(i.second) == t.index_mapping.end()) {
//                    t.unused_left[i.first] = i.second;
//                    cnt++;
//                }
//            }
//            if (cnt == 0) {
//                ind_r++;
//                if (t.interested_pairs.size() != 0)
//                    break;
//                if (ind_r == rCBT_size)
//                    break;
//                continue;
//            }
//            bool f = false;
//            for (size_t j = 0; j < max_ind_l; j++) {
//                auto t_neighbor = main_decode2(r, index * 2, data, comps, adds, j, depth + 1);
//                if (t.unused_left.find(t_neighbor) != t.unused_left.end()) {
//                    t.interested_pairs.insert(branch(t_neighbor, cos_r, t.unused_left[t_neighbor]));
//                    max_ind_l = j;
//                    f = true;
//                    break;
//                }
//            }
//            ind_r++;
//            if (ind_r == rCBT_size)
//                break;
//            if (!f)
//                break;
//        }
        //assert(x->interested_pairs.size() != 0);
        if (t.interested_pairs.size() == 0)
            return -1;
        branch best_pair;
        double mx = -1000000;
        for (auto it : t.interested_pairs) {
            size_t real_left_index = r[index * 2].index_mapping[it.ind_l];
            size_t real_right_index = r[index * 2 + 1].index_mapping[it.ind_r];
            double val = r[index * 2].final_CBT[real_left_index].first +
                         r[index * 2 + 1].final_CBT[real_right_index].first;
            adds++;
            comps++;
            if (val > mx) {
                mx = val;
                best_pair = it;
            }
        }
        comps--;
        t.best.push_back(best_pair.val);
        t.index_mapping[best_pair.val] = t.final_CBT.size();
        size_t real_left_index = r[index * 2].index_mapping[best_pair.ind_l];
        size_t real_right_index = r[index * 2 + 1].index_mapping[best_pair.ind_r];
        t.final_CBT.push_back({mx, concat(r[index * 2].final_CBT[real_left_index].second,
                                          r[index * 2 + 1].final_CBT[real_right_index].second)});
        t.interested_pairs.clear();
        return best_pair.val;
    }
}

void clear_structures(std::vector<CBT> &t, size_t index) {
    t[index].best.clear();
    t[index].final_CBT.clear();
    t[index].interested_pairs.clear();
    t[index].index_mapping.clear();
    if (t[index].x->is_leaf)
        return;
    clear_structures(t, index * 2);
    clear_structures(t, index * 2 + 1);
}

std::vector<bool> decode2(std::vector<CBT> &t, const std::vector<double> &data, long long &comps, long long &adds, size_t cnt1, size_t cnt2) {
    clear_structures(t, 1);
    auto ind = main_decode2(t, 1, data, comps, adds, 0, 0, cnt1, cnt2);
    if (ind == -1)
        return std::vector<bool>(t[1].x->r - t[1].x->l, false);
    return t[1].final_CBT[ind].second;
}

void fill_vector_pMatrix(std::vector<CBT> &t, pMatrix *x, size_t ind) {
    t[ind].x = x;
    if (x->is_leaf)
        return;
    fill_vector_pMatrix(t, x->fir, ind * 2);
    fill_vector_pMatrix(t, x->sec, ind * 2 + 1);
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
    std::vector<CBT> all_CBTs((1 << m) * 4);
    fill_vector_pMatrix(all_CBTs, ptr, 1);
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
            auto recieved = (f) ? decode2(all_CBTs, x, comps, adds, 100, 100) : decode(all_CBTs, 1, x, comps, adds);
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

void one_thread_work(std::vector<CBT> &q, std::vector<std::vector<bool>> &before_noise,
                     const std::vector<std::vector<double>> &codewords, long long &adds, long long &comps,
                     size_t &good, size_t cnt1, size_t cnt2) {
    for (size_t i = 0; i < codewords.size(); i++) {
        auto recieve = decode2(q, codewords[i], comps, adds, cnt1, cnt2);
        if (i % 2 == 0)
            std::cout << i << " " << (adds + comps) / (i + 1) << "\n";
        good += cmp(recieve, before_noise[i]);
    }
}

void check_polar(size_t n, size_t k, bool f, size_t cnt_iter, size_t threads_cnt, size_t cnt1, size_t cnt2) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::vector<long long> comps(threads_cnt), adds(threads_cnt);
    std::vector<size_t> goods(threads_cnt);
    std::vector<std::vector<std::vector<bool>>> before_coded(threads_cnt);
    std::vector<std::vector<std::vector<double>>> with_noise(threads_cnt);
    PolarEncoder q = PolarEncoder(n, k);
    std::vector<std::thread> threads;
//    auto action = [](std::vector<CBT> q,
//                     std::vector<std::vector<bool>> before_noise, const std::vector<std::vector<double>> codewords,
//                     long long &adds1, long long &comps1,
//                     size_t &goodd) {
//        return one_thread_work(q, before_noise, codewords, comps1, adds1, goodd);
//    };
    for (double Eb_N0_dB = 1.0; Eb_N0_dB <= 5.0; Eb_N0_dB += 0.5) {
        std::vector<std::vector<CBT>> all_CBT(threads_cnt, std::vector<CBT>(4 * n));
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
        for (size_t i = 0; i < threads_cnt; i++) {
            fill_vector_pMatrix(all_CBT[i], ptr, 1);
            std::fill(comps.begin(), comps.end(), 0);
            std::fill(adds.begin(), adds.end(), 0);
            std::fill(goods.begin(), goods.end(), 0);
            before_coded[i].clear();
            with_noise[i].clear();
        }
        threads.clear();

        std::cout << "Created \n";
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

            before_coded[i % threads_cnt].push_back(coded);
            with_noise[i % threads_cnt].push_back(x);
        }
        for (size_t i = 0; i < threads_cnt; i++)
            threads.push_back(std::thread(one_thread_work, std::ref(all_CBT[i]), std::ref(before_coded[i]), std::cref(with_noise[i]),
                              std::ref(adds[i]), std::ref(comps[i]), std::ref(goods[i]), std::ref(cnt1), std::ref(cnt2)));
        std::cout.precision(7);
        for (size_t i = 0; i < threads.size(); i++)
            threads[i].join();
        long long total_good = 0;
        long long total_ops = 0;
        for (size_t i = 0; i < threads.size(); i++) {
            total_good += goods[i];
            total_ops += comps[i] + adds[i];
        }
        std::cout << std::fixed << (double) total_good / cnt_iter << ' ' << total_ops / cnt_iter << "\n";

    }

}

int main() {
    srand(time(NULL));
//    check(1, 3, true);
//    check(2, 5, true);
//    check(2, 6, true);
//    check(3, 6, true);
//    check(3, 6, true);
//    check_polar(256, 128, true, 1000, 8);
    check_polar(512, 256, true, 24, 8, 50, 100);
//    check_polar(1024, 512, true, 0);

    return 0;
}
