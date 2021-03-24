#include "common.h"
#include "reed_muller.h"

const int MIN_MAKECBT = 2;

struct branch {
    long long ind_l, ind_r, val;
};

struct pMatrix {
    bool is_leaf;
    pMatrix *fir, *sec;
    int l, r, third_sz, fourth_sz;
    std::vector<std::vector<bool>> mt;
    std::vector<branch> rules;
    std::vector<std::pair<double, unsigned long long>> CBT;

    pMatrix(int _l, int _r) {
        l = _l;
        r = _r;
        is_leaf = ((r - l) <= MIN_MAKECBT);
        fir = nullptr;
        sec = nullptr;
        third_sz = 0;
        fourth_sz = 0;
    }

    void combCBT(const matrix &generator, const int &mid) {
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
//        std::cout << "KEKer: " << third_sz << ' ' << fourth_sz << "\n";
        std::sort(t.begin(), t.end());
        std::vector<std::vector<bool>> x;
        x.reserve(t.size());
        for (auto &i : t)
            x.emplace_back(std::move(i.second));
        mt = std::move(prepare_matrix(x, fourth_sz));
    }

    void makeCBT(const matrix &generator) {
        is_leaf = true;
        combCBT(generator, r);
    }

    void prepareLengthOneOrTwo() {
        if (r - l == 1) {
            rules.resize(2, {-1, -1, 0});
            CBT.resize(2);
            if (mt[0][0] == 1)
                rules[1] = {-1, -1, 1};
            return;
        }
        long long cosets_count = (1ll << (fourth_sz));
        long long masks_count = (1ll << (third_sz + fourth_sz));
        CBT.resize(cosets_count);
        rules.resize(1ll << (r - l));
        for (size_t i = 0; i < masks_count; i++) {
            auto x = get_ind(mt, i, r, false);
            rules[x] = {-1, -1, (long long) i % masks_count};
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
//        std::cout << "SYSTEMS SOLVED\n";
        long long cosets_count = (1ll << (fourth_sz));
        long long masks_count = (1ll << (third_sz + fourth_sz));
        CBT.resize(cosets_count, {-INF, 0});
        for (size_t i = 0; i < masks_count; i++) {
            long long fir_ind = get_ind(left_system_solutions, i, fir->fourth_sz, true);
            long long sec_ind = get_ind(right_system_solutions, i, sec->fourth_sz, true);
            long long rule = i % cosets_count;
            rules.push_back({fir_ind, sec_ind, rule});
        }
    }

};

void run(pMatrix *x, const matrix &generator) {
    if (x->r - x->l <= MIN_MAKECBT) {
        x->makeCBT(generator);
        x->prepareLengthOneOrTwo();
//        x->printMt();
        return;
    }
    int mid = (x->l + x->r) / 2;
    x->fir = new pMatrix(x->l, mid);
    x->sec = new pMatrix(mid, x->r);
    run(x->fir, generator);
    run(x->sec, generator);
    x->combCBT(generator, mid);
//    x->printMt();
    x->merge();
}


void Gray(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds) {
    double total_sum = data[x->l];
    for (size_t i = x->l + 1; i < x->r; i++) {
        total_sum += data[i];
        adds++;
    }

}

unsigned long long decode(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds) {
    if (x->is_leaf) {
        if (x->r - x->l == 1) {
            int curs = (data[x->l] > 0) ? 1 : 0;
            double value = (data[x->l] > 0) ? data[x->l] : -data[x->l];
            x->CBT.assign(x->CBT.size(), {0, 0});
            x->CBT[x->rules[curs].val] = {value, curs};
            if (x->mt[0][0]) {
                x->CBT[(x->rules[1 ^ curs].val)] = {-value, 1 ^ curs};
            }
            return 0;
        }

    } else {
        int mid = (x->l + x->r) / 2;
        x->CBT.assign(x->CBT.size(), {0, 0});
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

void check(int r, int m) {
    std::random_device rd{};
    std::mt19937 gen{rd()};
    ReedMuller reedMuller(r, m);
    matrix t(reedMuller.generated);
    auto *ptr = new pMatrix(0, t.m);
    t.to_span();
    run(ptr, t);
    std::cout << "RM(" << r << ", " << m << ") created\n";
    long long comps = 0, adds = 0;
    for (double Eb_N0_dB = 0.0; Eb_N0_dB <= 6.0; Eb_N0_dB += 1.) {
        int cnt = 0;
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
            auto ans = decode(ptr, x, comps, adds);
            to_vector(ans, recieve);
            auto decoded = get_message(t, recieve);
            cnt += cmp(decoded, word);
        }
        std::cout.precision(7);
        std::cout << std::fixed << (int) Eb_N0_dB << ' ' << (double) cnt / ITER << "\n";
    }
    std::cout << "Count of adds and cmps:" << (comps + adds) / (ITER * 7) << "\n";
    std::cout << "Count of adds:" << (adds) / (ITER * 7) << "\n";
    std::cout << "Count of cmps:" << (comps) / (ITER * 7) << "\n";
    delete ptr;
}

int main() {
    srand(time(NULL));
    check(1, 3);
    check(2, 5);
    check(2, 6);
    check(3, 6);
    check(4, 6);
    return 0;
}

