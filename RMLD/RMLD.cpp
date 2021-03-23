#include "common.h"
#include "reed_muller.h"

struct branch {
    int ind_l, ind_r, val;

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
        is_leaf = ((r - l) == 1);
        fir = nullptr;
        sec = nullptr;
        third_sz = 0;
        fourth_sz = 0;
    }

    void combCBT(const matrix &generator) {
        int mid = (l + r) / 2;
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
        bool f = false;
        mt.resize(1, std::vector<bool>(1, false));
        for (size_t i = 0; i < generator.n; i++)
            if (generator.mt[i][l] == 1) {
                f = true;
                break;
            }
        if (f) {
            mt[0][0] = true;
        }
        fourth_sz = 1;
    }

    void prepareLengthOne() {
        rules.resize(2, {-1, -1, 0});
        CBT.resize(2);
        if (mt[0][0] == 1)
            rules[1] = {-1, -1, 1};
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
        for (size_t i = 0; i < masks_count; i++) {
            int fir_ind = get_ind(left_system_solutions, i, fir->fourth_sz);
            int sec_ind = get_ind(right_system_solutions, i, sec->fourth_sz);
            int rule = i % cosets_count;
            rules.push_back({fir_ind, sec_ind, rule});
        }
    }

};

void run(pMatrix *x, const matrix &generator) {
    if (x->r - x->l == 1) {
        x->makeCBT(generator);
        x->prepareLengthOne();
//        x->printMt();
        return;
    }
    int mid = (x->l + x->r) / 2;
    x->fir = new pMatrix(x->l, mid);
    x->sec = new pMatrix(mid, x->r);
    run(x->fir, generator);
    run(x->sec, generator);
    x->combCBT(generator);
//    x->printMt();
    x->merge();
}

unsigned long long decode(pMatrix *x, const std::vector<double> &data, long long &comps, long long &adds) {
    if (x->is_leaf) {
        int curs = (data[x->l] > 0) ? 1 : 0;
        double value = abs(data[x->l]);
        x->CBT.assign(x->CBT.size(), {-INF, 0});
        x->CBT[x->rules[curs].val] = {value, curs};
        if (x->mt[0][0]) {
            x->CBT[(x->rules[1 ^ curs].val)] = {-value, 1 ^ curs};
        }
    } else {
        int mid = (x->l + x->r) / 2;
        x->CBT.assign(x->CBT.size(), {-INF, 0});
        decode(x->sec, data, comps, adds);
        decode(x->fir, data, comps, adds);
        for (auto &ind : x->rules) {
            double sum = x->fir->CBT[ind.ind_l].first + x->sec->CBT[ind.ind_r].first;
            adds++;
            comps++;
            if (x->CBT[ind.val].first < sum) {
                x->CBT[ind.val] = {sum,
                                   x->fir->CBT[ind.ind_l].second + ((x->sec->CBT[ind.ind_r].second) << (mid - x->l))};
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
    run(ptr, t);
    std::cout << "created\n";
    long long comps = 0, adds = 0;
    for (double Eb_N0_dB = -0.0; Eb_N0_dB <= 6.0; Eb_N0_dB += 1.) {
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
    std::cout << std::fixed << "Count of adds and cmps:" << (double) (comps + adds) / (ITER * 7) << "\n";
    std::cout << std::fixed << "Count of adds:" << (double) (adds) / (ITER * 7) << "\n";
    std::cout << std::fixed << "Count of cmps:" << (double) (comps) / (ITER * 7) << "\n";
    delete ptr;
}

int main() {
    srand(time(NULL));
//    check(1, 3);
//    check(2, 5);
//    check(2, 6);
    check(3, 6);
    check(4, 6);
    return 0;
}

