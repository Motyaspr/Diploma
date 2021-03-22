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
    std::vector<branch> cosets;
    std::vector<std::pair<double, unsigned long long>> cbt;

    pMatrix(int _l, int _r) {
        l = _l;
        r = _r;
        is_leaf = ((r - l) == 1);
        fir = nullptr;
        sec = nullptr;
        third_sz = 0;
        fourth_sz = 0;
    }

    void createBig(const matrix &generator) {
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

    void createSmall(const matrix &generator) {
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

    void smallCosets() {
        cosets.resize(2, {-1, -1, 0});
        cbt.resize(2);
        if (mt[0][0] == 1)
            cosets[1] = {-1, -1, 1};
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
        cbt.resize(cosets_count, {-INF, 0});
        for (size_t i = 0; i < masks_count; i++) {
            int fir_ind = get_ind(left_system_solutions, i, fir->fourth_sz);
            int sec_ind = get_ind(right_system_solutions, i, sec->fourth_sz);
            int coset = i % cosets_count;
            cosets.push_back({fir_ind, sec_ind, coset});
        }
    }

};

void run(pMatrix *x, const matrix &generator) {
    if (x->r - x->l == 1) {
        x->createSmall(generator);
        x->smallCosets();
//        x->printMt();
        return;
    }
    int mid = (x->l + x->r) / 2;
    x->fir = new pMatrix(x->l, mid);
    x->sec = new pMatrix(mid, x->r);
    run(x->fir, generator);
    run(x->sec, generator);
    x->createBig(generator);
//    x->printMt();
    x->merge();
}

unsigned long long decode(pMatrix *x, const std::vector<double> &data) {
    if (x->is_leaf) {
        int curs = (data[x->l] > 0) ? 1 : 0;
        double value = abs(data[x->l]);
        x->cbt.assign(x->cbt.size(), {-INF, 0});
        x->cbt[x->cosets[curs].val] = {value, curs};
        if (x->mt[0][0]) {
            x->cbt[(x->cosets[1 ^ curs].val)] = {-value, 1 ^ curs};
        }
    } else {
        int mid = (x->l + x->r) / 2;
        x->cbt.assign(x->cbt.size(), {-INF, 0});
        decode(x->sec, data);
        decode(x->fir, data);
        for (auto &ind : x->cosets) {
            double sum = x->fir->cbt[ind.ind_l].first + x->sec->cbt[ind.ind_r].first;
            if (x->cbt[ind.val].first < sum) {
                x->cbt[ind.val] = {sum,
                                   x->fir->cbt[ind.ind_l].second + ((x->sec->cbt[ind.ind_r].second) << (mid - x->l))};
            }
        }
    }
    return x->cbt[0].second;
}

int main() {
    srand(time(NULL));
    std::random_device rd{};
    std::mt19937 gen{rd()};
    ReedMuller reedMuller(1, 3);
    matrix t(reedMuller.generated);
    t.to_span();
    auto *ptr = new pMatrix(0, t.m);
    run(ptr, t);
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
            auto ans = decode(ptr, x);
            to_vector(ans, recieve);
            cnt += cmp(recieve, coded);
        }
        std::cout.precision(7);
        std::cout << std::fixed << (int) Eb_N0_dB << ' ' << (double) cnt / ITER << "\n";
    }
}
