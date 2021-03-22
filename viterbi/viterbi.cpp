//
// Created by Matvey Sprikut on 08.11.2020.
//

#include "common.h"
#include "reed_muller.h"


std::vector<bool>
decode(std::vector<double> &code, std::vector<std::vector<vertex>> trellis) {
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
        for (size_t j = 0; j < dp[i].size(); j++) {
            double w = (trellis[i][j].zer.w == 0) ? -code[i] : code[i];
            if (dp[i + 1][trellis[i][j].zer.nxt] > dp[i][j] + w) {
                dp[i + 1][trellis[i][j].zer.nxt] = dp[i][j] + w;
                p[i + 1][trellis[i][j].zer.nxt] = j;
            }
            if (trellis[i][j].one.w != -1) {
                w = (trellis[i][j].one.w == 0) ? -code[i] : code[i];
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
    //std::cout << dp.back().back() << "\n";
    return ans;
}


int main() {
    srand(time(NULL));
    std::random_device rd{};
    std::mt19937 gen{rd()};

    ReedMuller reedMuller(1, 3);
    matrix t(reedMuller.generated);
    t.print();
    std::cout << "\n";
    t.to_span();
    t.print();

    std::vector<std::vector<vertex>> trellis = create_graph(t);

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
            auto recieve = decode(x, trellis);
            auto decoded = get_message(t, recieve);
            cnt += cmp(decoded, word);
        }
        std::cout.precision(7);
        std::cout << std::fixed << (int) Eb_N0_dB << ' ' << (double) cnt / ITER << "\n";
    }

}


