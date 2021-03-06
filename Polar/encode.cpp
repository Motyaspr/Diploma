#include "TalVardyListDecoder.h"
#include "polar_encoder.h"
#include "common.h"


int main() {
    srand(time(0));
    int N = 512;
    int log_N = 9;
    int K = 256;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    PolarEncoder t(N, K);
    for (size_t l = 1; l <= 16; l *= 2) {
        std::cout << "Polar Code(" << N << ", " << K << ", " << l << ")\n";
        for (double Eb_N0_dB = 0.0; Eb_N0_dB <= 6.0; Eb_N0_dB += 0.5) {
            double sigma_square = 0.5 * ((double) N / K) * ((double) pow(10.0, -Eb_N0_dB / 10));
            std::normal_distribution<> d{0, sqrt(sigma_square)};
            t.reuse_frozen(sqrt(sigma_square));
            int cnt = 0;
            TalVardyListDecoder decoder(log_N, K, l, sqrt(sigma_square), t.frozen);
            for (size_t i = 0; i < ITER; i++) {
                std::vector<bool> word = gen_rand_vect(N);
                std::vector<bool> coded = t.encode(word, true);
                std::vector<double> noise;
                noise.reserve(coded.size());
                for (size_t j = 0; j < coded.size(); j++)
                    noise.push_back(d(gen));
                auto x = add_noise(coded, noise);
                auto decoded = decoder.decode(x);
                cnt += cmp(decoded, coded);
            }
            std::cout.precision(5);
            std::cout << std::fixed << Eb_N0_dB << ' ' << (double) cnt / ITER << "\n";
        }
        std::cout << "---------------\n";
    }
    return 0;
}