#include "TalVardyListDecoder.h"
#include "polar_encoder.h"


int main() {
    srand(time(0));
    int N = 1024;
    int log_N = 10;
    int K = 512;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    for (size_t l = 16; l <= 16; l *= 2) {
        std::cout << "Polar Code(" << N << ", " << K << ", " << l << ")\n";
        for (double Eb_N0_dB = 2.5; Eb_N0_dB <= 3.0; Eb_N0_dB += 0.5) {
            double sigma_square = 0.5 * ((double) N / K) * ((double) pow(10.0, -Eb_N0_dB / 10));
            std::normal_distribution<> d{0, sqrt(sigma_square)};
            PolarEncoder t(N, K, sqrt(sigma_square));
            int cnt = 0;
            TalVardyListDecoder decoder(log_N, K, l, sqrt(sigma_square), t.frozen);
            for (size_t i = 0; i < ITER; i++) {
                std::vector<bool> word = gen_rand_vect(N);
                std::vector<bool> coded = t.encode(word);
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