#include "TalVardyListDecoder.h"
#include "polar_encoder.h"


int main() {
    srand(time(0));
    int N = 16;
    int log_N = 4;
    int K = 8;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    for (double Eb_N0_dB = -0.0; Eb_N0_dB <= 0.0; Eb_N0_dB += 1.) {
        double sigma_square = 0.5 * ((double) N / K) * ((double) pow(10.0, -Eb_N0_dB / 10));
        std::normal_distribution<> d{0, sqrt(sigma_square)};
        PolarEncoder t(N, K, sqrt(sigma_square));
        TalVardyListDecoder decoder(log_N, K, 1, sqrt(sigma_square), t.frozen);
        for (size_t i = 0; i < 1; i++) {
            std::vector<bool> word = {0, 1, 1, 0, 0, 0, 0, 0 };
            std::vector<bool> coded = {0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0};
            std::cout << "CODED: ";
            for (size_t j = 0; j < coded.size(); j++)
                std::cout << coded[j] << ' ';
            std::cout << "\n--------------\n";
            std::vector<double> noise;
            noise.reserve(coded.size());
            for (size_t j = 0; j < coded.size(); j++)
                noise.push_back(d(gen));
            auto x = add_noise(coded, noise);
            auto decoded = decoder.decode(x);
            std::cout << "DECODED: ";
            for (size_t j = 0; j < decoded.size(); j++)
                std::cout << decoded[j] << ' ';
            std::cout << "\n--------------\n";
        }
    }
    return 0;
}