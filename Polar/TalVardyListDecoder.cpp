#include <iostream>
#include "TalVardyListDecoder.h"

void TalVardyListDecoder::init_data_structures(bool f) {
    while (!inactivePathIndices.empty())
        inactivePathIndices.pop();
    activePath.resize(L);
    arrayPointer_P.resize(M + 1);
    arrayPointer_C.resize(M + 1);
    pathIndexToArrayIndex.resize(M + 1);
    arrayReferenceCount.resize(M + 1);
    for (size_t i = 0; i < M + 1; i++) {
        arrayPointer_P[i].resize(L);
        arrayPointer_C[i].resize(L);
        pathIndexToArrayIndex[i].resize(L);
        arrayReferenceCount[i].resize(L);
    }
    inactiveArrayIndices.resize(M + 1);
    for (size_t i = 0; i < M + 1; i++) {
        while (!inactiveArrayIndices[i].empty())
            inactiveArrayIndices[i].pop();
    }
    for (size_t lam = 0; lam <= M; lam++) {
        for (size_t l = 0; l < L; l++) {
            size_t length = 1 << (M - lam);
            if (f) {
                arrayPointer_P[lam][l] = new double[2 * length];
                arrayPointer_C[lam][l] = new uint8_t[2 * length];
            }
            for (size_t i = 0; i < 2 * length; i++) {
                arrayPointer_P[lam][l][i] = 0;
                arrayPointer_C[lam][l][i] = 0;
            }
            arrayReferenceCount[lam][l] = 0;
            inactiveArrayIndices[lam].push(l);
        }
    }

    for (size_t l = 0; l < L; l++) {
        activePath[l] = false;
        inactivePathIndices.push(l);
    }


}

uint16_t TalVardyListDecoder::assignInitialPath() {
    auto l = inactivePathIndices.top();
    inactivePathIndices.pop();
    activePath[l] = true;
    for (size_t lam = 0; lam <= M; lam++) {
        auto s = inactiveArrayIndices[lam].top();
        inactiveArrayIndices[lam].pop();
        pathIndexToArrayIndex[lam][l] = s;
        arrayReferenceCount[lam][s] = 1;
    }
    return l;
}

uint16_t TalVardyListDecoder::clonePath(uint16_t l) {
    auto new_l = inactivePathIndices.top();
    inactivePathIndices.pop();
    activePath[new_l] = true;
    for (size_t lam = 0; lam <= M; lam++) {
        auto s = pathIndexToArrayIndex[lam][l];
        pathIndexToArrayIndex[lam][new_l] = s;
        arrayReferenceCount[lam][s]++;
    }
    return new_l;
}

void TalVardyListDecoder::killPath(uint16_t l) {
    activePath[l] = false;
    inactivePathIndices.push(l);
    for (size_t lam = 0; lam <= M; lam++) {
        auto s = pathIndexToArrayIndex[lam][l];
        arrayReferenceCount[lam][s]--;
        if (!arrayReferenceCount[lam][s])
            inactiveArrayIndices[lam].push(s);
    }
}

double *TalVardyListDecoder::getArrayPointer_P(uint16_t lam, uint16_t l) {
    auto s = pathIndexToArrayIndex[lam][l];
    auto new_s = s;
    if (arrayReferenceCount[lam][s] != 1) {
        new_s = inactiveArrayIndices[lam].top();
        inactiveArrayIndices[lam].pop();
        size_t length = (1 << (M - lam + 1));
        std::copy(arrayPointer_P[lam][s], arrayPointer_P[lam][s] + length, arrayPointer_P[lam][new_s]);
        std::copy(arrayPointer_C[lam][s], arrayPointer_C[lam][s] + length, arrayPointer_C[lam][new_s]);
        arrayReferenceCount[lam][s]--;
        arrayReferenceCount[lam][new_s] = 1;
        pathIndexToArrayIndex[lam][l] = new_s;
    }
    return arrayPointer_P[lam][new_s];
}

uint8_t *TalVardyListDecoder::getArrayPointer_C(uint16_t lam, uint16_t l) {
    auto s = pathIndexToArrayIndex[lam][l];
    auto new_s = s;
    if (arrayReferenceCount[lam][s] != 1) {
        new_s = inactiveArrayIndices[lam].top();
        inactiveArrayIndices[lam].pop();
        size_t length = (1 << (M - lam + 1));
        std::copy(arrayPointer_P[lam][s], arrayPointer_P[lam][s] + length, arrayPointer_P[lam][new_s]);
        std::copy(arrayPointer_C[lam][s], arrayPointer_C[lam][s] + length, arrayPointer_C[lam][new_s]);
        arrayReferenceCount[lam][s]--;
        arrayReferenceCount[lam][new_s] = 1;
        pathIndexToArrayIndex[lam][l] = new_s;
    }
    return arrayPointer_C[lam][new_s];
}

void TalVardyListDecoder::recursivelyCalcP(uint16_t lam, uint16_t phi) {
    if (lam == 0)
        return;
    uint16_t psi = (phi >> 1);
    if (phi % 2 == 0)
        recursivelyCalcP(lam - 1, psi);
    double sigma = 0;
    for (size_t l = 0; l < L; l++) {
        if (!activePath[l])
            continue;
        auto *pLambda = getArrayPointer_P(lam, l);
        auto *pLambdaPred = getArrayPointer_P(lam - 1, l);
        auto *cLambda = getArrayPointer_C(lam, l);
        size_t sz = (1 << (M - lam));
        for (size_t b = 0; b < sz; b++) {
            if (phi % 2 == 0) {
                pLambda[2 * b] = (pLambdaPred[2 * b] * pLambdaPred[2 * (b + sz)] +
                                  pLambdaPred[2 * b + 1] * pLambdaPred[2 * (b + sz) + 1]);
                pLambda[2 * b + 1] = (pLambdaPred[b * 2 + 1] * pLambdaPred[2 * (b + sz)] +
                                      pLambdaPred[b * 2] * pLambdaPred[2 * (b + sz) + 1]);
            } else {
                auto u1 = cLambda[2 * b];
                pLambda[2 * b] = (pLambdaPred[b * 2 + u1] * pLambdaPred[2 * (b + sz)]);
                pLambda[2 * b + 1] = (pLambdaPred[b * 2 + (1 - u1)] * pLambdaPred[2 * (b + sz) + 1]);
            }
            sigma = std::max(sigma, std::max(pLambda[2 * b], pLambda[2 * b + 1]));
        }
    }

    size_t len = (1 << (M - lam));
    for (size_t l = 0; l < L; l++) {
        if (!activePath[l])
            continue;
        auto *pLambda = getArrayPointer_P(lam, l);
        for (size_t bet = 0; bet < len; bet++) {
            pLambda[2 * bet] /= sigma;
            pLambda[2 * bet + 1] /= sigma;
        }
//        for (size_t i = 0; i < len; i++)
//            std::cout << std::fixed << pLambda[2 * i] << '-' << pLambda[2 * i + 1] << ' ';
//        std::cout << "\n";
    }
}

void TalVardyListDecoder::recursivelyUpdateC(uint16_t lam, uint16_t phi) {
    uint16_t psi = (phi >> 1);
    size_t sz = (1 << (M - lam));
    for (size_t l = 0; l < L; l++) {
        if (!activePath[l])
            continue;
        auto *cLambda = getArrayPointer_C(lam, l);
        auto *cLambdaPred = getArrayPointer_C(lam - 1, l);
        for (size_t b = 0; b < sz; b++) {
            cLambdaPred[2 * b + psi % 2] = (cLambda[2 * b] ^ cLambda[2 * b + 1]);
            cLambdaPred[2 * (b + sz) + psi % 2] = cLambda[2 * b + 1];
        }
    }
    if (psi % 2 == 1)
        recursivelyUpdateC(lam - 1, psi);
}

void TalVardyListDecoder::continuePaths_UnfrozenBit(uint16_t phi) {
    std::vector<double> probForks(2 * L, 0);
    std::vector<std::pair<double, size_t>> probs;
    std::vector<uint8_t> contForks(2 * L, 0);
    int i = 0;
    for (size_t l = 0; l < L; l++) {
        if (activePath[l]) {
            auto *pLambda = getArrayPointer_P(M, l);
            probForks[2 * l] = pLambda[0];
            probForks[2 * l + 1] = pLambda[1];
            i++;
        } else {
            probForks[2 * l] = NAN;
            probForks[2 * l + 1] = NAN;
        }
        probs.emplace_back(probForks[2 * l], 2 * l);
        probs.emplace_back(probForks[2 * l + 1], 2 * l + 1);
    }
    int ro = std::min(2 * i, L);
    std::sort(probs.begin(), probs.end());
    std::reverse(probs.begin(), probs.end());

    for (size_t j = 0; j < ro; j++) {
        contForks[probs[j].second] = 1;
    }

    for (size_t l = 0; l < L; l++) {
        if (!activePath[l])
            continue;
        if (!contForks[2 * l] && !contForks[2 * l + 1])
            killPath(l);
    }

    for (size_t l = 0; l < L; l++) {
        if (!contForks[2 * l] && !contForks[2 * l + 1])
            continue;
        auto *cLambda = getArrayPointer_C(M, l);
        if (contForks[2 * l] && contForks[2 * l + 1]) {
            cLambda[phi % 2] = 0;
            auto new_l = clonePath(l);
            cLambda = getArrayPointer_C(M, new_l);
            cLambda[phi % 2] = 1;
        } else {
            if (contForks[2 * l])
                cLambda[phi % 2] = 0;
            else
                cLambda[phi % 2] = 1;
        }
    }
}

void TalVardyListDecoder::continuePaths_FrozenBit(uint16_t phi) {
    for (size_t l = 0; l < L; l++) {
        if (!activePath[l])
            continue;
        auto *cLambda = getArrayPointer_C(M, l);
        cLambda[phi % 2] = 0;
    }

}

std::vector<bool> TalVardyListDecoder::decode(const std::vector<double> &word) {
    init_data_structures(false);
    size_t l = assignInitialPath();
    double *p0 = getArrayPointer_P(0, l);
    for (size_t i = 0; i < word.size(); i++) {
        p0[2 * i] = getProb(word[i], dist, -1);
        p0[2 * i + 1] = getProb(word[i], dist, 1);
    }
    int uk = 0;
    for (size_t phi = 0; phi < word.size(); phi++) {
        recursivelyCalcP(M, phi);
        if (uk == frozen.size() || frozen[uk] != phi) {
            continuePaths_UnfrozenBit(phi);
        } else {
            uk++;
            continuePaths_FrozenBit(phi);
        }
        if (phi % 2)
            recursivelyUpdateC(M, phi);
    }
    size_t l1 = 0;
    double p1 = 0;
    for (size_t l = 0; l < L; l++) {
        if (!activePath[l])
            continue;
        auto *cm = getArrayPointer_C(M, l);
        auto *pm = getArrayPointer_P(M, l);
        if (p1 < pm[cm[1]]) {
            l1 = l;
            p1 = pm[cm[1]];
        }
    }
    auto *c0 = getArrayPointer_C(0, l1);
    std::vector<bool> ans;
    for (size_t i = 0; i < (1 << M); i++)
        ans.push_back(c0[2 * i]);

    return ans;
}

double TalVardyListDecoder::getProb(double x, double dispersion, double expValue) {
    return 1 / (sqrt(2 * M_PI) * dispersion) * exp(-pow((x - expValue) / dispersion, 2) / 2);
}
