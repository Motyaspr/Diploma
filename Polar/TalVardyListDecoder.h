#ifndef DIPLOMA_TALVARDYLISTDECODER_H
#define DIPLOMA_TALVARDYLISTDECODER_H

#include "common.h"
#include <stack>

class TalVardyListDecoder {
    int M, K, L;

    double eps;

    std::vector<size_t> frozen;

    std::stack<uint16_t> inactivePathIndices;
    std::vector<bool> activePath;

    std::vector<std::vector<double *>> arrayPointer_P;
    std::vector<std::vector<uint8_t *>> arrayPointer_C;
    std::vector<std::vector<uint16_t>> pathIndexToArrayIndex;

    std::vector<std::stack<uint16_t>> inactiveArrayIndices;
    std::vector<std::vector<uint16_t>> arrayReferenceCount;

    void init_data_structures();

    TalVardyListDecoder(int _M, int _K, int _L, double _eps, std::vector<bool> _frozen) :
            M(_M), K(_K), L(_L), eps(_eps) {
        frozen.resize(M - K);
        for (size_t i = 0; i < _frozen.size(); i++)
            if (_frozen[i])
                frozen.push_back(i);
        init_data_structures();
    }

    uint16_t assignInitialPath();

    uint16_t clonePath(uint16_t l);

    void killPath(uint16_t l);

    double *getArrayPointer_P(uint16_t lam, uint16_t l);

    uint8_t *getArrayPointer_C(uint16_t lam, uint16_t l);

    void recursivelyCalcP(uint16_t lam, uint16_t phi);

    void recursivelyUpdateC(uint16_t lam, uint16_t phi);

    void continuePaths_FrozenBit(uint16_t phi);

    void continuePaths_UnfrozenBit(uint16_t phi);

    std::vector<bool> decode(const std::vector<bool> &word);

};


#endif //DIPLOMA_TALVARDYLISTDECODER_H
