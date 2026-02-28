#include <cassert>
#include "openfhe.h"
#include "encoding/z-encoding.h"
#include "scheme/ckksrns/z-leveledshe.h"
#include "scheme/ckksrns/z-advancedshe.h"

using namespace lbcrypto;
using CiphertextT        = ConstCiphertext<DCRTPoly>;
using MutableCiphertextT = Ciphertext<DCRTPoly>;
using CCParamsT          = CCParams<CryptoContextBFVRNS>;
using CryptoContextT     = CryptoContext<DCRTPoly>;
using EvalKeyT           = EvalKey<DCRTPoly>;
using PlaintextT         = Plaintext;
using PrivateKeyT        = PrivateKey<DCRTPoly>;
using PublicKeyT         = PublicKey<DCRTPoly>;

//=============================================================================
// Encryption Utils
//=============================================================================

#define WARMUP(x, str)                                                                     \
    {                                                                                      \
        auto startW = std::chrono::high_resolution_clock::now();                           \
        (x);                                                                               \
        auto endW                           = std::chrono::high_resolution_clock::now();   \
        std::chrono::duration<double> diffW = endW - startW;                               \
        std::cout << "Finished Warmup for " str " " << diffW.count() << " s" << std::endl; \
    }                                                                                      \
    while (0)

#define BENCHMARK(x, times, str)                                                                            \
    {                                                                                                       \
        auto start = std::chrono::high_resolution_clock::now();                                             \
        for (size_t i = 0; i != (times); ++i) {                                                             \
            (x);                                                                                            \
        }                                                                                                   \
        auto end                           = std::chrono::high_resolution_clock::now();                     \
        std::chrono::duration<double> diff = end - start;                                                   \
        if (times > 0) {                                                                                    \
            std::cout << "Time for " str " : " << diff.count() / (times) << " s" << std::endl << std::endl; \
        }                                                                                                   \
    }                                                                                                       \
    while (0)

