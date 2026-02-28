#include "encoding/z-encoding.h"
#include "math/hermite.h"
#include "math/z-constants.h"
#include "scheme/ckksrns/z-fhe.h"
#include "math/dftransform.h"
#include "math/dftransform-bigcomplex.h"

namespace lbcrypto {

//------------------------------------------------------------------------------
// Bootstrap Wrapper
//------------------------------------------------------------------------------

void FHEZImpl::EvalBootstrapSetup(const CryptoContextImpl<DCRTPoly>& cc, uint32_t zN, uint32_t zSlots,
                                  std::vector<uint32_t> levelBudget, std::vector<uint32_t> dim1, uint32_t w,
                                  int32_t arithToBooleanCutoff, uint32_t lutOrder) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());

    uint32_t N  = cc.GetRingDimension();
    uint32_t M  = cc.GetCyclotomicOrder();
    auto cSlots = zN * zSlots / 2;

    // Initialize transforms
    ZLinearTransform::Initialize(zN);
    auto cSlots2 = cSlots * 2;
    // For regular encoding
    DiscreteFourierTransform::Initialize(cSlots * 4, cSlots);
    //DiscreteFourierTransformBigComplex::Initialize(cSlots * 4, cSlots);
    // For encoding of t and tInv, and scalar ptxt in Z
    DiscreteFourierTransform::Initialize(zN * 2, zN / 2);
    //DiscreteFourierTransformBigComplex::Initialize(zN * 2, zN / 2);
    // For encoding of bootstrapping related plaintext for sparse bootstrapping
    DiscreteFourierTransform::Initialize(cSlots2 * 4, cSlots2);
    //DiscreteFourierTransformBigComplex::Initialize(cSlots2 * 4, cSlots2);
    auto& ZU = ZLinearTransform::GetZU(zN);
    auto& ZV = ZLinearTransform::GetZUInverse(zN);

    m_bootPrecomMap[cSlots] = std::make_shared<ZBootstrapPrecom>();

    auto& precom     = m_bootPrecomMap[cSlots];
    precom->m_cSlots = cSlots;

    // even for the case of a single slot we need one level for rescaling
    uint32_t logSlots = (cSlots < 3) ? 1 : std::log2(cSlots);
    // Perform some checks on the level budget and compute parameters
    uint32_t newBudget0 = levelBudget[0];
    if (newBudget0 > logSlots) {
        std::cerr << "\nWarning, the level budget for encoding is too large. Setting it to " << logSlots << std::endl;
        newBudget0 = logSlots;
    }
    if (newBudget0 < 1) {
        std::cerr << "\nWarning, the level budget for encoding can not be zero. Setting it to 1" << std::endl;
        newBudget0 = 1;
    }
    uint32_t newBudget1 = levelBudget[1];
    if (newBudget1 > logSlots) {
        std::cerr << "\nWarning, the level budget for decoding is too large. Setting it to " << logSlots << std::endl;
        newBudget1 = logSlots;
    }
    if (newBudget1 < 1) {
        std::cerr << "\nWarning, the level budget for decoding can not be zero. Setting it to 1" << std::endl;
        newBudget1 = 1;
    }

    precom->m_paramsEnc = GetCollapsedFFTParams(cSlots, newBudget0, dim1[0]);
    precom->m_paramsDec = GetCollapsedFFTParams(cSlots, newBudget1, dim1[1]);

    uint32_t m     = 4 * cSlots;
    uint32_t mmask = m - 1;  // assumes m is power of 2
    bool isSparse  = (M != m);

    // store isSparse for convenience
    precom->m_isSparse = isSparse;

    // computes indices for all primitive roots of unity
    std::vector<uint32_t> rotGroup(cSlots);
    uint32_t fivePows = 1;
    for (uint32_t i = 0; i < cSlots; ++i) {
        rotGroup[i] = fivePows;
        fivePows *= 5;
        fivePows &= mmask;
    }

    // computes all powers of a primitive root of unity exp(2 * M_PI/m)
    std::vector<BigComplex> ksiPows(m + 1);
    ksiPows[0] = BigComplex(BigFixedPoint::one(), BigFixedPoint::zero());
    ksiPows[1] = R_ROOT_MAP.at(m);
    for (uint32_t j = 2; j < m; ++j) {
        ksiPows[j] = ksiPows[j - 1] * ksiPows[1];
    }
    ksiPows[m] = ksiPows[0];

    bool isLTBootstrap = (precom->m_paramsEnc.lvlb == 1) && (precom->m_paramsDec.lvlb == 1);

    auto I = BigComplex(BigFixedPoint::zero(), BigFixedPoint::one());

    // Scale R2C by 1/N so we don't have to normalize after R2C
    BigFixedPoint scaleC2R   = BigFixedPoint::one();  // No scaling for C2R
    BigFixedPoint scaleR2C   = BigFixedPoint::one();  // No scaling for R2C in high precision case
    BigFixedPoint scaleR2CN  = BigFixedPoint::one() / BigFixedPoint::positive(N);
    BigFixedPoint scaleR2CNK = scaleR2CN / BigFixedPoint::positive(K_SPARSE_ENCAPSULATED);

    // For FFT-like, scale can be distributed among levels
    auto r2clvlb = precom->m_paramsEnc.lvlb;
    std::vector<BigFixedPoint> scaleR2CFFT(r2clvlb, BigFixedPoint::one());
    std::vector<BigFixedPoint> scaleR2CFFTN;
    std::vector<BigFixedPoint> scaleR2CFFTNK;
    uint32_t log2N = std::round(std::log2(N));
    // We require K to be power of 2
    uint32_t log2K = std::round(std::log2(K_SPARSE_ENCAPSULATED));
    // distribute log2N among levels
    uint32_t baseScale = log2N / r2clvlb;
    uint32_t remScale  = log2N % r2clvlb;
    for (uint32_t i = 0; i < r2clvlb; ++i) {
        if (i != r2clvlb - 1) {
            auto scaleBFP = BigFixedPoint::positive(1l << (baseScale));
            scaleR2CFFTN.push_back(BigFixedPoint::one() / scaleBFP);
        }
        else {
            auto scaleBFP = BigFixedPoint::positive(1l << (baseScale + remScale));
            scaleR2CFFTN.push_back(BigFixedPoint::one() / scaleBFP);
        }
    }
    uint32_t baseScaleK = (log2N + log2K) / r2clvlb;
    uint32_t remScaleK  = (log2N + log2K) % r2clvlb;
    for (uint32_t i = 0; i < r2clvlb; ++i) {
        if (i != r2clvlb - 1) {
            scaleR2CFFTNK.push_back(BigFixedPoint::one() / BigFixedPoint::positive(1l << (baseScaleK)));
        }
        else {
            scaleR2CFFTNK.push_back(BigFixedPoint::one() / BigFixedPoint::positive(1l << (baseScaleK + remScaleK)));
        }
    }

    if (isLTBootstrap) {
        if (isSparse) {
            BigCMatrix U0(cSlots, BigCVector(cSlots));
            BigCMatrix U0hatT(cSlots, BigCVector(cSlots));
            BigCMatrix U1(cSlots, BigCVector(cSlots));
            BigCMatrix U1hatT(cSlots, BigCVector(cSlots));
            for (uint32_t i = 0; i < cSlots; ++i) {
                for (uint32_t j = 0; j < cSlots; ++j) {
                    U0[i][j]     = ksiPows[(j * rotGroup[i]) & mmask];
                    U0hatT[j][i] = U0[i][j].conj();
                    U1[i][j]     = I * U0[i][j];
                    U1hatT[j][i] = U1[i][j].conj();
                }
            }
            precom->m_U0Pre             = EvalLinearTransformPrecompute(cc, U0, U1, 1, scaleC2R);
            precom->m_U0hatTPre         = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0, scaleR2C);
            precom->m_U0hatTPreScaledN  = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0, scaleR2CN);
            precom->m_U0hatTPreScaledNK = EvalLinearTransformPrecompute(cc, U0hatT, U1hatT, 0, scaleR2CNK);
        }
        else {
            BigCMatrix U0(cSlots, BigCVector(cSlots));
            BigCMatrix U0hatT(cSlots, BigCVector(cSlots));
            for (uint32_t i = 0; i < cSlots; ++i) {
                for (uint32_t j = 0; j < cSlots; ++j) {
                    U0[i][j]     = ksiPows[(j * rotGroup[i]) & mmask];
                    U0hatT[j][i] = U0[i][j].conj();
                }
            }
            precom->m_U0Pre             = EvalLinearTransformPrecompute(cc, U0, scaleC2R);
            precom->m_U0hatTPre         = EvalLinearTransformPrecompute(cc, U0hatT, scaleR2C);
            precom->m_U0hatTPreScaledN  = EvalLinearTransformPrecompute(cc, U0hatT, scaleR2CN);
            precom->m_U0hatTPreScaledNK = EvalLinearTransformPrecompute(cc, U0hatT, scaleR2CNK);
        }
    }
    else {
        precom->m_U0PreFFT             = EvalSlotsToCoeffsPrecompute(cc, ksiPows, rotGroup, false);
        precom->m_U0hatTPreFFT         = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false, scaleR2CFFT);
        precom->m_U0hatTPreFFTScaledN  = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false, scaleR2CFFTN);
        precom->m_U0hatTPreFFTScaledNK = EvalCoeffsToSlotsPrecompute(cc, ksiPows, rotGroup, false, scaleR2CFFTNK);
    }

    // Now deal with Z Linear Transform
    // We note that ZUInverse is scaled.
    BigCMatrix ZU0(zN / 2, BigCVector(zN / 2));
    BigCMatrix ZU1(zN / 2, BigCVector(zN / 2));
    BigCMatrix ZV0(zN / 2, BigCVector(zN / 2));
    BigCMatrix ZV1(zN / 2, BigCVector(zN / 2));
    for (size_t i = 0; i != zN / 2; ++i) {
        for (size_t j = 0; j != zN / 2; ++j) {
            ZU0[i][j] = ZU[i][j];
            ZU1[i][j] = ZU[i][j + zN / 2];
            ZV0[i][j] = ZV[i][j];
            ZV1[i][j] = ZV[i + zN / 2][j];
        }
    }

    // Special matrix for Arith-To-Boolean
    BigCMatrix ZV0SpecialB0 = ZV0;
    for (size_t i = 0; i != zN / 2; ++i) {
        ZV0SpecialB0[0][i] = ZV0SpecialB0[1][i] * BigFixedPoint::two();
    }

    // Special matrices for Arith-To-Arith Noise
    auto tInC    = ZPolynomial::getT(zN).toCSlots();
    auto tInvInC = ZPolynomial::getTInv(zN).toCSlots();
    // It is a diag(tInv encoded in C) * ZU
    BigCMatrix ZU0SpecialA2Ae(zN / 2, BigCVector(zN / 2));
    BigCMatrix ZU1SpecialA2Ae(zN / 2, BigCVector(zN / 2));
    // It is a ZV * diag(t encoded in C)
    BigCMatrix ZV0SpecialA2Ae(zN / 2, BigCVector(zN / 2));
    BigCMatrix ZV1SpecialA2Ae(zN / 2, BigCVector(zN / 2));
    for (size_t i = 0; i != zN / 2; ++i) {
        for (size_t j = 0; j != zN / 2; ++j) {
            ZU0SpecialA2Ae[i][j] = (tInvInC[i] * ZU0[i][j]).scaleTo(128);
            ZU1SpecialA2Ae[i][j] = (tInvInC[i] * ZU1[i][j]).scaleTo(128);
            ZV0SpecialA2Ae[i][j] = (ZV0[i][j] * tInC[j]).scaleTo(128);
            ZV1SpecialA2Ae[i][j] = (ZV1[i][j] * tInC[j]).scaleTo(128);
        }
    }

    BigFixedPoint scaleZU = BigFixedPoint::one();
    BigFixedPoint scaleZV = BigFixedPoint::one();

    if (isSparse) {
        precom->m_ZUPre          = EvalZLinearTransformPrecompute(cc, ZU0, ZU1, zSlots, scaleZU);
        precom->m_ZVPre          = EvalZLinearTransformPrecompute(cc, ZV0, ZV1, zSlots, scaleZV);
        precom->m_ZVSpecialB0Pre = EvalZLinearTransformPrecompute(cc, ZV0SpecialB0, ZV1, zSlots, scaleZV);
        precom->m_ZUSpecialA2AePre =
            EvalZLinearTransformPrecompute(cc, ZU0SpecialA2Ae, ZU1SpecialA2Ae, zSlots, scaleZU);
        precom->m_ZVSpecialA2AePre =
            EvalZLinearTransformPrecompute(cc, ZV0SpecialA2Ae, ZV1SpecialA2Ae, zSlots, scaleZV);
    }
    else {
        precom->m_ZU0Pre            = EvalZLinearTransformPrecompute(cc, ZU0, zSlots, scaleZU);
        precom->m_ZU1Pre            = EvalZLinearTransformPrecompute(cc, ZU1, zSlots, scaleZU);
        precom->m_ZV0Pre            = EvalZLinearTransformPrecompute(cc, ZV0, zSlots, scaleZV);
        precom->m_ZV1Pre            = EvalZLinearTransformPrecompute(cc, ZV1, zSlots, scaleZV);
        precom->m_ZV0SpecialB0Pre   = EvalZLinearTransformPrecompute(cc, ZV0SpecialB0, zSlots, scaleZV);
        precom->m_ZU0SpecialA2AePre = EvalZLinearTransformPrecompute(cc, ZU0SpecialA2Ae, zSlots, scaleZU);
        precom->m_ZU1SpecialA2AePre = EvalZLinearTransformPrecompute(cc, ZU1SpecialA2Ae, zSlots, scaleZU);
        precom->m_ZV0SpecialA2AePre = EvalZLinearTransformPrecompute(cc, ZV0SpecialA2Ae, zSlots, scaleZV);
        precom->m_ZV1SpecialA2AePre = EvalZLinearTransformPrecompute(cc, ZV1SpecialA2Ae, zSlots, scaleZV);
    }

    // LUTs for ArithToBoolean
    uint64_t p            = 1l << w;
    precom->m_w           = w;
    precom->m_cutoff      = arithToBooleanCutoff;
    precom->m_lutIDOrder  = lutOrder;
    precom->m_lutMSBOrder = lutOrder;

    auto lutMSBCoeffs = GetHermiteTrigCoefficients(
        [&](int64_t x) -> int64_t {
            // Input x is in {0, 1, ..., p-1}
            // Extract the negated MSB
            // Because we work in (-1, 0] after Flatten
            if (x <= static_cast<int64_t>(p) / 2 && x != 0) {
                return 1;
            }
            else {
                return 0;
            }
        },
        p, lutOrder, 1);  // We do not rescale here
    auto lutIDCoeffs = GetHermiteTrigCoefficients(
        [&](int64_t x) -> int64_t {
            // Input x is in {0, 1, ..., p-1}
            // Because we work in (-1, 0] after Flatten
            if (x == 0)
                return 0;
            else
                return x - p;
        },
        p, lutOrder, p);

    // Temporary: directly convert to BigComplex
    std::vector<BigComplex> lutMSBCoeffsBC(lutMSBCoeffs.size());
    for (size_t i = 0; i != lutMSBCoeffsBC.size(); ++i) {
        lutMSBCoeffsBC[i] = BigComplex(BigFixedPoint::fromDouble(lutMSBCoeffs[i].real()),
                                       BigFixedPoint::fromDouble(lutMSBCoeffs[i].imag()));
    }
    std::vector<BigComplex> lutIDCoeffsBC(lutIDCoeffs.size());
    for (size_t i = 0; i != lutIDCoeffsBC.size(); ++i) {
        lutIDCoeffsBC[i] = BigComplex(BigFixedPoint::fromDouble(lutIDCoeffs[i].real()),
                                      BigFixedPoint::fromDouble(lutIDCoeffs[i].imag()));
    }
    precom->m_lutMSBCoeffs = lutMSBCoeffsBC;
    precom->m_lutIDCoeffs  = lutIDCoeffsBC;
}

//------------------------------------------------------------------------------
// Precomputations for ZCoeffsToSlots and SlotsToZCoeffs
//------------------------------------------------------------------------------

std::vector<ZBootstrapPlaintextCache> FHEZImpl::EvalZLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                               const BigCMatrix& A,
                                                                               uint32_t zSlotsUnsigned,
                                                                               BigFixedPoint scale) const {
    int32_t zSlots = zSlotsUnsigned;
    int32_t zNDiv2 = A.size();
    if (zNDiv2 != static_cast<int32_t>(A[0].size()))
        OPENFHE_THROW("The matrix passed to EvalLTPrecompute is not square");

    const int32_t cSlots = zSlots * zNDiv2;

    const int32_t step = std::ceil(std::sqrt(zNDiv2));

    // digonal, 1st digonal, ..., zN/2-1th diagonal, -1th diagonal, ..., -zN/2+1 th diagonal
    std::vector<ZBootstrapPlaintextCache> result(2 * zNDiv2 - 1);
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(zNDiv2))
#endif
    for (int32_t ji = 0; ji < zNDiv2; ++ji) {
        auto diag = ExtractShiftedDiagonal(A, ji);
        for (int32_t k = zNDiv2 - ji; k < zNDiv2; ++k) {
            diag[k] = BigFixedPoint::zero();
        }
        auto repeatedDiag = BigCVector(cSlots, BigFixedPoint::zero());
        for (int32_t r = 0; r < zSlots; ++r) {
            for (int32_t k = 0; k < zNDiv2; ++k) {
                repeatedDiag[r * zNDiv2 + k] = diag[k];
            }
        }
        for (auto& d : repeatedDiag)
            d *= scale;
        result[ji] = std::make_shared<ZBootstrapPlaintextCacheImpl>(Rotate(repeatedDiag, -step * (ji / step)));
    }
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(zNDiv2))
#endif
    for (int32_t ji = -1; ji > -zNDiv2; --ji) {
        auto diag = ExtractShiftedDiagonal(A, zNDiv2 + ji);
        for (int32_t k = 0; k < -ji; ++k) {
            diag[k] = BigFixedPoint::zero();
        }
        auto repeatedDiag = BigCVector(cSlots, BigFixedPoint::zero());
        for (int32_t r = 0; r < zSlots; ++r) {
            for (int32_t k = 0; k < zNDiv2; ++k) {
                repeatedDiag[r * zNDiv2 + k] = diag[k];
            }
        }
        for (auto& d : repeatedDiag)
            d *= scale;
        result[zNDiv2 - 1 - ji] =
            std::make_shared<ZBootstrapPlaintextCacheImpl>(Rotate(repeatedDiag, -step * (ji / step)));
    }
    return result;
}

std::vector<ZBootstrapPlaintextCache> FHEZImpl::EvalZLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                               const BigCMatrix& A, const BigCMatrix& B,
                                                                               uint32_t zSlotsUnsigned,
                                                                               BigFixedPoint scale) const {
    int32_t zSlots = zSlotsUnsigned;
    int32_t zNDiv2 = A.size();
    if (zNDiv2 != static_cast<int32_t>(A[0].size()))
        OPENFHE_THROW("The matrix passed to EvalLTPrecompute is not square");

    const int32_t cSlots = zSlots * zNDiv2;

    const int32_t step = std::ceil(std::sqrt(zNDiv2));

    // digonal, 1st digonal, ..., zN/2-1th diagonal, -1th diagonal, ..., -zN/2+1 th diagonal
    std::vector<ZBootstrapPlaintextCache> result(2 * zNDiv2 - 1);

#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(zNDiv2))
#endif
    for (int32_t ji = 0; ji < zNDiv2; ++ji) {
        auto vecA = ExtractShiftedDiagonal(A, ji);
        for (int32_t k = zNDiv2 - ji; k < zNDiv2; ++k) {
            vecA[k] = BigFixedPoint::zero();
        }
        auto vecB = ExtractShiftedDiagonal(B, ji);
        for (int32_t k = zNDiv2 - ji; k < zNDiv2; ++k) {
            vecB[k] = BigFixedPoint::zero();
        }
        //vecA.insert(vecA.end(), vecB.begin(), vecB.end());
        auto repeatedDiag = BigCVector(2 * cSlots, BigFixedPoint::zero());
        for (int32_t r = 0; r < zSlots; ++r) {
            for (int32_t k = 0; k < zNDiv2; ++k) {
                repeatedDiag[r * zNDiv2 + k]          = vecA[k];
                repeatedDiag[r * zNDiv2 + k + cSlots] = vecB[k];
            }
        }
        for (auto& d : repeatedDiag)
            d *= scale;
        result[ji] = std::make_shared<ZBootstrapPlaintextCacheImpl>(Rotate(repeatedDiag, -step * (ji / step)));
    }
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(zNDiv2))
#endif
    for (int32_t ji = -1; ji > -zNDiv2; --ji) {
        auto vecA = ExtractShiftedDiagonal(A, zNDiv2 + ji);
        for (int32_t k = 0; k < -ji; ++k) {
            vecA[k] = BigFixedPoint::zero();
        }
        auto vecB = ExtractShiftedDiagonal(B, zNDiv2 + ji);
        for (int32_t k = 0; k < -ji; ++k) {
            vecB[k] = BigFixedPoint::zero();
        }
        //vecA.insert(vecA.end(), vecB.begin(), vecB.end());
        auto repeatedDiag = BigCVector(2 * cSlots, BigFixedPoint::zero());
        for (int32_t r = 0; r < zSlots; ++r) {
            for (int32_t k = 0; k < zNDiv2; ++k) {
                repeatedDiag[r * zNDiv2 + k]          = vecA[k];
                repeatedDiag[r * zNDiv2 + k + cSlots] = vecB[k];
            }
        }
        for (auto& d : repeatedDiag)
            d *= scale;
        result[zNDiv2 - 1 - ji] =
            std::make_shared<ZBootstrapPlaintextCacheImpl>(Rotate(repeatedDiag, -step * (ji / step)));
    }
    return result;
}

//------------------------------------------------------------------------------
// Precomputations for CoeffsToSlots and SlotsToCoeffs
//------------------------------------------------------------------------------

std::vector<ZBootstrapPlaintextCache> FHEZImpl::EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                              const BigCMatrix& A,
                                                                              BigFixedPoint scale) const {
    const int32_t slots = A.size();
    if (slots != static_cast<int32_t>(A[0].size()))
        OPENFHE_THROW("The matrix passed to EvalLTPrecompute is not square");

    auto g = GetBootPrecom(slots).m_paramsEnc.g;

    const int32_t step = (g == 0) ? std::ceil(std::sqrt(slots)) : g;

    std::vector<ZBootstrapPlaintextCache> result(slots);
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(slots))
#endif
    for (int32_t ji = 0; ji < slots; ++ji) {
        auto diag = ExtractShiftedDiagonal(A, ji);
        for (auto& d : diag)
            d *= scale;
        result[ji] = std::make_shared<ZBootstrapPlaintextCacheImpl>(Rotate(diag, -step * (ji / step)));
    }
    return result;
}

std::vector<ZBootstrapPlaintextCache> FHEZImpl::EvalLinearTransformPrecompute(const CryptoContextImpl<DCRTPoly>& cc,
                                                                              const BigCMatrix& A, const BigCMatrix& B,
                                                                              uint32_t orientation,
                                                                              BigFixedPoint scale) const {
    const int32_t slots = static_cast<int32_t>(A.size());

    auto g = GetBootPrecom(slots).m_paramsEnc.g;

    const int32_t step = (g == 0) ? std::ceil(std::sqrt(slots)) : g;

    std::vector<ZBootstrapPlaintextCache> result(slots);

    if (orientation == 0) {
        // vertical concatenation - used during homomorphic encoding
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(slots))
#endif
        for (int32_t ji = 0; ji < slots; ++ji) {
            auto vecA = ExtractShiftedDiagonal(A, ji);
            auto vecB = ExtractShiftedDiagonal(B, ji);
            vecA.insert(vecA.end(), vecB.begin(), vecB.end());
            for (auto& d : vecA)
                d *= scale;
            result[ji] = std::make_shared<ZBootstrapPlaintextCacheImpl>(Rotate(vecA, -step * (ji / step)));
        }
    }
    else {
        // horizontal concatenation - used during homomorphic decoding
        BigCMatrix newA(slots);

        //  A and B are concatenated horizontally
        for (int32_t i = 0; i < slots; ++i) {
            newA[i].reserve(A[i].size() + B[i].size());
            newA[i].insert(newA[i].end(), A[i].begin(), A[i].end());
            newA[i].insert(newA[i].end(), B[i].begin(), B[i].end());
        }

#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(slots))
#endif
        for (int32_t ji = 0; ji < slots; ++ji) {
            // shifted diagonal is computed for rectangular map newA of dimension
            // slots x 2*slots
            auto vec = ExtractShiftedDiagonal(newA, ji);
            auto res = Rotate(vec, -step * (ji / step));
            for (auto& d : vec)
                d *= scale;
            result[ji] = std::make_shared<ZBootstrapPlaintextCacheImpl>(res);
        }
    }

    return result;
}

std::vector<std::vector<ZBootstrapPlaintextCache>> FHEZImpl::EvalCoeffsToSlotsPrecompute(
    const CryptoContextImpl<DCRTPoly>& cc, const BigCVector& A, const std::vector<uint32_t>& rotGroup, bool flag_i,
    std::vector<BigFixedPoint> scales) const {
    const uint32_t slots = rotGroup.size();

    const auto& p = GetBootPrecom(slots).m_paramsEnc;

    // result is the rotated plaintext version of the coefficients
    std::vector<std::vector<ZBootstrapPlaintextCache>> result(p.lvlb,
                                                              std::vector<ZBootstrapPlaintextCache>(p.numRotations));

    int32_t stop    = -1;
    int32_t flagRem = 0;
    if (p.remCollapse != 0) {
        stop    = 0;
        flagRem = 1;

        // remainder corresponds to index 0 in encoding and to last index in decoding
        result[0].resize(p.numRotationsRem);
    }

    auto M = cc.GetCyclotomicOrder();

    if (uint32_t M4 = M / 4; slots == M4) {
        //------------------------------------------------------------------------------
        // fully-packed mode
        //------------------------------------------------------------------------------

        auto coeff = CoeffEncodingCollapse(A, rotGroup, p.lvlb, flag_i);

        for (int32_t s = -1 + p.lvlb; s > stop; --s) {
            const int32_t rotScale = (1 << ((s - flagRem) * p.layersCollapse + p.remCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    // Maybe we need it one day???
                    //if ((flagRem == 0) && (s == stop + 1)) {
                    //    // do the scaling only at the last set of coefficients
                    //}
                    // We do the scaling at every level
                    for (auto& c : coeff[s][ij])
                        c *= scales[s];

                    auto rot = Rotate(coeff[s][ij], ReduceRotation(-rotScale * (ij / p.g), slots));

                    result[s][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }

        if (flagRem == 1) {
            const uint32_t limit = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    for (auto& c : coeff[stop][ij])
                        c *= scales[stop];

                    auto rot = Rotate(coeff[stop][ij], ReduceRotation(-p.gRem * (ij / p.gRem), slots));

                    result[stop][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }
    }
    else {
        //------------------------------------------------------------------------------
        // sparsely-packed mode
        //------------------------------------------------------------------------------

        auto coeff  = CoeffEncodingCollapse(A, rotGroup, p.lvlb, false);
        auto coeffi = CoeffEncodingCollapse(A, rotGroup, p.lvlb, true);

        for (int32_t s = -1 + p.lvlb; s > stop; --s) {
            const int32_t rotScale = (1 << ((s - flagRem) * p.layersCollapse + p.remCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    // concatenate the coefficients horizontally on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[s][ij];
                    auto& clearTmpi = coeffi[s][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    // Maybe we need it one day
                    //if ((flagRem == 0) && (s == stop + 1)) {
                    //    // do the scaling only at the last set of coefficients
                    //    //for (auto& c : clearTmp)
                    //    //    c *= scale;
                    //}
                    // We do the scaling at every level
                    for (auto& c : clearTmp)
                        c *= scales[s];

                    auto rot = Rotate(clearTmp, ReduceRotation(-rotScale * (ij / p.g), M4));

                    result[s][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }

        if (flagRem == 1) {
            const uint32_t limit = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    // concatenate the coefficients on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[stop][ij];
                    auto& clearTmpi = coeffi[stop][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    for (auto& c : clearTmp)
                        c *= scales[stop];

                    auto rot = Rotate(clearTmp, ReduceRotation(-p.gRem * (ij / p.gRem), M4));

                    result[stop][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }
    }
    return result;
}

std::vector<std::vector<ZBootstrapPlaintextCache>> FHEZImpl::EvalSlotsToCoeffsPrecompute(
    const CryptoContextImpl<DCRTPoly>& cc, const BigCVector& A, const std::vector<uint32_t>& rotGroup,
    bool flag_i) const {
    const uint32_t slots = rotGroup.size();

    const auto& p = GetBootPrecom(slots).m_paramsDec;

    const int32_t flagRem = (p.remCollapse == 0) ? 0 : 1;

    // result is the rotated plaintext version of coeff
    std::vector<std::vector<ZBootstrapPlaintextCache>> result(p.lvlb,
                                                              std::vector<ZBootstrapPlaintextCache>(p.numRotations));
    if (flagRem == 1) {
        // remainder corresponds to index 0 in encoding and to last index in decoding
        result[p.lvlb - 1].resize(p.numRotationsRem);
    }

    // make sure the plaintext is created only with the necessary amount of moduli

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc.GetCryptoParameters());
    auto elementParams      = *(cryptoParams->GetElementParams());

    if (uint32_t M4 = cc.GetCyclotomicOrder() / 4; M4 == slots) {
        // fully-packed
        auto coeff          = CoeffDecodingCollapse(A, rotGroup, p.lvlb, flag_i);
        const uint32_t smax = p.lvlb - flagRem;
        for (uint32_t s = 0; s < smax; ++s) {
            const int32_t rotScale = (1 << (s * p.layersCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    if ((flagRem == 0) && (s + 1 == smax)) {
                        // do the scaling only at the last set of coefficients
                        //for (auto& c : coeff[s][ij])
                        //    c *= scale;
                    }

                    auto rot = Rotate(coeff[s][ij], ReduceRotation(-rotScale * (ij / p.g), slots));

                    result[s][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }

        if (flagRem == 1) {
            const int32_t rotScale = (1 << (smax * p.layersCollapse)) * p.gRem;
            const uint32_t limit   = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    //for (auto& c : coeff[smax][ij])
                    //    c *= scale;

                    auto rot = Rotate(coeff[smax][ij], ReduceRotation(-rotScale * (ij / p.gRem), slots));

                    result[smax][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }
    }
    else {
        //------------------------------------------------------------------------------
        // sparsely-packed mode
        //------------------------------------------------------------------------------

        auto coeff  = CoeffDecodingCollapse(A, rotGroup, p.lvlb, false);
        auto coeffi = CoeffDecodingCollapse(A, rotGroup, p.lvlb, true);

        const uint32_t smax = p.lvlb - flagRem;
        for (uint32_t s = 0; s < smax; ++s) {
            const int32_t rotScale = (1 << (s * p.layersCollapse)) * p.g;
            const uint32_t limit   = p.b * p.g;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotations) {
                    // concatenate the coefficients horizontally on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[s][ij];
                    auto& clearTmpi = coeffi[s][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    if ((flagRem == 0) && (s + 1 == smax)) {
                        // do the scaling only at the last set of coefficients
                        //for (auto& c : clearTmp)
                        //    c *= scale;
                    }

                    auto rot = Rotate(clearTmp, ReduceRotation(-rotScale * (ij / p.g), M4));

                    result[s][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }

        if (flagRem == 1) {
            const int32_t rotScale = (1 << (smax * p.layersCollapse)) * p.g;
            const uint32_t limit   = p.bRem * p.gRem;
#if !defined(__MINGW32__) && !defined(__MINGW64__)
    #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(limit))
#endif
            for (uint32_t ij = 0; ij < limit; ++ij) {
                if (ij != p.numRotationsRem) {
                    // concatenate the coefficients on their third dimension, which corresponds to the # of slots
                    auto clearTmp   = coeff[smax][ij];
                    auto& clearTmpi = coeffi[smax][ij];
                    clearTmp.insert(clearTmp.end(), clearTmpi.begin(), clearTmpi.end());
                    //for (auto& c : clearTmp)
                    //    c *= scale;

                    auto rot = Rotate(clearTmp, ReduceRotation(-rotScale * (ij / p.gRem), M4));

                    result[smax][ij] = std::make_shared<ZBootstrapPlaintextCacheImpl>(rot);
                }
            }
        }
    }
    return result;
}

Plaintext ZBootstrapPlaintextCacheImpl::GetPlaintext(const BigFixedPoint& scalingFactor,
                                                     const std::shared_ptr<typename DCRTPoly::Params>& elementParams) {
    auto q   = elementParams->GetModulus();
    auto key = std::make_tuple(scalingFactor, q);
    auto it  = m_cache.find(key);
    if (it != m_cache.end()) {
        return it->second;
    }
    // Construct plaintext
    auto n = m_value.size() * 2;

    ZEncodingParams params(CMode, n);
    Plaintext ptxt = ZEncodingImpl::encodeC(CSlots(params, m_value), elementParams, scalingFactor);
    // Store in cache
    // FIXME: is this thread safe???
    m_cache[key] = ptxt;
    return ptxt;
}

Plaintext FHEZImpl::getATBMaskSparsePacking(uint32_t iter, uint32_t w, uint32_t zN, uint32_t zSlots,
                                            const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                            const BigFixedPoint& scalingFactor) {
    auto cSlots = zN * zSlots / 2;
    // Encode low-hot vector
    std::vector<BigComplex> oneHotVec(2 * cSlots, BigFixedPoint::zero());
    for (uint32_t j = 0; j != zSlots; ++j) {
        auto index = j * (zN / 2) + (iter * w);
        if (iter * w >= zN / 2) {
            index += (zSlots - 1) * (zN / 2);
        }
        for (uint32_t b = 0; b != w; ++b) {
            oneHotVec[index + b] = BigFixedPoint::one();
        }
    }

    ZEncodingParams oneHotZEncodeParams(CMode, zN * zSlots * 2);  // sparse packing
    RPolynomial oneHotPoly = CSlots(oneHotZEncodeParams, oneHotVec).toRPolynomial();

    Plaintext oneHotPtxt = ZEncodingImpl::encodeR(oneHotPoly, elementParams, scalingFactor);
    return oneHotPtxt;
}

Plaintext FHEZImpl::getATBMaskFullPacking(uint32_t iter, uint32_t w, uint32_t zN, uint32_t zSlots,
                                          const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                          const BigFixedPoint& scalingFactor) {
    auto cSlots = zN * zSlots / 2;
    // Encode low-hot vector
    std::vector<BigComplex> oneHotVec(cSlots, BigFixedPoint::zero());
    for (uint32_t j = 0; j != zSlots; ++j) {
        auto index = j * (zN / 2);
        if (iter * w >= zN / 2) {
            index += (iter * w) - zN / 2;
        }
        else {
            index += (iter * w);
        }
        for (uint32_t b = 0; b != w; ++b) {
            oneHotVec[index + b] = BigFixedPoint::one();
        }
    }

    ZEncodingParams oneHotZEncodeParams(CMode, zN * zSlots);  // full packing
    RPolynomial oneHotPoly = CSlots(oneHotZEncodeParams, oneHotVec).toRPolynomial();

    Plaintext oneHotPtxt = ZEncodingImpl::encodeR(oneHotPoly, elementParams, scalingFactor);
    return oneHotPtxt;
}

Plaintext FHEZImpl::getATBSubtractMaskFullPacking(uint32_t iter, uint32_t w, uint32_t zN, uint32_t zSlots,
                                                  const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                                  const BigFixedPoint& scalingFactor, BigFixedPoint scaleDown) {
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));
    // wrap around if iter exceeds the number of iterations needed to cover all slots
    // This is for B-ATB where we need to reuse the same masks for multiple ciphertexts
    iter = iter % numIter;

    auto key = std::make_tuple(iter, w, zN, zSlots, elementParams->GetModulus(), scalingFactor, scaleDown);
    auto it  = m_atbSubtractMaskPtxtCache.find(key);
    if (it != m_atbSubtractMaskPtxtCache.end()) {
        return it->second;
    }

    auto cSlots = zN * zSlots / 2;
    // Encode low-hot vector
    std::vector<BigComplex> oneHotVec(cSlots, BigFixedPoint::zero());
    for (uint32_t j = 0; j != zSlots; ++j) {
        auto index = j * (zN / 2);
        if (iter * w >= zN / 2) {
            index += (iter * w) - zN / 2;
        }
        else {
            index += (iter * w);
        }
        for (uint32_t b = 0; b != w; ++b) {
            oneHotVec[index + b] = BigFixedPoint::one() * scaleDown;
        }
    }

    ZEncodingParams oneHotZEncodeParams(CMode, zN * zSlots);  // full packing
    RPolynomial oneHotPoly = CSlots(oneHotZEncodeParams, oneHotVec).toRPolynomial();

    Plaintext oneHotPtxt = ZEncodingImpl::encodeR(oneHotPoly, elementParams, scalingFactor);

    m_atbSubtractMaskPtxtCache[key] = oneHotPtxt;
    return oneHotPtxt;
}

Plaintext FHEZImpl::getATBMask(uint32_t iter, uint32_t w, uint32_t zN, uint32_t zSlots,
                               const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                               const BigFixedPoint& scalingFactor) {
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));
    // wrap around if iter exceeds the number of iterations needed to cover all slots
    // This is for B-ATB where we need to reuse the same masks for multiple ciphertexts
    iter = iter % numIter;

    auto ringDim         = elementParams->GetRingDimension();
    bool isSparse        = (zN * zSlots != ringDim);
    MaskPlaintextKey key = std::make_tuple(iter, w, zN, zSlots, elementParams->GetModulus(), scalingFactor);
    auto it              = m_atbMaskPtxtCache.find(key);
    if (it != m_atbMaskPtxtCache.end()) {
        return it->second;
    }
    Plaintext ptxt;
    if (isSparse) {
        ptxt = getATBMaskSparsePacking(iter, w, zN, zSlots, elementParams, scalingFactor);
    }
    else {
        ptxt = getATBMaskFullPacking(iter, w, zN, zSlots, elementParams, scalingFactor);
    }
    m_atbMaskPtxtCache[key] = ptxt;
    return ptxt;
}

}  // namespace lbcrypto