#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/ckksrns-fhe.h"

namespace lbcrypto {
//------------------------------------------------------------------------------
// Bootstrap KeyGen
//------------------------------------------------------------------------------

namespace {
EvalKey<DCRTPoly> ConjugateKeyGen(const PrivateKey<DCRTPoly> privateKey) {
    uint32_t N = privateKey->GetPrivateElement().GetRingDimension();
    std::vector<uint32_t> vec(N);
    PrecomputeAutoMap(N, 2 * N - 1, &vec);
    const auto cc   = privateKey->GetCryptoContext();
    auto pkPermuted = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
    pkPermuted->SetPrivateElement(privateKey->GetPrivateElement().AutomorphismTransform(2 * N - 1, vec));
    pkPermuted->SetKeyTag(privateKey->GetKeyTag());
    return cc->GetScheme()->KeySwitchGen(privateKey, pkPermuted);
}
}  // namespace

void FHEZImpl::EvalBootstrapKeyGen(const PrivateKey<DCRTPoly> privateKey, uint32_t zN, uint32_t zSlots) {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(privateKey->GetCryptoParameters());

    if (cryptoParams->GetKeySwitchTechnique() != HYBRID)
        OPENFHE_THROW("CKKS Bootstrapping is only supported for the Hybrid key switching method.");

    auto cc   = privateKey->GetCryptoContext();
    auto algo = cc->GetScheme();
    auto M    = cc->GetCyclotomicOrder();

    // computing all indices for baby-step giant-step procedure
    auto evalKeys = algo->EvalAtIndexKeyGen(privateKey, FindBootstrapRotationIndices(zN, zSlots, M));

    (*evalKeys)[M - 1] = ConjugateKeyGen(privateKey);

    if (cryptoParams->GetSecretKeyDist() == SPARSE_ENCAPSULATED) {
        DCRTPoly::TugType tug;

        // sparse key used for the modraising step
        auto skNew = std::make_shared<PrivateKeyImpl<DCRTPoly>>(cc);
        skNew->SetPrivateElement(DCRTPoly(tug, cryptoParams->GetElementParams(), Format::EVALUATION, 32));

        // we reserve M-4 and M-2 for the sparse encapsulation switching keys
        // Even autorphism indices are not possible, so there will not be any conflict
        (*evalKeys)[M - 4] = FHECKKSRNS::KeySwitchGenSparse(privateKey, skNew);
        (*evalKeys)[M - 2] = algo->KeySwitchGen(skNew, privateKey);
    }

    cc->InsertEvalAutomorphismKey(evalKeys, privateKey->GetKeyTag());
}

//------------------------------------------------------------------------------
// Find Rotation Indices
//------------------------------------------------------------------------------

std::vector<int32_t> FHEZImpl::FindBootstrapRotationIndices(uint32_t zN, uint32_t zSlots, uint32_t M) {
    // This is cSlots
    auto slots    = zN * zSlots / 2;
    const auto& p = GetBootPrecom(slots);

    // Remove possible duplicates and remove automorphisms corresponding to 0 and M/4 by using std::set
    std::set<uint32_t> s;
    if (p.m_paramsEnc.lvlb == 1 && p.m_paramsDec.lvlb == 1) {
        auto tmp = FindLinearTransformRotationIndices(slots, M);
        s.insert(tmp.begin(), tmp.end());
    }
    else {
        auto tmp = FindCoeffsToSlotsRotationIndices(slots, M);
        s.insert(tmp.begin(), tmp.end());
        tmp = FindSlotsToCoeffsRotationIndices(slots, M);
        s.insert(tmp.begin(), tmp.end());
    }

    // Now this is Z specific
    // First, all rotations for Z LT
    {
        // Computing the baby-step g and the giant-step h.
        auto zNDiv2 = zN / 2;
        // FUNNY that bStep = g...
        const int32_t g = std::ceil(std::sqrt(zNDiv2));
        const int32_t h = std::ceil(static_cast<double>(zNDiv2) / g);

        // We have both positive and negative rotations
        for (int32_t i = 1; i <= g; ++i) {
            s.insert(i);
            s.insert(-i);
        }
        for (int32_t i = 2; i < h; ++i) {
            s.insert(i * g);
            s.insert(-i * g);
        }
    }
    // Then all rotations for ArithToBoolean remove low
    {
        auto w       = p.m_w;
        auto numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));
        for (uint32_t iter = 0; iter != numIter; ++iter) {
            for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
                int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
                // Note the rotation index is negative here
                int32_t rotationIndex = diff * w;
                if (rotationIndex <= p.m_cutoff) {
                    // We do not remove them any more
                    // Just treat the lower parts as noises
                    break;
                }
                // If cross the half-way point, need to rotate more
                if (nextIter * w >= zN / 2 && iter * w < zN / 2) {
                    rotationIndex -= static_cast<int32_t>((zSlots - 1) * zN / 2);
                }
                s.insert(rotationIndex);
            }
        }
    }
    // For all rotations for ArithToBoolean batched
    {
        auto w         = p.m_w;
        auto numIter   = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));
        auto batchSize = numIter / 2;

        // For pre-rotation and LUT.MSB result combination
        for (int32_t i = 0; i != static_cast<int32_t>(batchSize); ++i) {
            s.insert(i * w);
            s.insert(-i * w);
        }

        for (uint32_t iter = 0; iter != numIter; ++iter) {
            for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
                int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
                // Note the rotation index is negative here
                int32_t rotationIndex = diff * w;
                if (rotationIndex <= p.m_cutoff) {
                    // We do not remove them any more
                    // Just treat the lower parts as noises
                    break;
                }
                if (nextIter * w >= zN / 2 && iter * w < zN / 2) {
                    rotationIndex += zN / 2;
                }
                s.insert(rotationIndex);
            }
        }
    }

    s.erase(0);
    s.erase(M / 4);
    return std::vector<int32_t>(s.begin(), s.end());
}

// ATTN: This function is a helper methods to be called in FindBootstrapRotationIndices() only.
// so it DOES NOT remove possible duplicates and automorphisms corresponding to 0 and M/4.
// This method completely depends on FindBootstrapRotationIndices() to do that.
std::vector<uint32_t> FHEZImpl::FindLinearTransformRotationIndices(uint32_t slots, uint32_t M) {
    // Computing the baby-step g and the giant-step h.
    const auto& p    = GetBootPrecom(slots);
    const uint32_t g = (p.m_paramsEnc.g == 0) ? std::ceil(std::sqrt(slots)) : p.m_paramsEnc.g;
    const uint32_t h = std::ceil(static_cast<double>(slots) / g);

    // To avoid overflowing uint32_t variables, we do some math operations below in a specific order
    // computing all indices for baby-step giant-step procedure
    const int32_t indexListSz = static_cast<int32_t>(g) + h + M - 2;
    if (indexListSz < 0)
        OPENFHE_THROW("indexListSz can not be negative");

    std::vector<uint32_t> indexList;
    indexList.reserve(indexListSz);

    for (uint32_t i = 1; i <= g; ++i)
        indexList.emplace_back(i);
    for (uint32_t i = 2; i < h; ++i)
        indexList.emplace_back(i * g);

    // additional automorphisms are needed for sparse bootstrapping
    if (uint32_t m = slots * 4; m != M) {
        for (uint32_t j = 1; j < M / m; j <<= 1)
            indexList.emplace_back(j * slots);
    }

    return indexList;
}

// ATTN: This function is a helper methods to be called in FindBootstrapRotationIndices() only.
// so it DOES NOT remove possible duplicates and automorphisms corresponding to 0 and M/4.
// This method completely depends on FindBootstrapRotationIndices() to do that.
std::vector<uint32_t> FHEZImpl::FindCoeffsToSlotsRotationIndices(uint32_t slots, uint32_t M) {
    const auto& p = GetBootPrecom(slots).m_paramsEnc;

    // To avoid overflowing uint32_t variables, we do some math operations below in a specific order
    // Computing all indices for baby-step giant-step procedure for encoding and decoding
    const int32_t indexListSz = static_cast<int32_t>(p.b) + p.g - 2 + p.bRem + p.gRem - 2 + 1 + M;
    if (indexListSz < 0)
        OPENFHE_THROW("indexListSz can not be negative");

    std::vector<uint32_t> indexList;
    indexList.reserve(indexListSz);

    // additional automorphisms are needed for sparse bootstrapping
    if (uint32_t m = slots * 4; m != M) {
        for (uint32_t j = 1; j < M / m; j <<= 1)
            indexList.emplace_back(j * slots);
    }

    M >>= 2;
    const int32_t flagRem   = (p.remCollapse == 0) ? 0 : 1;
    const int32_t halfRots  = 1 - (p.numRotations + 1) / 2;
    const int32_t halfRotsg = halfRots + p.g;
    for (int32_t s = -1 + p.lvlb; s >= flagRem; --s) {
        const uint32_t scalingFactor = 1U << ((s - flagRem) * p.layersCollapse + p.remCollapse);
        for (int32_t j = halfRots; j < halfRotsg; ++j)
            indexList.emplace_back(ReduceRotation(j * scalingFactor, slots));
        for (uint32_t i = 0; i < p.b; ++i)
            indexList.emplace_back(ReduceRotation(i * p.g * scalingFactor, M));
    }

    if (flagRem == 1) {
        const int32_t halfRotsRem  = (1 - (p.numRotationsRem + 1) / 2);
        const int32_t halfRotsRemg = halfRotsRem + p.gRem;
        for (int32_t j = halfRotsRem; j < halfRotsRemg; ++j)
            indexList.emplace_back(ReduceRotation(j, slots));
        for (uint32_t i = 0; i < p.bRem; ++i)
            indexList.emplace_back(ReduceRotation(i * p.gRem, M));
    }

    return indexList;
}

std::vector<uint32_t> FHEZImpl::FindSlotsToCoeffsRotationIndices(uint32_t slots, uint32_t M) {
    const auto& p = GetBootPrecom(slots).m_paramsDec;

    // To avoid overflowing uint32_t variables, we do some math operations below in a specific order
    // Computing all indices for baby-step giant-step procedure for encoding and decoding
    const int32_t indexListSz = static_cast<int32_t>(p.b) + p.g - 2 + p.bRem + p.gRem - 2 + 1 + M;
    if (indexListSz < 0)
        OPENFHE_THROW("indexListSz can not be negative");

    std::vector<uint32_t> indexList;
    indexList.reserve(indexListSz);

    // additional automorphisms are needed for sparse bootstrapping
    if (uint32_t m = slots * 4; m != M) {
        for (uint32_t j = 1; j < M / m; j <<= 1)
            indexList.emplace_back(j * slots);
    }

    M >>= 2;
    const uint32_t flagRem  = (p.remCollapse == 0) ? 0 : 1;
    const uint32_t smax     = p.lvlb - flagRem;
    const int32_t halfRots  = (1 - (p.numRotations + 1) / 2);
    const int32_t halfRotsg = halfRots + p.g;
    for (uint32_t s = 0; s < smax; ++s) {
        const uint32_t scalingFactor = 1U << (s * p.layersCollapse);
        for (int32_t j = halfRots; j < halfRotsg; ++j)
            indexList.emplace_back(ReduceRotation(j * scalingFactor, M));
        for (uint32_t i = 0; i < p.b; ++i)
            indexList.emplace_back(ReduceRotation(i * p.g * scalingFactor, M));
    }

    if (flagRem == 1) {
        const uint32_t scalingFactor = 1U << (smax * p.layersCollapse);
        const int32_t halfRotsRem    = (1 - (p.numRotationsRem + 1) / 2);
        const int32_t halfRotsRemg   = halfRotsRem + p.gRem;
        for (int32_t j = halfRotsRem; j < halfRotsRemg; ++j)
            indexList.emplace_back(ReduceRotation(j * scalingFactor, M));
        for (uint32_t i = 0; i < p.bRem; ++i)
            indexList.emplace_back(ReduceRotation(i * p.gRem * scalingFactor, M));
    }

    return indexList;
}

}  // namespace lbcrypto