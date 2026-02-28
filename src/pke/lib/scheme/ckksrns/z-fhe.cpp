#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/z-pke.h"
#include "scheme/ckksrns/ckksrns-fhe.h"
#include "encoding/z-encoding.h"

namespace lbcrypto {

// Not used
Ciphertext<DCRTPoly> FHEZImpl::EvalTruncate(ConstCiphertext<DCRTPoly>& ct) const {
    auto q     = ct->GetElements()[0].GetModulus();
    auto sfNow = ct->GetScalingFactorBFP();
    auto qBFP  = BigFixedPoint(q, 0, false).scaleTo(128);
    // q / Delta
    auto divScalar = (qBFP / sfNow).getRoundedBigInteger();
    auto ct2       = z->EvalMultScalar(ct, divScalar);
    ct2->SetScalingFactorBFP(qBFP);
    // Reduce all the way to the bottom
    z->ModReduceInPlace(ct2, ct2->GetElements()[0].GetNumOfElements() - 1);
    // To make sure the scaling factor is exactly q0 / 2
    auto q0     = ct2->GetElements()[0].GetModulus();
    auto q0BFP  = BigFixedPoint(q0, 0, false).scaleTo(128);
    auto sfNow2 = q0BFP;
    ct2->SetScalingFactorBFP(sfNow2);
    return ct2;
}

Ciphertext<DCRTPoly> FHEZImpl::EvalModRaise(ConstCiphertext<DCRTPoly>& ct) const {
    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(ct->GetCryptoParameters());

    auto paramsQ   = cryptoParams->GetElementParams()->GetParams();
    uint32_t sizeQ = paramsQ.size();
    std::vector<NativeInteger> moduli(sizeQ);
    std::vector<NativeInteger> roots(sizeQ);
    for (uint32_t i = 0; i < sizeQ; ++i) {
        moduli[i] = paramsQ[i]->GetModulus();
        roots[i]  = paramsQ[i]->GetRootOfUnity();
    }

    auto cc = ct->GetCryptoContext();
    auto M  = cc->GetCyclotomicOrder();
    auto N  = cc->GetRingDimension();

    auto elementParamsRaisedPtr = std::make_shared<ILDCRTParams<DCRTPoly::Integer>>(M, moduli, roots);

    auto raised = ct->Clone();
    auto algo   = cc->GetScheme();

    uint32_t L0 = cryptoParams->GetElementParams()->GetParams().size();

    if (cryptoParams->GetSecretKeyDist() == SPARSE_ENCAPSULATED) {
        auto evalKeyMap = cc->GetEvalAutomorphismKeyMap(raised->GetKeyTag());

        // transform from a denser secret to a sparser one
        raised = FHECKKSRNS::KeySwitchSparse(raised, evalKeyMap.at(2 * N - 4));

        // Only level 0 ciphertext used here. Other towers ignored to make CKKS bootstrapping faster.
        auto& ctxtDCRTs = raised->GetElements();

        for (auto& dcrt : ctxtDCRTs) {
            dcrt.SetFormat(COEFFICIENT);
            DCRTPoly tmp(dcrt.GetElementAtIndex(0), elementParamsRaisedPtr);
            tmp.SetFormat(EVALUATION);
            dcrt = std::move(tmp);
        }
        raised->SetLevel(L0 - ctxtDCRTs[0].GetNumOfElements());

        // go back to a denser secret
        algo->KeySwitchInPlace(raised, evalKeyMap.at(2 * N - 2));
    }
    else {
        OPENFHE_THROW("Other Key unsupported");
    }

    raised->SetScalingFactorBFP(ct->GetScalingFactorBFP());
    return raised;
}

void FHEZImpl::EvalPartialSumInPlace(Ciphertext<DCRTPoly>& ct) const {
    auto cc     = ct->GetCryptoContext();
    auto cSlots = ct->GetZEncodingParams().getCSlots();
    auto N      = cc->GetRingDimension();

    const uint32_t limit = N / (cSlots * 2);
    for (uint32_t j = 1; j < limit; j <<= 1) {
        cc->EvalAddInPlace(ct, cc->EvalRotate(ct, j * (cSlots)));
    }
    // Now the message is multplied by N/(rN)
}

Ciphertext<DCRTPoly> FHEZImpl::EvalModRaisePartialSum(ConstCiphertext<DCRTPoly>& ct) const {
    // Truncate is fused into C2R
    auto ctNew = EvalModRaise(ct);
    EvalPartialSumInPlace(ctNew);
    return ctNew;
}

CiphertextGroup FHEZImpl::EvalZ2C(ConstCiphertext<DCRTPoly>& ct, Z2COption z2cOption) const {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto isSparse = precomp.m_isSparse;

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    // adjust ciphertext to level just above z2c + (lMask) + c2r
    uint32_t lZ2C    = 1;
    uint32_t lMask   = (z2cOption == Z2C_SPECIAL_B0);  // when specialB0, leave one level for masking.
    uint32_t lC2R    = precomp.m_paramsDec.lvlb;
    auto mulDepth    = cryptoParams->GetMultiplicativeDepth();
    auto levelForZ2C = mulDepth - lZ2C - lMask - lC2R;
    auto ctNew       = z->AdjustCiphertextToLevel(ct, levelForZ2C);

    const CiphertextGroup::MapFunc postProcess = [&, this](ConstCiphertext<DCRTPoly>& z2c) -> Ciphertext<DCRTPoly> {
        // Take the Real
        auto ct = z->EvalAdd(z2c, z->EvalConjugateInC(z2c));
        z->ModReduceInPlace(ct);
        return ct;
    };

    if (isSparse) {
        auto z2c = EvalZLinearTransform(z2cOption == Z2C_SPECIAL_A2AE ?
                                            precomp.m_ZVSpecialA2AePre :
                                            (z2cOption == Z2C_SPECIAL_B0 ? precomp.m_ZVSpecialB0Pre : precomp.m_ZVPre),
                                        ctNew);
        return postProcess(z2c);
    }
    else {
        auto z2c0 =
            EvalZLinearTransform(z2cOption == Z2C_SPECIAL_A2AE ?
                                     precomp.m_ZV0SpecialA2AePre :
                                     (z2cOption == Z2C_SPECIAL_B0 ? precomp.m_ZV0SpecialB0Pre : precomp.m_ZV0Pre),
                                 ctNew);
        auto z2c1 =
            EvalZLinearTransform(z2cOption == Z2C_SPECIAL_A2AE ? precomp.m_ZV1SpecialA2AePre : precomp.m_ZV1Pre, ctNew);

        CiphertextGroup z2cGroup({z2c0, z2c1});
        return z2cGroup.map(postProcess);
    }
}

// The m0 + i * m1 preprocess is done elsewhere
Ciphertext<DCRTPoly> FHEZImpl::EvalC2R(ConstCiphertext<DCRTPoly>& ct) const {
    auto cc            = ct->GetCryptoContext();
    auto cSlots        = ct->GetZEncodingParams().getCSlots();
    auto precomp       = GetBootPrecom(cSlots);
    bool isSparse      = precomp.m_isSparse;
    bool isLTBootstrap = (precomp.m_paramsEnc.lvlb == 1) && (precomp.m_paramsDec.lvlb == 1);

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    // adjust ciphertext to level just above c2r + trunc
    uint32_t lC2R    = precomp.m_paramsDec.lvlb;
    auto mulDepth    = cryptoParams->GetMultiplicativeDepth();
    auto levelForC2R = mulDepth - lC2R;
    auto ctNew       = z->AdjustCiphertextToLevel(ct, levelForC2R);

    if (isSparse) {
        auto c2r =
            isLTBootstrap ? EvalLinearTransform(precomp.m_U0Pre, ctNew) : EvalSlotsToCoeffs(precomp.m_U0PreFFT, ctNew);
        // Trace
        z->EvalAddInPlace(c2r, cc->EvalRotate(c2r, cSlots));
        z->ModReduceInPlace(c2r);
        return c2r;
    }
    else {
        auto c2r = (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0Pre, ctNew) :
                                     EvalSlotsToCoeffs(precomp.m_U0PreFFT, ctNew);
        z->ModReduceInPlace(c2r);
        return c2r;
    }
}

CiphertextGroup FHEZImpl::EvalR2C(ConstCiphertext<DCRTPoly>& ct, R2CScalingOption scalingOption) const {
    auto cc            = ct->GetCryptoContext();
    auto cSlots        = ct->GetZEncodingParams().getCSlots();
    auto precomp       = GetBootPrecom(cSlots);
    bool isLTBootstrap = (precomp.m_paramsEnc.lvlb == 1) && (precomp.m_paramsDec.lvlb == 1);
    auto isSparse      = precomp.m_isSparse;

    auto scaleDown = [&](Ciphertext<DCRTPoly> target, BigFixedPoint scale) {
        // Manually scale down by N
        auto sfNow  = target->GetScalingFactorBFP();
        auto scalar = (sfNow * scale).getRoundedBigInteger();
        //{
        //    std::cout << "R2C scaleDown factor: " << scale.log2Norm() << std::endl;
        //    auto scalarWhole = sfNow * scale;
        //    std::cout << "scaleWhole: " << scalarWhole.toHexString() << std::endl;
        //    auto scalarFrac = scalarWhole - scalarWhole.round();
        //    auto precLoss   = scalarFrac / scalarWhole;
        //    std::cout << "R2C Scale down precision loss: " << precLoss.log2Norm() << std::endl;
        //}
        z->EvalMultScalarInPlace(target, scalar);
        target->SetScalingFactorBFP(sfNow * sfNow);
        z->ModReduceInPlace(target);
    };

    BigFixedPoint scaleN  = BigFixedPoint::one() / BigFixedPoint::positive(cc->GetRingDimension());
    BigFixedPoint scaleNK = BigFixedPoint::one() / (BigFixedPoint::positive(cc->GetRingDimension()) *
                                                    BigFixedPoint::positive(K_SPARSE_ENCAPSULATED));

    Ciphertext<DCRTPoly> r2c;

    // only one linear transform is needed as the other one can be derived
    if (scalingOption == SCALE_N || scalingOption == SCALE_NK) {
        r2c = (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0hatTPre, ct) :
                                EvalCoeffsToSlots(precomp.m_U0hatTPreFFT, ct);
    }
    else if (scalingOption == SCALE_N_PRE) {
        r2c = (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0hatTPreScaledN, ct) :
                                EvalCoeffsToSlots(precomp.m_U0hatTPreFFTScaledN, ct);
    }
    else {  // scalingOption == SCALE_NK_PRE
        r2c = (isLTBootstrap) ? EvalLinearTransform(precomp.m_U0hatTPreScaledNK, ct) :
                                EvalCoeffsToSlots(precomp.m_U0hatTPreFFTScaledNK, ct);
    }

    if (isSparse) {
        z->EvalAddInPlace(r2c, z->EvalConjugateInC(r2c));
        z->ModReduceInPlace(r2c);

        if (scalingOption == SCALE_N || scalingOption == SCALE_NK) {
            if (scalingOption == SCALE_N) {
                scaleDown(r2c, scaleN);
            }
            else {
                scaleDown(r2c, scaleNK);
            }
        }
        return r2c;
    }
    else {
        auto r2cConj = z->EvalConjugateInC(r2c);
        auto r2cI    = z->EvalSub(r2c, r2cConj);
        z->EvalAddInPlace(r2c, r2cConj);
        cc->GetScheme()->MultByMonomialInPlace(r2cI, 3 * cSlots);
        z->ModReduceInPlace(r2c);
        z->ModReduceInPlace(r2cI);

        if (scalingOption == SCALE_N || scalingOption == SCALE_NK) {
            if (scalingOption == SCALE_N) {
                scaleDown(r2c, scaleN);
                scaleDown(r2cI, scaleN);
            }
            else {
                scaleDown(r2c, scaleNK);
                scaleDown(r2cI, scaleNK);
            }
        }
        return std::vector{r2c, r2cI};
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalC2Z(CiphertextGroup ct, C2ZOption c2zOption) const {
    auto cc       = ct[0]->GetCryptoContext();
    auto cSlots   = ct[0]->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    bool isSparse = precomp.m_isSparse;

    if (isSparse) {
        auto c2z =
            EvalZLinearTransform(c2zOption == C2Z_SPECIAL_A2AE ? precomp.m_ZUSpecialA2AePre : precomp.m_ZUPre, ct[0]);
        z->EvalAddInPlace(c2z, cc->EvalRotate(c2z, cSlots));
        z->ModReduceInPlace(c2z);
        return c2z;
    }
    else {
        auto c2z0 =
            EvalZLinearTransform(c2zOption == C2Z_SPECIAL_A2AE ? precomp.m_ZU0SpecialA2AePre : precomp.m_ZU0Pre, ct[0]);
        auto c2z1 =
            EvalZLinearTransform(c2zOption == C2Z_SPECIAL_A2AE ? precomp.m_ZU1SpecialA2AePre : precomp.m_ZU1Pre, ct[1]);
        z->EvalAddInPlace(c2z0, c2z1);
        z->ModReduceInPlace(c2z0);
        return c2z0;
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalZ2R(ConstCiphertext<DCRTPoly>& ct, Z2COption z2cOption) const {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto isSparse = precomp.m_isSparse;

    if (isSparse) {
        return EvalC2R(EvalZ2C(ct, z2cOption)[0]);
    }
    else {
        auto z2cGroup = EvalZ2C(ct, z2cOption);
        auto z2c0     = z2cGroup[0];
        auto z2c1     = z2cGroup[1];
        cc->GetScheme()->MultByMonomialInPlace(z2c1, cSlots);
        cc->EvalAddInPlaceNoCheck(z2c0, z2c1);
        return EvalC2R(z2c0);
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalR2Z(ConstCiphertext<DCRTPoly>& ct, R2CScalingOption scalingOption,
                                       C2ZOption c2zOption) const {
    return EvalC2Z(EvalR2C(ct, scalingOption), c2zOption);
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToArithHigh(ConstCiphertext<DCRTPoly>& ct) const {
    auto cc      = ct->GetCryptoContext();
    auto cSlots  = ct->GetZEncodingParams().getCSlots();
    auto precomp = GetBootPrecom(cSlots);

    auto z2r = EvalZ2R(ct, Z2C_NORMAL);

    //------------------------------------------------------------------------------
    // ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalModRaisePartialSum(z2r);

    //------------------------------------------------------------------------------
    // R-To-Z
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N. Carried out in high precision
    auto r2z = EvalR2Z(raised, R2CScalingOption::SCALE_N_PRE, C2Z_NORMAL);
    return r2z;
}

void FHEZImpl::ApplyDoubleAngleIterations(Ciphertext<DCRTPoly>& ct, uint32_t numIter) const {
    auto cc = ct->GetCryptoContext();
    for (uint32_t i = 0; i != numIter; ++i) {
        ct = z->EvalMult(ct, ct);
        z->ModReduceInPlace(ct);
        z->EvalAddInPlace(ct, z->EvalAddInC(ct, r_sparse_scalars[i]));
    }
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToArithNoise(ConstCiphertext<DCRTPoly>& ct) const {
    auto cc        = ct->GetCryptoContext();
    auto cSlots    = ct->GetZEncodingParams().getCSlots();
    auto precomp   = GetBootPrecom(cSlots);
    auto elemParam = ct->GetElements()[0].GetParams();

    auto z2r = EvalZ2R(ct, Z2C_SPECIAL_A2AE);

    //------------------------------------------------------------------------------
    // ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalModRaisePartialSum(z2r);
    // Now the message is multplied by N/(rN)
    //__heir_debug2(raised, "Raised");

    //------------------------------------------------------------------------------
    // R-To-C
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N * K; carried out in high precision
    auto r2c = EvalR2C(raised, R2CScalingOption::SCALE_NK_PRE);
    //pkeZ_global->debug(r2c[0], "R2C0");
    //__heir_debug2(r2c, "R2C");

    //------------------------------------------------------------------------------
    // Approximate Mod Reduction
    //------------------------------------------------------------------------------

    // Evaluate Chebyshev series for the sine wave
    auto& coeff_g0 = coeff_g0_big_complex_32;
    auto g0        = r2c.map([&](ConstCiphertext<DCRTPoly>& input) -> Ciphertext<DCRTPoly> {
        auto g0 = advZ->EvalChebyshevSeriesPS(input, coeff_g0);

        // Double-angle iterations
        uint32_t numIter = FHEZImpl::R_SPARSE;
        ApplyDoubleAngleIterations(g0, numIter);
        return g0;
    });

    //__heir_debug2(g0, "Sine");

    //------------------------------------------------------------------------------
    // C-To-Z
    //------------------------------------------------------------------------------

    auto r2z = EvalC2Z(g0, C2Z_SPECIAL_A2AE);

    //------------------------------------------------------------------------------
    // Subtract
    //------------------------------------------------------------------------------
    return z->EvalSubWithAdjust(ct, r2z);
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToArith(ConstCiphertext<DCRTPoly>& ct) const {
    // We first call resetI then clean noise, as large I will affect *linear transformation precision*!!!!
    auto resetICt   = EvalArithToArithHigh(ct);
    auto lowNoiseCt = EvalArithToArithNoise(resetICt);
    return lowNoiseCt;
}

Ciphertext<DCRTPoly> FHEZImpl::EvalArithToBooleanSparse(ConstCiphertext<DCRTPoly>& ct) {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    bool isSparse = precomp.m_isSparse;
    if (!isSparse) {
        OPENFHE_THROW("ArithToBooleanSparse called on non-sparse ciphertext");
    }

    // Our core ct
    auto core = EvalZ2C(ct, Z2C_SPECIAL_B0)[0];
    // TODO: We need to keep core ct at bottom; rescale if necessary

    auto zN     = ct->GetZEncodingParams().getZN();
    auto zSlots = ct->GetZEncodingParams().getZSlots();
    uint32_t w  = precomp.m_w;
    // We ask zN to be multiple of w now...
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));

    std::vector<Ciphertext<DCRTPoly>> msbs;

    // Iteratively process each low bits
    for (uint32_t iter = 0; iter != numIter; ++iter) {
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

        // core may change over time
        auto elemParam       = core->GetElements()[0].GetParams();
        auto sf              = core->GetScalingFactorBFP();
        Plaintext oneHotPtxt = ZEncodingImpl::encodeR(oneHotPoly, elemParam, sf);

        auto coreMasked = z->EvalMult(core, oneHotPtxt);
        z->ModReduceInPlace(coreMasked);

        auto lutGroup =
            internalBooleanToBooleanCustomTwoLUTSparse(coreMasked, precomp.m_lutIDCoeffs, precomp.m_lutMSBCoeffs);
        auto lut = lutGroup[0];
        auto msb = lutGroup[1];

        //------------------------------------------------------------------------------
        // Store the parts and remove it from core
        //------------------------------------------------------------------------------

        msbs.push_back(msb);

        // remove part from core
        for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
            int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
            // Note the rotation index is negative here
            int32_t rotationIndex = diff * w;
            if (rotationIndex <= precomp.m_cutoff) {
                // We do not remove them any more
                // Just treat the lower parts as noises
                break;
            }
            // Multiply by 1 / (2^{nextIter - iter} * p)
            auto scaled = z->EvalMultInC(lut, BigFixedPoint::pow2(diff * w));
            // If cross the half-way point, need to rotate more
            if (nextIter * w >= zN / 2 && iter * w < zN / 2) {
                rotationIndex -= static_cast<int32_t>((zSlots - 1) * zN / 2);
            }
            scaled = cc->EvalRotate(scaled, rotationIndex);
            z->ModReduceInPlace(scaled);
            z->EvalSubWithAdjustInPlace(core, scaled);
        }
        // Remove itself from core
        z->EvalSubWithAdjustInPlace(core, lut);
    }

    //------------------------------------------------------------------------------
    // Combine all parts
    //------------------------------------------------------------------------------

    for (size_t i = 1; i != msbs.size(); ++i) {
        // They may have different scaling factors...
        z->EvalAddInPlace(msbs[0], msbs[i]);
    }
    // Properly set the Z-encoding params
    msbs[0]->SetZEncodingParams(ZEncodingParams(BModeSparse, zN, zSlots));
    return msbs[0];
}

CiphertextGroup FHEZImpl::EvalArithToBooleanFull(ConstCiphertext<DCRTPoly>& ct) {
    auto cc       = ct->GetCryptoContext();
    auto cSlots   = ct->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    bool isSparse = precomp.m_isSparse;
    if (isSparse) {
        OPENFHE_THROW("ArithToBooleanFull called on sparse ciphertext");
    }

    // Our core ct
    auto cores = EvalZ2C(ct, Z2C_SPECIAL_B0);
    auto core0 = cores[0];
    auto core1 = cores[1];
    // TODO: We need to keep core ct at bottom; rescale if necessary

    auto zN     = ct->GetZEncodingParams().getZN();
    auto zSlots = ct->GetZEncodingParams().getZSlots();
    uint32_t w  = precomp.m_w;
    // We ask zN to be multiple of w now...
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));

    std::vector<Ciphertext<DCRTPoly>> msbs;

    // Iteratively process each low bits
    for (uint32_t iter = 0; iter != numIter; ++iter) {
        auto core       = core0;
        bool secondHalf = (iter * w >= zN / 2);
        if (secondHalf) {
            core = core1;
        }

        auto elemParam  = core->GetElements()[0].GetParams();
        auto sf         = core->GetScalingFactorBFP();
        auto maskPtxt   = getATBMaskFullPacking(iter, w, zN, zSlots, elemParam, sf);
        auto coreMasked = z->EvalMult(core, maskPtxt);
        z->ModReduceInPlace(coreMasked);

        auto lutGroup =
            internalBooleanToBooleanCustomTwoLUTSparse(coreMasked, precomp.m_lutIDCoeffs, precomp.m_lutMSBCoeffs);
        auto lut = lutGroup[0];
        auto msb = lutGroup[1];

        //------------------------------------------------------------------------------
        // Store the parts and remove it from core
        //------------------------------------------------------------------------------

        msbs.push_back(msb);

        // remove part from core
        for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
            int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
            // Note the rotation index is negative here
            int32_t rotationIndex = diff * w;
            if (rotationIndex <= precomp.m_cutoff) {
                // We do not remove them any more
                // Just treat the lower parts as noises
                break;
            }
            // Multiply by 1 / (2^{nextIter - iter} * p)
            auto scaled = z->EvalMultInC(lut, BigFixedPoint::pow2(diff * w));

            bool targetSecondHalf = (nextIter * w >= zN / 2);

            // If cross the half-way point, switch to second ciphertext
            auto targetCore = core;
            if (targetSecondHalf && !secondHalf) {
                rotationIndex += zN / 2;
                targetCore = core1;
            }
            if (rotationIndex != 0) {
                scaled = cc->EvalRotate(scaled, rotationIndex);
            }
            z->ModReduceInPlace(scaled);
            z->EvalSubWithAdjustInPlace(targetCore, scaled);
        }
        // Remove itself from core
        z->EvalSubWithAdjustInPlace(core, lut);
    }

    //------------------------------------------------------------------------------
    // Combine all parts
    //------------------------------------------------------------------------------

    auto msb0 = msbs[0];
    auto msb1 = msbs[msbs.size() / 2];
    for (size_t i = 1; i != msbs.size() / 2; ++i) {
        // They may have different scaling factors...
        z->EvalAddInPlace(msb0, msbs[i]);
        z->EvalAddInPlace(msb1, msbs[msbs.size() / 2 + i]);
    }
    // Properly set the Z-encoding params
    msb0->SetZEncodingParams(ZEncodingParams(BModeFull, zN, zSlots));
    msb1->SetZEncodingParams(ZEncodingParams(BModeFull, zN, zSlots));

    CiphertextGroup msbGroup({msb0, msb1});
    return msbGroup;
}

CiphertextGroup FHEZImpl::EvalArithToBooleanBatched(CiphertextGroup ctxts) {
    auto cc       = ctxts[0]->GetCryptoContext();
    auto cSlots   = ctxts[0]->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto isSparse = precomp.m_isSparse;

    auto zN     = ctxts[0]->GetZEncodingParams().getZN();
    auto zSlots = ctxts[0]->GetZEncodingParams().getZSlots();
    uint32_t w  = precomp.m_w;
    // We ask zN to be multiple of w now...
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));

    if (isSparse) {
        OPENFHE_THROW("Batched EvalArithToBooleanBatched only supports full packing");
    }
    if (ctxts.getParts().size() != numIter) {
        // TODO: support slightly smaller batch size..
        OPENFHE_THROW("Batch size mismatch for batched EvalArithToBooleanBatched");
    }

    auto core = ctxts.mapWide([&](ConstCiphertext<DCRTPoly>& ct) { return EvalZ2C(ct, Z2C_SPECIAL_B0); });

    std::vector<Ciphertext<DCRTPoly>> coreFirstHalfVec;
    std::vector<Ciphertext<DCRTPoly>> coreSecondHalfVec;
    for (size_t i = 0; i != ctxts.getParts().size(); ++i) {
        coreFirstHalfVec.push_back(core[2 * i]);
        coreSecondHalfVec.push_back(core[2 * i + 1]);
    }
    // Now arrange first halves (and second halves) with proper rotations so that they can be processed in batch
    // With example parameter of w = 4, zN = 32, numIter = 8, the first halves will be arranged as below, where each row is a ciphertext and each column corresponds to the same part (same low bits):
    // First halves:
    // ----
    //  ----
    //   ----
    //    ----
    // ----
    //  ----
    //   ----
    //    ----

    for (size_t i = 1; i != ctxts.getParts().size() / 2; ++i) {
        int32_t rotationIndex = static_cast<int32_t>(-i * w);
        coreFirstHalfVec[i]   = cc->EvalRotate(coreFirstHalfVec[i], rotationIndex);
        coreFirstHalfVec[i + ctxts.getParts().size() / 2] =
            cc->EvalRotate(coreFirstHalfVec[i + ctxts.getParts().size() / 2], rotationIndex);

        coreSecondHalfVec[i] = cc->EvalRotate(coreSecondHalfVec[i], rotationIndex);
        coreSecondHalfVec[i + ctxts.getParts().size() / 2] =
            cc->EvalRotate(coreSecondHalfVec[i + ctxts.getParts().size() / 2], rotationIndex);

        // Update core as well; this is pointer update...
        // We must update it other wise later subtraction will be wrong
        core[2 * i]                                     = coreFirstHalfVec[i];
        core[2 * i + 1]                                 = coreSecondHalfVec[i];
        core[2 * (i + ctxts.getParts().size() / 2)]     = coreFirstHalfVec[i + ctxts.getParts().size() / 2];
        core[2 * (i + ctxts.getParts().size() / 2) + 1] = coreSecondHalfVec[i + ctxts.getParts().size() / 2];
    }

    CiphertextGroup coreFirstHalf(coreFirstHalfVec);
    CiphertextGroup coreSecondHalf(coreSecondHalfVec);

    std::vector<Ciphertext<DCRTPoly>> msbs;

    // Iteratively process each low bits
    for (uint32_t iter = 0; iter != numIter; ++iter) {
        bool secondHalf = (iter * w >= zN / 2);

        CiphertextGroup targetGroup = coreFirstHalf;
        if (secondHalf) {
            targetGroup = coreSecondHalf;
        }

        // core may change over time
        auto elemParam = targetGroup[0]->GetElements()[0].GetParams();
        auto sf        = targetGroup[0]->GetScalingFactorBFP();

        std::vector<Ciphertext<DCRTPoly>> masked(targetGroup.size());
        for (size_t i = 0; i != targetGroup.size(); ++i) {
            // Here, the mask is (iter + i % (halfSize)), for accounting the previous rotations.
            auto maskPtxt = getATBMask(iter + (i % (targetGroup.size() / 2)), w, zN, zSlots, elemParam, sf);
            auto maskedCt = z->EvalMult(targetGroup[i], maskPtxt);
            z->ModReduceInPlace(maskedCt);
            masked[i] = maskedCt;
        }

        auto targetCombined0 = masked[0];
        auto targetCombined1 = masked[masked.size() / 2];
        for (size_t i = 1; i != masked.size() / 2; ++i) {
            z->EvalAddInPlace(targetCombined0, masked[i]);
            z->EvalAddInPlace(targetCombined1, masked[masked.size() / 2 + i]);
        }

        auto lut  = internalBooleanToBooleanCustomTwoLUTFull(std::vector{targetCombined0, targetCombined1},
                                                             precomp.m_lutIDCoeffs, precomp.m_lutMSBCoeffs);
        auto lut0 = lut[0];
        auto lut1 = lut[1];
        auto msb0 = lut[2];
        auto msb1 = lut[3];

        //------------------------------------------------------------------------------
        // Store the parts and remove it from core
        //------------------------------------------------------------------------------

        msbs.push_back(msb0);
        msbs.push_back(msb1);

        elemParam = lut0->GetElements()[0].GetParams();
        sf        = lut0->GetScalingFactorBFP();

        // remove part from core
        for (uint32_t nextIter = iter + 1; nextIter != numIter; ++nextIter) {
            int32_t diff = static_cast<int32_t>(iter) - static_cast<int32_t>(nextIter);
            // Note the rotation index is negative here
            int32_t rotationIndex = diff * w;
            if (rotationIndex <= precomp.m_cutoff) {
                // We do not remove them any more
                // Just treat the lower parts as noises
                break;
            }
            if (nextIter * w >= zN / 2 && iter * w < zN / 2) {
                rotationIndex += zN / 2;
            }
            for (size_t j = 0; j != targetGroup.size() / 2; ++j) {
                // Multiply by 1 / (2^{nextIter - iter} * p)
                auto scaleDown = BigFixedPoint::pow2(diff * w);
                // iter + j same reason as above for maskPtxt
                auto subtractMaskPtxt =
                    getATBSubtractMaskFullPacking(iter + j, w, zN, zSlots, elemParam, sf, scaleDown);

                {
                    auto scaled = z->EvalMult(lut0, subtractMaskPtxt);
                    scaled      = cc->EvalRotate(scaled, rotationIndex);
                    z->ModReduceInPlace(scaled);

                    bool targetSecondHalf = (nextIter * w >= zN / 2);

                    // Here we used pointer....changes to core will be reflected in coreGroup
                    auto coreCt = core[2 * j + targetSecondHalf];
                    z->EvalSubWithAdjustInPlace(coreCt, scaled);
                }
                {
                    auto scaled = z->EvalMult(lut1, subtractMaskPtxt);
                    scaled      = cc->EvalRotate(scaled, rotationIndex);
                    z->ModReduceInPlace(scaled);

                    bool targetSecondHalf = (nextIter * w >= zN / 2);
                    auto coreCt           = core[2 * (targetGroup.size() / 2 + j) + targetSecondHalf];
                    z->EvalSubWithAdjustInPlace(coreCt, scaled);
                }
            }
        }
        // No need to remove itself from core unlike the sparse case
    }

    //------------------------------------------------------------------------------
    // Combine all MSBs
    //------------------------------------------------------------------------------

    std::vector<Ciphertext<DCRTPoly>> finalMSBs;
    std::vector<Ciphertext<DCRTPoly>> finalMSBsSec;
    for (size_t j = 0; j != ctxts.size() / 2; ++j) {
        Ciphertext<DCRTPoly> MSBjFirst, MSBjSecond, MSBSecjFirst, MSBSecjSecond;
        for (size_t iter = 0; iter != numIter / 2; ++iter) {
            auto mask = getATBMask(iter + j, w, zN, zSlots, msbs[0]->GetElements()[0].GetParams(),
                                   msbs[0]->GetScalingFactorBFP());
            {
                auto maskedMSB = z->EvalMult(msbs[2 * iter], mask);
                z->ModReduceInPlace(maskedMSB);
                if (iter == 0) {
                    MSBjFirst = maskedMSB;
                }
                else {
                    z->EvalAddInPlace(MSBjFirst, maskedMSB);
                }
            }
            {
                auto maskedMSB = z->EvalMult(msbs[2 * iter + 1], mask);
                z->ModReduceInPlace(maskedMSB);
                if (iter == 0) {
                    MSBSecjFirst = maskedMSB;
                }
                else {
                    z->EvalAddInPlace(MSBSecjFirst, maskedMSB);
                }
            }
        }
        if (j != 0) {
            MSBjFirst    = cc->EvalRotate(MSBjFirst, static_cast<int32_t>(j * w));
            MSBSecjFirst = cc->EvalRotate(MSBSecjFirst, static_cast<int32_t>(j * w));
        }
        for (size_t iter = numIter / 2; iter != numIter; ++iter) {
            auto mask = getATBMask(iter + j, w, zN, zSlots, msbs[0]->GetElements()[0].GetParams(),
                                   msbs[0]->GetScalingFactorBFP());
            {
                auto maskedMSB = z->EvalMult(msbs[2 * iter], mask);
                z->ModReduceInPlace(maskedMSB);
                if (iter == numIter / 2) {
                    MSBjSecond = maskedMSB;
                }
                else {
                    z->EvalAddInPlace(MSBjSecond, maskedMSB);
                }
            }
            {
                auto maskedMSB = z->EvalMult(msbs[2 * iter + 1], mask);
                z->ModReduceInPlace(maskedMSB);
                if (iter == numIter / 2) {
                    MSBSecjSecond = maskedMSB;
                }
                else {
                    z->EvalAddInPlace(MSBSecjSecond, maskedMSB);
                }
            }
        }
        if (j != 0) {
            MSBjSecond    = cc->EvalRotate(MSBjSecond, static_cast<int32_t>(j * w));
            MSBSecjSecond = cc->EvalRotate(MSBSecjSecond, static_cast<int32_t>(j * w));
        }
        finalMSBs.push_back(MSBjFirst);
        finalMSBs.push_back(MSBjSecond);
        finalMSBsSec.push_back(MSBSecjFirst);
        finalMSBsSec.push_back(MSBSecjSecond);
    }

    finalMSBs.insert(finalMSBs.end(), finalMSBsSec.begin(), finalMSBsSec.end());

    // Properly set the Z-encoding params
    for (auto ct : finalMSBs) {
        ct->SetZEncodingParams(ZEncodingParams(BModeFull, zN, zSlots));
    }
    return finalMSBs;
}

Ciphertext<DCRTPoly> FHEZImpl::internalBooleanToBooleanLTsSparse(ConstCiphertext<DCRTPoly>& ct) const {
    auto cc      = ct->GetCryptoContext();
    auto cSlots  = ct->GetZEncodingParams().getCSlots();
    auto precomp = GetBootPrecom(cSlots);

    auto c2r = EvalC2R(ct);

    //------------------------------------------------------------------------------
    // ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalModRaisePartialSum(c2r);
    // Now the message is multplied by N/(rN)

    //------------------------------------------------------------------------------
    // R-To-C
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N * K because of scaling in pre-compute
    // can do so in low precision as we do not need high precision here
    // Since we call C2R without imaginary part, we only need the real part here
    auto r2c = EvalR2C(raised, SCALE_NK_PRE)[0];
    return r2c;
}

CiphertextGroup FHEZImpl::internalBooleanToBooleanLTsFull(CiphertextGroup ct) const {
    auto cc      = ct[0]->GetCryptoContext();
    auto cSlots  = ct[0]->GetZEncodingParams().getCSlots();
    auto precomp = GetBootPrecom(cSlots);

    // Put another ct in imaginary part
    auto ctComb = ct[0]->Clone();
    auto ct1I   = cc->GetScheme()->MultByMonomial(ct[1], cSlots);
    z->EvalAddInPlace(ctComb, ct1I);
    auto c2r = EvalC2R(ctComb);

    //------------------------------------------------------------------------------
    // ModRaise and PartialSum
    //------------------------------------------------------------------------------

    auto raised = EvalModRaisePartialSum(c2r);

    //------------------------------------------------------------------------------
    // R-To-C
    //------------------------------------------------------------------------------

    // R2C then will multiply by rN then divide by N * K because of scaling in pre-compute
    // can do so in low precision as we do not need high precision here
    auto r2cGroup = EvalR2C(raised, SCALE_NK_PRE);
    return r2cGroup;
}

CiphertextGroup FHEZImpl::internalBooleanToBooleanCustomTwoLUTSparse(ConstCiphertext<DCRTPoly>& ct,
                                                                     const std::vector<BigComplex>& lutCoeffs,
                                                                     const std::vector<BigComplex>& lutCoeffs2) const {
    if (lutCoeffs.size() != lutCoeffs2.size()) {
        OPENFHE_THROW("MVB LUT size mismatch");
    }

    auto r2c = internalBooleanToBooleanLTsSparse(ct);

    //------------------------------------------------------------------------------
    // Exp
    //------------------------------------------------------------------------------

    auto& coeff_exp = coeff_exp_16_big_complex_46;
    auto res        = advZ->EvalChebyshevSeriesPS(r2c, coeff_exp);

    // Double angle-iterations to get exp(2*Pi*i*x)
    res = z->EvalSquare(res);
    z->ModReduceInPlace(res);
    res = z->EvalSquare(res);
    z->ModReduceInPlace(res);

    //------------------------------------------------------------------------------
    // Running LUT
    //------------------------------------------------------------------------------

    auto powers = advZ->EvalPowers(res, lutCoeffs);
    auto lut    = advZ->EvalPolyWithPrecomp(powers, lutCoeffs);
    // Take the real part
    z->EvalAddInPlace(lut, z->EvalConjugateInC(lut));

    auto lut2 = advZ->EvalPolyWithPrecomp(powers, lutCoeffs2);
    // Take the real part
    z->EvalAddInPlace(lut2, z->EvalConjugateInC(lut2));
    return std::vector{lut, lut2};
}

CiphertextGroup FHEZImpl::internalBooleanToBooleanCustomTwoLUTFull(CiphertextGroup ctGroup,
                                                                   const std::vector<BigComplex>& lutCoeffs,
                                                                   const std::vector<BigComplex>& lutCoeffs2) const {
    if (lutCoeffs.size() != lutCoeffs2.size()) {
        OPENFHE_THROW("MVB LUT size mismatch");
    }

    auto r2c = internalBooleanToBooleanLTsFull(ctGroup);

    std::vector<std::shared_ptr<seriesPowers<DCRTPoly>>> powersVec;

    auto expCt = r2c.map([&](ConstCiphertext<DCRTPoly> ct) {
        //------------------------------------------------------------------------------
        // Exp
        //------------------------------------------------------------------------------

        auto& coeff_exp = coeff_exp_16_big_complex_46;
        auto res        = advZ->EvalChebyshevSeriesPS(ct, coeff_exp);

        // Double angle-iterations to get exp(2*Pi*i*x)
        res = z->EvalSquare(res);
        z->ModReduceInPlace(res);
        res = z->EvalSquare(res);
        z->ModReduceInPlace(res);

        auto powers = advZ->EvalPowers(res, lutCoeffs);
        powersVec.push_back(powers);
        // Unused
        return res;
    });

    //------------------------------------------------------------------------------
    // Running LUT
    //------------------------------------------------------------------------------

    std::vector<Ciphertext<DCRTPoly>> resVec;

    for (size_t i = 0; i != 2; ++i) {
        auto lut = advZ->EvalPolyWithPrecomp(powersVec[i], lutCoeffs);
        // Take the real part
        z->EvalAddInPlace(lut, z->EvalConjugateInC(lut));
        resVec.push_back(lut);
    }

    for (size_t i = 0; i != 2; ++i) {
        auto lut = advZ->EvalPolyWithPrecomp(powersVec[i], lutCoeffs2);
        // Take the real part
        z->EvalAddInPlace(lut, z->EvalConjugateInC(lut));
        resVec.push_back(lut);
    }
    return resVec;
}

Ciphertext<DCRTPoly> FHEZImpl::EvalBooleanToBooleanSparse(ConstCiphertext<DCRTPoly>& ct) const {
    auto ctNew = ct->Clone();
    // Double scaling factor to divide by 2
    ctNew->SetScalingFactorBFP(ctNew->GetScalingFactorBFP() * BigFixedPoint::two());
    auto r2c = internalBooleanToBooleanLTsSparse(ctNew);

    //------------------------------------------------------------------------------
    // Cosine
    //------------------------------------------------------------------------------

    auto& coeff_cos = coeff_cos_16_big_complex_50;
    auto res        = advZ->EvalChebyshevSeriesPS(r2c, coeff_cos);

    // Double angle-iterations to get cos(pi*x)
    res = z->EvalSquare(res);
    z->ModReduceInPlace(res);
    z->EvalAddInPlace(res, res);
    z->EvalAddInPlaceInC(res, -BigFixedPoint::one());  // cos(pi x)
    res = z->EvalSquare(res);
    z->ModReduceInPlace(res);  // cos^2(pi x)

    //------------------------------------------------------------------------------
    // The LUT (1 - cos^2(pi x))
    // This is p == 2 and order == 1 of AKP25 for the identity map
    // f(0) = 0, f(1) = 1
    //------------------------------------------------------------------------------

    z->EvalNegateInPlace(res);                        // -cos^2(pi x)
    z->EvalAddInPlaceInC(res, BigFixedPoint::one());  // 1 - cos^2(pi x)

    return res;
}

CiphertextGroup FHEZImpl::EvalBooleanToBooleanFull(CiphertextGroup ct) const {
    auto ctNew0 = ct[0]->Clone();
    auto ctNew1 = ct[1]->Clone();
    // Double scaling factor to divide by 2
    ctNew0->SetScalingFactorBFP(ctNew0->GetScalingFactorBFP() * BigFixedPoint::two());
    ctNew1->SetScalingFactorBFP(ctNew1->GetScalingFactorBFP() * BigFixedPoint::two());
    auto r2c = internalBooleanToBooleanLTsFull(std::vector{ctNew0, ctNew1});

    return r2c.map([&](ConstCiphertext<DCRTPoly> ct) {
        //------------------------------------------------------------------------------
        // Cosine
        //------------------------------------------------------------------------------

        auto& coeff_cos = coeff_cos_16_big_complex_50;
        auto res        = advZ->EvalChebyshevSeriesPS(ct, coeff_cos);

        // Double angle-iterations to get cos(pi*x)
        res = z->EvalSquare(res);
        z->ModReduceInPlace(res);
        z->EvalAddInPlace(res, res);
        z->EvalAddInPlaceInC(res, -BigFixedPoint::one());  // cos(pi x)
        res = z->EvalSquare(res);
        z->ModReduceInPlace(res);  // cos^2(pi x)

        //------------------------------------------------------------------------------
        // The LUT (1 - cos^2(pi x))
        // This is p == 2 and order == 1 of AKP25 for the identity map
        // f(0) = 0, f(1) = 1
        //------------------------------------------------------------------------------

        z->EvalNegateInPlace(res);                        // -cos^2(pi x)
        z->EvalAddInPlaceInC(res, BigFixedPoint::one());  // 1 - cos^2(pi x)

        return res;
    });
}

Ciphertext<DCRTPoly> FHEZImpl::EvalBooleanToArith(CiphertextGroup ct) const {
    // which will multiply by t^{-1} in Z
    auto res = EvalC2Z(ct, C2Z_SPECIAL_A2AE);

    // Properly set the Z-encoding params
    auto zN     = ct[0]->GetZEncodingParams().getZN();
    auto zSlots = ct[0]->GetZEncodingParams().getZSlots();
    res->SetZEncodingParams(ZEncodingParams(ZMode, zN, zSlots));
    return res;
}

CiphertextGroup FHEZImpl::EvalArithToBoolean(CiphertextGroup ct) {
    auto cSlots   = ct[0]->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto zN       = ct[0]->GetZEncodingParams().getZN();
    auto isSparse = precomp.m_isSparse;

    // For batched
    uint32_t w = precomp.m_w;
    // We ask zN to be multiple of w now...
    uint32_t numIter = static_cast<uint32_t>(std::ceil(static_cast<double>(zN) / (static_cast<double>(w))));

    if (isSparse) {
        auto res = EvalArithToBooleanSparse(ct[0]);
        return CiphertextGroup({res});
    }
    if (ct.size() == 1) {
        auto res = EvalArithToBooleanFull(ct[0]);
        return res;
    }
    else {
        if (ct.size() > numIter) {
            // We can do batched EvalArithToBoolean
            OPENFHE_THROW("Batch size mismatch for batched EvalArithToBoolean");
        }
        if (ct.size() < numIter) {
            // Pad with dummy ciphertexts
            std::cout
                << "Warning: Padding ciphertexts for batched EvalArithToBoolean (TODO: implement batched mode for less ciphertexts)"
                << std::endl;
            auto dummyCt = ct[0]->Clone();
            std::vector<Ciphertext<DCRTPoly>> newCts;
            for (size_t i = 0; i != ct.size(); ++i) {
                newCts.push_back(ct[i]);
            }
            for (size_t i = ct.size(); i != numIter; ++i) {
                newCts.push_back(dummyCt);
            }
            ct = CiphertextGroup(newCts);
        }
        auto res = EvalArithToBooleanBatched(ct);
        return res;
    }
}

CiphertextGroup FHEZImpl::EvalBooleanToBoolean(CiphertextGroup ct) const {
    auto cSlots   = ct[0]->GetZEncodingParams().getCSlots();
    auto precomp  = GetBootPrecom(cSlots);
    auto isSparse = precomp.m_isSparse;

    if (isSparse) {
        auto res = EvalBooleanToBooleanSparse(ct[0]);
        return CiphertextGroup({res});
    }
    else {
        if (ct.size() != 2) {
            OPENFHE_THROW("Full packing BooleanToBoolean requires 2 ciphertexts");
        }
        auto res = EvalBooleanToBooleanFull(ct);
        return res;
    }
}

}  // namespace lbcrypto