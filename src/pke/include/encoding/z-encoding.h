#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_H_

#include <cassert>
#include "ciphertext-fwd.h"
#include "math/z-polynomial.h"
#include "plaintext.h"

namespace lbcrypto {

class ZEncodingImpl;

using ZEncoding = std::shared_ptr<ZEncodingImpl>;

// This is only interplay between Plaintext and RPolynomial
// RPolynomial to CSlots to ZPolynomial is handled elsewhere
class ZEncodingImpl : public PlaintextImpl {
private:
    BigFixedPoint m_scalingFactorBFP;
    ZEncodingParams m_zEncodingParams;

public:
    ZEncodingImpl(std::shared_ptr<DCRTPoly::Params> vp, const DCRTPoly& elements, BigFixedPoint scalingFactorBFP,
                  const ZEncodingParams& params)
        : PlaintextImpl(vp, nullptr, INVALID_ENCODING, INVALID_SCHEME) {
        encodedVectorDCRT  = elements;
        m_scalingFactorBFP = scalingFactorBFP;
        m_zEncodingParams  = params;
    }

    // Just to fullfill the PlaintextImpl interface
    virtual bool Encode() override {
        OPENFHE_THROW("not implemented");
    }
    virtual bool Decode() override {
        OPENFHE_THROW("not implemented");
    }
    virtual void PrintValue(std::ostream& out) const override {
        OPENFHE_THROW("Not implemented");
    }
    virtual bool CompareTo(const PlaintextImpl& other) const override {
        OPENFHE_THROW("Not implemented");
    }
    virtual size_t GetLength() const override {
        return m_zEncodingParams.getLength();
    }
    BigFixedPoint GetScalingFactorBFP() const {
        return m_scalingFactorBFP;
    }
    ZEncodingParams GetZEncodingParams() const {
        return m_zEncodingParams;
    }

    static ZEncoding encodeR(const RPolynomial& input, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor);

    static RPolynomial decodeR(ZEncoding input) {
        auto dcrtPoly      = input->GetElement<DCRTPoly>();
        auto outputSize    = input->GetLength();
        auto scalingFactor = input->GetScalingFactorBFP();
        dcrtPoly.SetFormat(Format::COEFFICIENT);
        auto bigPoly = dcrtPoly.CRTInterpolate();
        std::vector<BigFixedPoint> output;
        auto ringDim = dcrtPoly.GetParams()->GetRingDimension();
        auto q       = dcrtPoly.GetParams()->GetModulus();
        for (size_t i = 0; i != outputSize; ++i) {
            auto valInteger = bigPoly[i * (ringDim / outputSize)];
            bool neg        = false;
            if (valInteger > q.DividedBy(2)) {
                neg        = true;
                valInteger = q - valInteger;
            }
            auto valFixedPoint = BigFixedPoint(valInteger, 0, neg).scaleTo(128);
            output.push_back(valFixedPoint / scalingFactor);
        }
        return RPolynomial(input->GetZEncodingParams(), output);
    }

    static std::vector<uint64_t> decodeArithSmall(ZEncoding input);

    static ZEncoding encodeC(const BigComplex& input, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor);

    static ZEncoding encodeC(const CSlots& cSlots, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor) {
        return encodeR(cSlots.toRPolynomial(), elementParams, scalingFactor);
    }

    static ZEncoding encodeZ(std::vector<ZPolynomial> input, uint32_t zN, uint32_t zSlots,
                             const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                             const BigFixedPoint& scalingFactor) {
        if (input.size() < zSlots) {
            input.resize(zSlots, ZPolynomial::encodeZeros(zN));
        }
        std::vector<BigComplex> mergedSlots(zSlots * (zN / 2));
        //#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(zSlots))
        for (size_t i = 0; i != zSlots; ++i) {
            auto singleCSlots = input[i].toCSlots();
            for (size_t j = 0; j != zN / 2; ++j) {
                mergedSlots[i * (zN / 2) + j] = singleCSlots[j];
            }
        }
        ZEncodingParams params(ZMode, zN, zSlots);
        CSlots mergedCSlots(params, mergedSlots);
        return encodeC(mergedCSlots, elementParams, scalingFactor);
    }

    static ZEncoding encodeArith(std::vector<uint64_t> input, uint32_t zN, uint32_t zSlots,
                                 const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                 const BigFixedPoint& scalingFactor);

    // This is too slow
    static ZEncoding encodeArith(std::vector<BigInteger> input, uint32_t zN, uint32_t zSlots,
                                 const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                 const BigFixedPoint& scalingFactor) {
        std::vector<ZPolynomial> zPolys;
        for (size_t i = 0; i != input.size(); ++i) {
            zPolys.push_back(ZPolynomial::encode(zN, input[i]));
        }
        if (zPolys.size() < zSlots) {
            zPolys.resize(zSlots, ZPolynomial::encodeZeros(zN));
        }
        return encodeZ(zPolys, zN, zSlots, elementParams, scalingFactor);
    }
    //

    static ZEncoding encodeArithSingle(BigInteger input, uint32_t zN,
                                       const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                       const BigFixedPoint& scalingFactor) {
        // Here zSlots = 1 can be automatically broadcasted to actual zSlots by sparse packing
        return encodeArith({input}, zN, 1, elementParams, scalingFactor);
    }

    static ZEncoding encodeTInZ(uint32_t zN, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                const BigFixedPoint& scalingFactor) {
        std::vector<ZPolynomial> zPolys(1, ZPolynomial::getT(zN));
        return encodeZ(zPolys, zN, 1, elementParams, scalingFactor);
    }

    static ZEncoding encodeTInvInZ(uint32_t zN, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                   const BigFixedPoint& scalingFactor) {
        std::vector<ZPolynomial> zPolys(1, ZPolynomial::getTInv(zN));
        return encodeZ(zPolys, zN, 1, elementParams, scalingFactor);
    }

    static std::vector<ZEncoding> encodeBooleanFull(std::vector<BigInteger> input, uint32_t zN, uint32_t zSlots,
                                                    const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                                    const BigFixedPoint& scalingFactor) {
        if (input.size() < zSlots) {
            input.resize(zSlots, 0);
        }
        auto N = elementParams->GetRingDimension();
        if (zSlots * zN != N) {
            OPENFHE_THROW("zSlots mismatch with ring dimension in encodeBooleanFull");
        }
        std::vector<BigComplex> lowHalf, highHalf;
        for (size_t i = 0; i != input.size(); ++i) {
            auto zPoly = ZPolynomial::encodeBinary(zN, input[i]);
            for (size_t i = 0; i != zN / 2; ++i) {
                lowHalf.push_back(zPoly[i]);
                highHalf.push_back(zPoly[i + zN / 2]);
            }
        }
        ZEncodingParams params(BModeFull, zN, zSlots);
        auto lowHalfEncoded  = encodeC(CSlots(params, lowHalf), elementParams, scalingFactor);
        auto highHalfEncoded = encodeC(CSlots(params, highHalf), elementParams, scalingFactor);
        return {lowHalfEncoded, highHalfEncoded};
    }

    static ZEncoding encodeBooleanSparse(std::vector<BigInteger> input, uint32_t zN, uint32_t zSlots,
                                         const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
                                         const BigFixedPoint& scalingFactor) {
        if (input.size() < zSlots) {
            input.resize(zSlots, 0);
        }
        auto N = elementParams->GetRingDimension();
        if (zSlots * zN > N / 2) {
            OPENFHE_THROW("zSlots mismatch with ring dimension in encodeBooleanSparse");
        }
        std::vector<BigComplex> lowHalf, highHalf;
        for (size_t i = 0; i != input.size(); ++i) {
            auto zPoly = ZPolynomial::encodeBinary(zN, input[i]);
            for (size_t i = 0; i != zN / 2; ++i) {
                lowHalf.push_back(zPoly[i]);
                highHalf.push_back(zPoly[i + zN / 2]);
            }
        }
        // append high half to low half
        lowHalf.insert(lowHalf.end(), highHalf.begin(), highHalf.end());
        ZEncodingParams params(BModeSparse, zN, zSlots);
        return encodeC(CSlots(params, lowHalf), elementParams, scalingFactor);
    }

    static std::vector<ZEncoding> encodeBooleanFullSingle(
        BigInteger input, uint32_t zN, const std::shared_ptr<typename DCRTPoly::Params>& elementParams,
        const BigFixedPoint& scalingFactor) {
        std::vector<BigComplex> lowHalf, highHalf;
        auto zPoly = ZPolynomial::encodeBinary(zN, input);
        for (size_t i = 0; i != zN / 2; ++i) {
            lowHalf.push_back(zPoly[i]);
            highHalf.push_back(zPoly[i + zN / 2]);
        }
        ZEncodingParams params(BModeFull, zN, 1);
        auto lowHalfEncoded  = encodeC(CSlots(params, lowHalf), elementParams, scalingFactor);
        auto highHalfEncoded = encodeC(CSlots(params, highHalf), elementParams, scalingFactor);
        return {lowHalfEncoded, highHalfEncoded};
    }
};

}  // namespace lbcrypto

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_H_