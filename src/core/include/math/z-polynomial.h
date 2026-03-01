#ifndef SRC_PKE_EXAMPLES_Z_ENCODE_UTILS_H_
#define SRC_PKE_EXAMPLES_Z_ENCODE_UTILS_H_

#include <cassert>
#include "math/z-constants.h"

namespace lbcrypto {

enum ZEncodingType {
    INVALID = 0,
    ZMode,        // Representing ZPolynomial, at arithmetic mode
    BModeSparse,  // Representing bit vector, at boolean mode in sparse packing
    BModeFull,    // Representing bit vector, at boolean mode in full packing
    CMode,        // Representing Complex values. Used at bootstrapping and C
};

std::ostream& operator<<(std::ostream& os, const ZEncodingType& type);

struct ZEncodingParams {
private:
    ZEncodingType m_encodingType;

    // Meaningful only for ZPolynomial/Boolean representation
    uint32_t m_zN;
    uint32_t m_zSlots;
    // Meaningful only for C mode
    uint32_t m_rN;
    uint32_t m_cSlots;

public:
    ZEncodingParams() : m_encodingType(INVALID) {}
    ZEncodingParams(const ZEncodingParams&) = default;

    ZEncodingParams(ZEncodingType type, uint32_t N, uint32_t slots = 0) {
        m_encodingType = type;
        if (type == ZMode || type == BModeSparse || type == BModeFull) {
            m_zN     = N;
            m_zSlots = slots;
            m_rN     = 0;
            m_cSlots = 0;
        }
        if (type == CMode) {
            m_zN     = 0;
            m_zSlots = 0;
            m_rN     = N;
            m_cSlots = N / 2;
        }
    }

    uint32_t getLength() const {
        if (m_encodingType == ZMode || m_encodingType == BModeFull) {
            return m_zN * m_zSlots;
        }
        if (m_encodingType == BModeSparse) {
            // all bits in complex slots. Number of coefficients is twice the number of complex slots.
            return m_zN * m_zSlots * 2;
        }
        if (m_encodingType == CMode) {
            return m_rN;
        }
        OPENFHE_THROW("getLength called for invalid encoding");
    }

    uint32_t getSlots() const {
        if (m_encodingType == ZMode || m_encodingType == BModeFull || m_encodingType == BModeSparse) {
            return m_zSlots;
        }
        else {
            return m_cSlots;
        }
    }

    uint32_t getCSlots() const {
        if (m_encodingType == ZMode || m_encodingType == BModeFull) {
            // Equivalent number of C slots
            // Useful inside bootstrapping
            return m_zSlots * m_zN / 2;
        }
        else if (m_encodingType == BModeSparse) {
            return m_zSlots * m_zN / 2;
        }
        else {
            return m_cSlots;
        }
    }

    uint32_t getZN() const {
        if (m_encodingType == ZMode || m_encodingType == BModeFull || m_encodingType == BModeSparse) {
            return m_zN;
        }
        else {
            OPENFHE_THROW("getZN called for CMode");
        }
    }

    uint32_t getZSlots() const {
        if (m_encodingType == ZMode || m_encodingType == BModeFull || m_encodingType == BModeSparse) {
            return m_zSlots;
        }
        else {
            OPENFHE_THROW("getZSlots called for CMode");
        }
    }

    bool isZMode() const {
        return m_encodingType == ZMode;
    }
    bool isCMode() const {
        return m_encodingType == CMode;
    }
    bool isBModeSparse() const {
        return m_encodingType == BModeSparse;
    }
    bool isBModeFull() const {
        return m_encodingType == BModeFull;
    }

    ZEncodingType getEncodingType() const {
        return m_encodingType;
    }

    bool compatible(ZEncodingParams rhs) {
        if (m_encodingType == rhs.m_encodingType) {
            if (m_encodingType == ZMode && m_zN == rhs.getZN() && m_zSlots == rhs.getZSlots()) {
                return true;
            }
            if (m_encodingType == CMode) {  // in CMode, we may take advantage of trace/repeating
                return true;
            }
        }
        return false;
    }

    std::string toString() const {
        std::string result;
        if (m_encodingType == ZMode) {
            result = "ZMode: zN=" + std::to_string(m_zN) + " zSlots=" + std::to_string(m_zSlots);
        }
        else if (m_encodingType == BModeSparse) {
            result = "BModeSparse: zN=" + std::to_string(m_zN) + " zSlots=" + std::to_string(m_zSlots);
        }
        else if (m_encodingType == BModeFull) {
            result = "BModeFull: zN=" + std::to_string(m_zN) + " zSlots=" + std::to_string(m_zSlots);
        }
        else if (m_encodingType == CMode) {
            result = "CMode: rN=" + std::to_string(m_rN) + " cSlots=" + std::to_string(m_cSlots);
        }
        else {
            result = "Invalid Encoding";
        }
        return result;
    }

    template <class Archive>
    void serialize(Archive& ar) {
        ar(m_encodingType, m_zN, m_zSlots, m_rN, m_cSlots);
    }
};

struct RPolynomial;
struct CSlots;

struct ZPolynomial {
public:
    ZPolynomial(uint32_t zN) : zN(zN), coefficients(zN) {}
    ZPolynomial(const std::vector<BigFixedPoint>& coeffs) : zN(coeffs.size()), coefficients(coeffs) {}
    std::vector<BigFixedPoint> getCoefficients() const {
        return coefficients;
    }

    uint32_t getZN() const {
        return zN;
    }

    double getLog2Norm() {
        double maxNorm = -1e30;
        for (const auto& coeff : coefficients) {
            double coeffNorm = coeff.log2Norm();
            if (coeffNorm > maxNorm) {
                maxNorm = coeffNorm;
            }
        }
        return maxNorm;
    }

    BigFixedPoint& operator[](size_t index) {
        return coefficients[index];
    }

    const BigFixedPoint& operator[](size_t index) const {
        return coefficients[index];
    }

    static ZPolynomial getT(uint32_t zN) {
        std::vector<BigFixedPoint> t(zN, BigFixedPoint::zero());
        auto two = BigFixedPoint::two();
        auto one = BigFixedPoint::one();
        t[0]     = -two;
        t[1]     = one;
        return ZPolynomial(t);
    }

    static ZPolynomial getTInv(uint32_t zN) {
        std::vector<BigFixedPoint> tInv;
        // denominator
        auto one = BigFixedPoint::one();
        auto d   = BigFixedPoint::pow2(zN);
        for (size_t i = 0; i != zN; ++i) {
            auto powerOf2 = BigFixedPoint::pow2(zN - 1 - i);
            if (i == 0) {
                tInv.push_back((one - powerOf2) / d);
            }
            else {
                tInv.push_back(-powerOf2 / d);
            }
        }
        return tInv;
    }

    static ZPolynomial addRaw(ZPolynomial a, ZPolynomial b) {
        auto zN = a.getZN();
        std::vector<BigFixedPoint> result(zN, BigFixedPoint::zero());
        for (size_t i = 0; i != zN; ++i) {
            result[i] = a[i] + b[i];
        }
        return result;
    }

    static ZPolynomial multiplyRaw(ZPolynomial a, ZPolynomial b) {
        auto zN = a.getZN();
        std::vector<BigFixedPoint> result(2 * zN - 1, BigFixedPoint::zero());
        for (size_t i = 0; i != zN; ++i) {
            for (size_t j = 0; j != zN; ++j) {
                result[i + j] += a[i] * b[j];
            }
        }
        auto two = BigFixedPoint::two();
        // now euclidean reduction mod X^zN - X + 2
        for (size_t i = result.size() - 1; i >= zN; --i) {
            result[i - zN + 1] += result[i];
            result[i - zN] += -two * result[i];
            result.pop_back();
        }
        assert(result.size() == zN);
        return result;
    }

    // Binary
    static ZPolynomial encodeBinary(uint32_t zN, BigInteger input) {
        auto one  = BigFixedPoint::one();
        auto zero = BigFixedPoint::zero();
        std::vector<BigFixedPoint> bits;
        // get bits of input in {0, 1}
        for (size_t i = 0; i != zN; ++i) {
            // OpenFHE have LSB = 1...
            auto flag = input.GetBitAtIndex(i + 1);
            //auto flag = (input & (1 << i)) >> i;
            if (flag) {
                bits.push_back(one);
            }
            else {
                bits.push_back(zero);
            }
        }
        return bits;
    }

    // Standard
    static ZPolynomial encode(uint32_t zN, BigInteger input) {
        return multiplyRaw(encodeBinary(zN, input), getTInv(zN));
    }

    static ZPolynomial encodeZeros(uint32_t zN) {
        return ZPolynomial(zN);
    }

    // Round to [-1, 1)
    static ZPolynomial roundNOneToOne(ZPolynomial input) {
        std::vector<BigFixedPoint> output;
        auto two = BigFixedPoint::positive(2);
        for (auto i : input.getCoefficients()) {
            auto j      = i / two;
            auto jRound = j.round();
            auto jFrac  = j - jRound;
            auto iFrac  = jFrac * two;
            output.push_back(iFrac);
        }
        return ZPolynomial(output);
    }

    // Round to (-1, 0]
    static ZPolynomial roundNOneToZero(ZPolynomial input) {
        // Add small epsilon to ensure range in (-1, 0)
        // epsilon = 2^{-zN-1}
        auto zN  = input.getZN();
        auto eps = BigFixedPoint::pow2(128 - zN - 1);
        // do [\cdot ]_1
        std::vector<BigFixedPoint> output;
        for (auto i : input.getCoefficients()) {
            i = i - eps;
            // then do ceil
            auto iCeil = i.ceil();
            output.push_back(i - iCeil + eps);
        }
        return output;
    }

    static ZPolynomial roundBigI(ZPolynomial input) {
        return roundNOneToZero(input);
    }

    static BigInteger decode(ZPolynomial input) {
        auto poly = multiplyRaw(input, getT(input.getZN()));

        // For each coefficient, do rounding
        std::vector<BigFixedPoint> result = poly.getCoefficients();
        for (size_t i = 0; i != result.size(); ++i) {
            auto rounded = result[i].round();
            result[i]    = rounded;
        }
        // Do Euclidean division by X-2
        auto two = BigFixedPoint::positive(2);
        for (size_t i = result.size() - 1; i >= 1; --i) {
            result[i - 1] += two * result[i];
            result.pop_back();
        }
        auto rounded = result[0].round();
        // Now do mod 2^zN
        auto modValueBFP = BigFixedPoint::pow2(input.getZN());
        rounded          = rounded - (rounded / modValueBFP).floor() * modValueBFP;
        //std::cout << "Decoded X-2: " << result[0].toHexString() << "\n";
        return rounded.getRoundedBigInteger();
    }

    static ZPolynomial add(ZPolynomial lhs, ZPolynomial rhs) {
        return roundBigI(addRaw(lhs, rhs));
    }

    static ZPolynomial multiply(ZPolynomial lhs, ZPolynomial rhs) {
        return roundBigI(multiplyRaw(multiplyRaw(lhs, rhs), getT(lhs.getZN())));
    }

    static ZPolynomial extractError(ZPolynomial input) {
        // Input: [m]_t / t + I + e
        // recoded: [m_t] / t
        auto zN      = input.getZN();
        auto decoded = decode(input);
        auto recoded = encode(zN, decoded);
        // error: I + e
        ZPolynomial error(zN);
        for (size_t i = 0; i != input.coefficients.size(); ++i) {
            error[i] = input[i] - recoded[i];
        }
        for (size_t i = 0; i != error.coefficients.size(); ++i) {
            auto errorRounded = error[i].round();
            error[i] -= errorRounded;
        }
        // Now we get pure e
        return error;
    }

    static ZPolynomial extractI(ZPolynomial input) {
        ZPolynomial Ipoly(input.getZN());
        // Input: [m]_t / t + I + e
        // recoded: [m_t] / t
        auto decoded = decode(input);
        auto recoded = encode(input.getZN(), decoded);
        // error: I + e
        for (size_t i = 0; i != input.coefficients.size(); ++i) {
            auto diff      = input[i] - recoded[i];
            auto diffRound = diff.round();
            Ipoly[i]       = diffRound;
        }
        return Ipoly;
    }

    static double getMaxLogI(ZPolynomial input) {
        auto I         = extractI(input);
        auto maxIValue = 1;  // for log2 calculation
        for (size_t i = 0; i != I.getCoefficients().size(); ++i) {
            auto Idouble = std::abs(I[i].round().convertToDouble());
            if (Idouble > maxIValue) {
                maxIValue = static_cast<int>(Idouble);
            }
        }
        return std::log2(maxIValue);
    }

    static ZPolynomial truncError(ZPolynomial input) {
        ZPolynomial output(input.getZN());
        auto error = extractError(input);
        for (size_t i = 0; i != error.coefficients.size(); ++i) {
            output[i] = input[i] - error[i];
        }
        return output;
    }

    static ZPolynomial extractStandardMessage(ZPolynomial input) {
        auto noError = truncError(input);
        auto decoded = decode(noError);
        auto recoded = encode(input.getZN(), decoded);
        return recoded;
    }

    CSlots toCSlots() const;

private:
    uint32_t zN;
    std::vector<BigFixedPoint> coefficients;
};

struct RPolynomial {
public:
    RPolynomial(const ZEncodingParams& p) : params(p), coefficients(p.getLength()) {}
    RPolynomial(const ZEncodingParams& p, const std::vector<BigFixedPoint>& coeffs) : params(p), coefficients(coeffs) {}
    std::vector<BigFixedPoint> getCoefficients() const {
        return coefficients;
    }

    BigFixedPoint& operator[](size_t index) {
        return coefficients[index];
    }
    const BigFixedPoint& operator[](size_t index) const {
        return coefficients[index];
    }

    ZEncodingParams getZEncodingParams() const {
        return params;
    }

    CSlots toCSlots() const;

private:
    ZEncodingParams params;
    std::vector<BigFixedPoint> coefficients;
};

struct CSlots {
public:
    CSlots(const ZEncodingParams& p) : params(p), slots(p.getLength() / 2) {}
    CSlots(const ZEncodingParams& p, const std::vector<BigComplex>& slotVec) : params(p), slots(slotVec) {}
    std::vector<BigComplex> getSlots() const {
        return slots;
    }

    BigComplex& operator[](size_t index) {
        return slots[index];
    }
    const BigComplex& operator[](size_t index) const {
        return slots[index];
    }

    uint32_t size() const {
        return slots.size();
    }

    ZEncodingParams getZEncodingParams() const {
        return params;
    }

    ZPolynomial getZPolynomial(size_t slotIndex) const;

    RPolynomial toRPolynomial() const;

    // These three functions requires sparse packing
    // Or if at full packing, concat cSlots of two ciphertexts
    std::pair<uint64_t, double> getIntegerAndErrorAtBooleanMode(size_t slotIndex) const;

    uint64_t getIntegerAtBooleanMode(size_t slotIndex) const;

    double getIntegerErrorAtBooleanMode(size_t slotIndex) const;

private:
private:
    ZEncodingParams params;
    std::vector<BigComplex> slots;
};

}  // namespace lbcrypto

#endif  // SRC_PKE_EXAMPLES_Z_ENCODE_UTILS_H_