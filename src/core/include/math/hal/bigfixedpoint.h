#ifndef SRC_CORE_INCLUDE_MATH_HAL_BIGFIXEDPOINT_H_
#define SRC_CORE_INCLUDE_MATH_HAL_BIGFIXEDPOINT_H_

#include <cassert>
#include "math/math-hal.h"

namespace lbcrypto {

struct BigFixedPoint {
public:
    BigFixedPoint() : value(0), log2Scale(0), neg(false) {}
    BigFixedPoint(const BigInteger& val, int log2S, bool neg) : value(val), log2Scale(log2S), neg(neg) {}

    // Should not use.
    //explicit BigFixedPoint(int a, int log2S = 128) {
    //    log2Scale = 0;
    //    if (a < 0) {
    //        value = -a;
    //        neg   = true;
    //    }
    //    else {
    //        value = a;
    //        neg   = false;
    //    }
    //}

    static BigFixedPoint zero(int log2S = 128) {
        return BigFixedPoint(BigInteger(0), log2S, false);
    }
    static BigFixedPoint one(int log2S = 128) {
        return BigFixedPoint(BigInteger(1) << log2S, log2S, false);
    }
    static BigFixedPoint two(int log2S = 128) {
        return BigFixedPoint(BigInteger(2) << log2S, log2S, false);
    }
    // Requires power > -128
    static BigFixedPoint pow2(int power, int log2S = 128) {
        return BigFixedPoint(BigInteger(1) << (log2S + power), log2S, false);
    }
    static BigFixedPoint half(int log2S = 128) {
        return BigFixedPoint(BigInteger(1) << (log2S - 1), log2S, false);
    }
    static BigFixedPoint positive(uint64_t a, int log2S = 128) {
        return BigFixedPoint(BigInteger(a) << log2S, log2S, false);
    }
    static BigFixedPoint fromDouble(double a, int log2S = 128) {
        auto neg = false;
        if (a < 0) {
            neg = true;
            a   = -a;
        }
        // find log2 scale to represent fractional part
        auto log2Scale = 0;
        double intpart;
        double fracpart = std::modf(a, &intpart);
        while (fracpart != 0.0 && log2Scale < 128) {
            a *= 2.0;
            fracpart = std::modf(a, &intpart);
            log2Scale++;
        }
        auto value = BigInteger(static_cast<int64_t>(std::round(a)));
        return BigFixedPoint(value, log2Scale, neg).scaleTo(log2S);
    }

    BigInteger getValue() const {
        return value;
    }
    int getLog2Scale() const {
        return log2Scale;
    }
    bool getNeg() const {
        return neg;
    }

    BigFixedPoint scaleUp(int log2S) const {
        return BigFixedPoint(getValue() << log2S, getLog2Scale() + log2S, getNeg());
    }
    BigFixedPoint scaleDown(int log2S) const {
        return BigFixedPoint(getValue() >> log2S, getLog2Scale() - log2S, getNeg());
    }
    BigFixedPoint scaleTo(int log2S) const {
        if (log2S > log2Scale) {
            return scaleUp(log2S - log2Scale);
        }
        if (log2S < log2Scale) {
            return scaleDown(log2Scale - log2S);
        }
        return *this;
    }

    bool equalZero() const {
        return (value == BigInteger(0));
    }

    // -40 is based on the consideration we often choose Delta = 2^40 (?TODO)
    bool almostEqual(const BigFixedPoint& other, int maxLog2Diff = -40) const {
        auto a = this->scaleTo(std::max(this->getLog2Scale(), other.getLog2Scale()));
        auto b = other.scaleTo(std::max(this->getLog2Scale(), other.getLog2Scale()));
        if (a.getNeg() != b.getNeg()) {
            return false;
        }
        auto scale  = a.getLog2Scale();
        auto aValue = a.getValue();
        auto bValue = b.getValue();
        auto diff   = aValue > bValue ? aValue - bValue : bValue - aValue;
        return diff <= (BigInteger(1) << (scale + maxLog2Diff));
    }

    BigFixedPoint floor() const {
        auto intPart  = getValue() >> getLog2Scale();
        auto intBFP   = BigFixedPoint(intPart, 0, getNeg()).scaleTo(log2Scale);
        auto fracPart = *this - intBFP;
        if (equalZero()) {
            return BigFixedPoint(0, 0, false).scaleTo(log2Scale);
        }
        if (getNeg()) {
            if (fracPart.equalZero()) {
                return BigFixedPoint(intPart, 0, getNeg()).scaleTo(log2Scale);
            }
            else {
                return BigFixedPoint(intPart.Add(1), 0, getNeg()).scaleTo(log2Scale);
            }
        }
        else {
            return BigFixedPoint(intPart, 0, getNeg()).scaleTo(log2Scale);
        }
    }

    BigFixedPoint ceil() const {
        auto intPart  = getValue() >> getLog2Scale();
        auto intBFP   = BigFixedPoint(intPart, 0, getNeg()).scaleTo(log2Scale);
        auto fracPart = *this - intBFP;
        if (equalZero()) {
            return BigFixedPoint(0, 0, false).scaleTo(log2Scale);
        }
        if (getNeg()) {
            return BigFixedPoint(intPart, 0, getNeg()).scaleTo(log2Scale);
        }
        else {
            if (fracPart.equalZero()) {
                return BigFixedPoint(intPart, 0, getNeg()).scaleTo(log2Scale);
            }
            else {
                return BigFixedPoint(intPart.Add(1), 0, getNeg()).scaleTo(log2Scale);
            }
        }
    }

    BigFixedPoint round() const {
        auto one  = BigFixedPoint(1, 0, false).scaleTo(log2Scale);
        auto half = BigFixedPoint(1, 1, false).scaleTo(log2Scale);
        // first do mod 1
        auto i      = *this;
        auto iFloor = i.floor();
        // In [0, 1)
        auto iFracPos = i - iFloor;
        // then do round
        BigFixedPoint frac;
        auto diff = iFracPos - half;
        if (diff.equalZero() || !diff.getNeg()) {
            frac = iFracPos - one;
        }
        else {
            frac = iFracPos;
        }
        return i - frac;
    }

    uint64_t getRoundedInteger() const {
        auto rounded = this->round();
        return (rounded.getValue() >> rounded.getLog2Scale()).ConvertToInt();
    }

    BigInteger getRoundedBigInteger() const {
        auto rounded = this->round();
        return rounded.getValue() >> rounded.getLog2Scale();
    }

    double log2Norm() const {
        auto valDouble = value.ConvertToDouble();
        auto norm      = std::log2(valDouble) - static_cast<double>(log2Scale);
        return norm;
    }

    double convertToDouble() const {
        double val = value.ConvertToDouble();
        val /= std::pow(2.0, log2Scale);
        if (neg) {
            val = -val;
        }
        return val;
    }

    std::string toString(int prec = 32) const {
        std::string result;

        if (neg && value.Compare(0) != 0) {
            result += "-";
        }

        // Calculate the integer and fractional parts
        BigInteger scaledValue    = value;
        BigInteger integerPart    = scaledValue >> log2Scale;                  // Divide by 2^log2Scale
        BigInteger fractionalPart = scaledValue - (integerPart << log2Scale);  // Remainder

        // Convert integer part to string
        result += integerPart.ToString();

        // Convert fractional part to decimal
        if (fractionalPart.Compare(0) != 0 && log2Scale > 0) {
            result += ".";

            // Convert fractional part (which is in base-2) to base-10
            BigInteger frac            = fractionalPart;
            const int maxDecimalDigits = prec;  // Maximum decimal digits to display

            for (int i = 0; i < maxDecimalDigits && frac.Compare(0) != 0; i++) {
                frac.MulEq(10);                        // Multiply by 10 to get next decimal digit
                BigInteger digit = frac >> log2Scale;  // Extract the digit
                result += digit.ToString();
                frac = frac - (digit << log2Scale);  // Remainder for next iteration

                // Early exit if remainder becomes zero
                if (frac.Compare(0) == 0) {
                    break;
                }
            }
        }
        else if (log2Scale > 0) {
            // Add .0 if there's fractional precision but value is integer
            result += ".0";
        }

        return result;
    }

    std::string toHexString(int prec = 32) const {
        std::string result;

        if (neg && value.Compare(0) != 0) {
            result += "-";
        }

        // Calculate the integer and fractional parts
        BigInteger scaledValue    = value;
        BigInteger integerPart    = scaledValue >> log2Scale;                  // Divide by 2^log2Scale
        BigInteger fractionalPart = scaledValue - (integerPart << log2Scale);  // Remainder

        // Convert integer part to string
        if (integerPart.Compare(0) != 0) {
            // Manual binary conversion for integer part
            BigInteger temp = integerPart;
            std::string hexStr;
            while (temp.Compare(0) != 0) {
                // OpenFHE count from 1???
                auto rem = temp.Mod(16);
                std::stringstream ss;
                ss << std::hex << static_cast<int64_t>(rem.ConvertToInt());
                hexStr = ss.str() + hexStr;
                temp >>= 4;
            }
            result += hexStr;
        }
        else {
            result += "0";
        }

        // Convert fractional part to decimal
        if (fractionalPart.Compare(0) != 0 && log2Scale > 0) {
            result += ".";

            // Convert fractional part (which is in base-2) to base-10
            BigInteger frac            = fractionalPart;
            const int maxDecimalDigits = prec;  // Maximum decimal digits to display

            std::string fracStr;

            for (int i = 0; i < maxDecimalDigits && frac.Compare(0) != 0; i++) {
                frac.MulEq(16);                        // Multiply by 10 to get next decimal digit
                BigInteger digit = frac >> log2Scale;  // Extract the digit
                std::stringstream ss;
                ss << std::hex << static_cast<int64_t>(digit.ConvertToInt());
                fracStr += ss.str();
                frac = frac - (digit << log2Scale);  // Remainder for next iteration

                // Early exit if remainder becomes zero
                if (frac.Compare(0) == 0) {
                    break;
                }
                // If there is still more digit
                // And we have 8 char (32 bits)
                // Add a '_'
                if ((i + 1) % 8 == 0 && (i + 1) < maxDecimalDigits && frac.Compare(0) != 0) {
                    fracStr += "_";
                }
            }
            result += fracStr;
        }
        else if (log2Scale > 0) {
            // Add .0 if there's fractional precision but value is integer
            result += ".0";
        }

        return result;
    }

    std::string toBinary(int prec = 64) const {
        std::string result;

        if (neg && value.Compare(0) != 0) {
            result += "-";
        }

        if (value.Compare(0) == 0) {
            return "0";
        }

        // Convert integer part to binary manually
        BigInteger integerPart = value >> log2Scale;
        if (integerPart.Compare(0) != 0) {
            // Manual binary conversion for integer part
            BigInteger temp = integerPart;
            std::string binaryStr;
            while (temp.Compare(0) != 0) {
                // OpenFHE count from 1???
                if (temp.GetBitAtIndex(1)) {
                    binaryStr = "1" + binaryStr;
                }
                else {
                    binaryStr = "0" + binaryStr;
                }
                temp >>= 1;
            }
            result += binaryStr;
        }
        else {
            result += "0";
        }

        // Convert fractional part
        if (log2Scale > 0) {
            BigInteger fractionalPart = value - (integerPart << log2Scale);

            if (fractionalPart.Compare(0) != 0) {
                result += ".";

                BigInteger frac = fractionalPart;
                for (int i = 0; i < log2Scale && i < prec && frac.Compare(0) != 0; i++) {
                    frac <<= 1;  // Shift left to check next bit
                    if (frac >= (BigInteger(1) << log2Scale)) {
                        result += "1";
                        frac = frac - (BigInteger(1) << log2Scale);
                    }
                    else {
                        result += "0";
                    }

                    if (frac.Compare(0) == 0) {
                        break;
                    }
                }
            }
            else {
                result += ".0";
            }
        }
        return result;
    }

    friend BigFixedPoint operator+(const BigFixedPoint& a, const BigFixedPoint& b);

    BigFixedPoint& operator+=(const BigFixedPoint& other) {
        *this = *this + other;
        return *this;
    }

    friend BigFixedPoint operator-(const BigFixedPoint& a, const BigFixedPoint& b);

    friend BigFixedPoint operator-(const BigFixedPoint& a);
    BigFixedPoint& operator-=(const BigFixedPoint& other) {
        *this = *this - other;
        return *this;
    }

    template <class Archive>
    void serialize(Archive& ar) {
        ar(CEREAL_NVP(value), CEREAL_NVP(log2Scale), CEREAL_NVP(neg));
    }

private:
    BigInteger value;  // unsigned big
    int log2Scale;     // scale = 2^log2Scale
    bool neg;
    // we actually encode neg * value * 2^{-log2Scale}
};

BigFixedPoint operator+(const BigFixedPoint& a, const BigFixedPoint& b);

BigFixedPoint operator-(const BigFixedPoint& a, const BigFixedPoint& b);

BigFixedPoint operator-(const BigFixedPoint& a);

// make it scale preserving
BigFixedPoint operator*(const BigFixedPoint& a, const BigFixedPoint& b);

BigFixedPoint operator/(const BigFixedPoint& a, const BigFixedPoint& b);

bool operator<(const BigFixedPoint& a, const BigFixedPoint& b);

struct BigComplex {
public:
    BigComplex() : real(), imag() {}
    BigComplex(const BigFixedPoint& re, const BigFixedPoint& im) : real(re), imag(im) {}
    BigComplex(const BigFixedPoint& re) : real(re), imag(BigFixedPoint(0, 0, false).scaleTo(re.getLog2Scale())) {}
    BigFixedPoint getReal() const {
        return real;
    }
    BigFixedPoint getImag() const {
        return imag;
    }

    BigComplex scaleTo(int log2S) const {
        return BigComplex(real.scaleTo(log2S), imag.scaleTo(log2S));
    }
    BigComplex conj() const {
        return BigComplex(real, -imag);
    }

    double log2Norm() const {
        auto length = real * real + imag * imag;
        auto norm   = 0.5 * length.log2Norm();
        return norm;
    }

    std::complex<double> convertToComplex() const {
        return {real.convertToDouble(), imag.convertToDouble()};
    }

    std::string toString(int prec = 32) const {
        return "(" + real.toString(prec) + ", " + imag.toString(prec) + ")";
    }

    std::string toHexString(int prec = 32) const {
        return "(" + real.toHexString(prec) + ", " + imag.toHexString(prec) + ")";
    }

    std::string toBinary(int prec = 64) const {
        return "(" + real.toBinary(prec) + ", " + imag.toBinary(prec) + ")";
    }

    BigComplex& operator+=(const BigComplex& a);
    BigComplex& operator-=(const BigComplex& a);
    BigComplex& operator*=(const BigComplex& a);
    BigComplex& operator/=(const BigComplex& a);

    bool equalZero() const {
        return real.equalZero() && imag.equalZero();
    }

private:
    BigFixedPoint real;
    BigFixedPoint imag;
};

BigComplex operator+(const BigComplex& a, const BigComplex& b);
BigComplex operator-(const BigComplex& a, const BigComplex& b);
BigComplex operator*(const BigComplex& a, const BigComplex& b);
BigComplex operator/(const BigComplex& a, const BigComplex& b);

BigComplex operator-(const BigComplex& a);

struct BigComplexCompare {
    bool operator()(const BigComplex& a, const BigComplex& b) const {
        if (a.getReal().almostEqual(b.getReal())) {
            return a.getImag() < b.getImag();
        }
        return a.getReal() < b.getReal();
    }
};

using BigCVector = std::vector<BigComplex>;
using BigCMatrix = std::vector<BigCVector>;

using BigFPVector = std::vector<BigFixedPoint>;

BigCVector ToCVector(const BigFPVector& input);
BigFPVector ToReal(const BigCVector& input);

std::vector<std::complex<double>> ToStdComplexVector(const BigCVector& input);
BigCVector FromStdComplexVector(const std::vector<std::complex<double>>& input);

}  // namespace lbcrypto

#endif  // SRC_CORE_INCLUDE_MATH_HAL_BIGFIXEDPOINT_H_