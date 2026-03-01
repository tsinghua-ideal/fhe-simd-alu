#include "math/z-constants.h"

namespace lbcrypto {

std::unordered_map<uint32_t, ZLinearTransform::PrecomputedValues> ZLinearTransform::precomputedValues;

void ZLinearTransform::Initialize(uint32_t zN) {
#pragma omp critical
    if (precomputedValues.find(zN) == precomputedValues.end()) {
        precomputedValues.insert({zN, PrecomputedValues(zN)});
    }
}

ZLinearTransform::PrecomputedValues::PrecomputedValues(uint32_t zN) {
    m_zN    = zN;
    m_ZU    = GetZU(zN);
    m_ZUInv = GetZUInverse(zN);

    for (size_t i = 0; i != m_ZU.size(); ++i) {
        std::vector<std::complex<double>> row;
        for (size_t j = 0; j != m_ZU[0].size(); ++j) {
            row.push_back(
                std::complex<double>(m_ZU[i][j].getReal().convertToDouble(), m_ZU[i][j].getImag().convertToDouble()));
        }
        m_ZU_lowprec.push_back(row);
    }
    for (size_t i = 0; i != m_ZUInv.size(); ++i) {
        std::vector<std::complex<double>> row;
        for (size_t j = 0; j != m_ZUInv[0].size(); ++j) {
            row.push_back(std::complex<double>(m_ZUInv[i][j].getReal().convertToDouble(),
                                               m_ZUInv[i][j].getImag().convertToDouble()));
        }
        m_ZUInv_lowprec.push_back(row);
    }
}

BigCMatrix ZLinearTransform::PrecomputedValues::GetZU(uint32_t zN) {
    // Vandermond matrix
    BigCMatrix zu(zN / 2, std::vector<BigComplex>(zN));
    for (size_t i = 0; i != zN / 2; ++i) {
        zu[i][0] = BigFixedPoint::one();
        for (size_t j = 1; j != zN; ++j) {
            zu[i][j] = zu[i][j - 1] * Z_ROOTS_MAP.at(zN)[i];
        }
    }
    return zu;
}

BigCMatrix ZLinearTransform::PrecomputedValues::GetZUInverse(uint32_t zN) {
    std::vector<std::vector<BigComplex>> zUInv(zN, std::vector<BigComplex>(zN / 2));
    // Build the inverse Vandermond matrix
    for (size_t j = 0; j != zN / 2; ++j) {
        auto xi = Z_ROOTS_MAP.at(zN)[j];
        std::vector<BigComplex> xiPowers;
        auto one   = BigFixedPoint::one();
        auto zNbig = BigFixedPoint::positive(zN);
        xiPowers.push_back(one);
        for (size_t p = 1; p != zN; ++p) {
            xiPowers.push_back(xiPowers[p - 1] * xi);
        }
        // build power map
        for (size_t i = 0; i != zN; ++i) {
            // d = zN * xi^(zN-1) - 1
            auto d = (zNbig * xiPowers[zN - 1]) - one;
            if (i == 0) {
                zUInv[i][j] = (xiPowers[zN - 1] - one) / d;
            }
            else {
                zUInv[i][j] = xiPowers[zN - 1 - i] / d;
            }
        }
    }
    return zUInv;
}

const BigCMatrix& ZLinearTransform::GetZU(uint32_t zN) {
    auto it = precomputedValues.find(zN);
    if (it == precomputedValues.end()) {
        OPENFHE_THROW("ZLinearTransform::Initialize not called for " + std::to_string(zN));
    }
    return it->second.m_ZU;
}

const BigCMatrix& ZLinearTransform::GetZUInverse(uint32_t zN) {
    auto it = precomputedValues.find(zN);
    if (it == precomputedValues.end()) {
        OPENFHE_THROW("ZLinearTransform::Initialize not called for " + std::to_string(zN));
    }
    return it->second.m_ZUInv;
}

#define LOW_PREC
#ifndef LOW_PREC
BigCVector ZLinearTransform::MultZU(uint32_t zN, const BigFPVector& input) {
    auto& U = GetZU(zN);
    assert(input.size() == U[0].size() && "Input size does not match expected size for U");

    // This is the vandermond matrix
    std::vector<BigComplex> result;

    size_t halfSize = input.size() / 2;
    for (size_t i = 0; i < halfSize; ++i) {
        BigComplex sum;
        for (size_t j = 0; j < input.size(); ++j) {
            sum = sum + U[i][j] * input[j];
        }
        result.push_back(sum);
    }
    return result;
}

BigFPVector ZLinearTransform::MultZUInverse(uint32_t zN, const BigCVector& input) {
    auto& UInv = GetZUInverse(zN);
    std::vector<BigFixedPoint> result;
    assert(input.size() == UInv[0].size() && "Input size does not match expected size for multiplyByZUInverse");

    for (size_t i = 0; i != zN; ++i) {
        BigComplex sum;
        for (size_t j = 0; j != zN / 2; ++j) {
            sum = sum + UInv[i][j] * input[j];
        }
        // z + conj(z) = 2*real(z)
        auto two = BigFixedPoint::positive(2);
        result.push_back(two * sum.getReal());
    }
    return result;
}
#else

const std::vector<std::vector<std::complex<double>>>& ZLinearTransform::GetZULowPrec(uint32_t zN) {
    auto it = precomputedValues.find(zN);
    if (it == precomputedValues.end()) {
        OPENFHE_THROW("ZLinearTransform::Initialize not called for " + std::to_string(zN));
    }
    return it->second.m_ZU_lowprec;
}

const std::vector<std::vector<std::complex<double>>>& ZLinearTransform::GetZUInverseLowPrec(uint32_t zN) {
    auto it = precomputedValues.find(zN);
    if (it == precomputedValues.end()) {
        OPENFHE_THROW("ZLinearTransform::Initialize not called for " + std::to_string(zN));
    }
    return it->second.m_ZUInv_lowprec;
}

BigCVector ZLinearTransform::MultZU(uint32_t zN, const BigFPVector& input) {
    auto& U = GetZULowPrec(zN);
    assert(input.size() == U[0].size() && "Input size does not match expected size for U");

    // This is the vandermond matrix
    std::vector<BigComplex> result;

    std::vector<double> inputDouble(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        inputDouble[i] = input[i].convertToDouble();
    }

    size_t halfSize = input.size() / 2;
    for (size_t i = 0; i < halfSize; ++i) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j < input.size(); ++j) {
            sum += U[i][j] * inputDouble[j];
        }
        result.push_back(BigComplex(BigFixedPoint::fromDouble(sum.real()), BigFixedPoint::fromDouble(sum.imag())));
    }
    return result;
}

BigFPVector ZLinearTransform::MultZUInverse(uint32_t zN, const BigCVector& input) {
    auto& UInv = GetZUInverseLowPrec(zN);
    std::vector<BigFixedPoint> result;
    assert(input.size() == UInv[0].size() && "Input size does not match expected size for multiplyByZUInverse");

    std::vector<std::complex<double>> inputDouble(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        inputDouble[i] =
            std::complex<double>(input[i].getReal().convertToDouble(), input[i].getImag().convertToDouble());
    }

    for (size_t i = 0; i != zN; ++i) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t j = 0; j != zN / 2; ++j) {
            sum += UInv[i][j] * inputDouble[j];
        }
        // z + conj(z) = 2*real(z)
        auto two = BigFixedPoint::positive(2);
        result.push_back(two * BigFixedPoint::fromDouble(sum.real()));
    }
    return result;
}

#endif

}  // namespace lbcrypto