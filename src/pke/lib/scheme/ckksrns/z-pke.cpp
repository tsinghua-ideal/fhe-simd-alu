#include "scheme/ckksrns/z-pke.h"
#include "encoding/z-encoding.h"

namespace lbcrypto {

PKEZ pkeZ_global;  // Temporary

//=============================================================================
// Encryption Utils
//=============================================================================

// DecryptCore not accessible from CryptoContext
// so copy from @openfhe//src/pke/lib/schemerns/rns-pke.cpp
static DCRTPoly DecryptCore(const std::vector<DCRTPoly>& cv, const PrivateKey<DCRTPoly> privateKey) {
    const DCRTPoly& s = privateKey->GetPrivateElement();

    size_t sizeQ  = s.GetParams()->GetParams().size();
    size_t sizeQl = cv[0].GetParams()->GetParams().size();

    size_t diffQl = sizeQ - sizeQl;

    auto scopy(s);
    scopy.DropLastElements(diffQl);

    DCRTPoly sPower(scopy);

    DCRTPoly b(cv[0]);
    b.SetFormat(Format::EVALUATION);

    DCRTPoly ci;
    for (size_t i = 1; i < cv.size(); i++) {
        ci = cv[i];
        ci.SetFormat(Format::EVALUATION);

        b += sPower * ci;
        sPower *= scopy;
    }
    return b;
}

static std::shared_ptr<std::vector<DCRTPoly>> EncryptZeroCore(const PublicKey<DCRTPoly> publicKey) {
    const auto cryptoParams =
        std::dynamic_pointer_cast<CryptoParametersRLWE<DCRTPoly>>(publicKey->GetCryptoParameters());

    const auto ns = cryptoParams->GetNoiseScale();

    auto elementParams = cryptoParams->GetElementParams();

    auto pk = publicKey->GetPublicElements();

    auto p0 = pk[0];
    auto p1 = pk[1];

    uint32_t sizeQ  = elementParams->GetParams().size();
    uint32_t sizePK = p0.GetParams()->GetParams().size();

    if (sizePK > sizeQ) {
        p0.DropLastElements(sizePK - sizeQ);
        p1.DropLastElements(sizePK - sizeQ);
    }

    DCRTPoly::TugType tug;
    DCRTPoly v = DCRTPoly(tug, elementParams, Format::EVALUATION);

    // noise generation with the discrete gaussian generator dgg
    auto& dgg = cryptoParams->GetDiscreteGaussianGenerator();
    DCRTPoly e0(dgg, elementParams, Format::EVALUATION);
    DCRTPoly e1(dgg, elementParams, Format::EVALUATION);

    DCRTPoly b(elementParams);
    DCRTPoly a(elementParams);

    b = p0 * v + ns * e0;
    a = p1 * v + ns * e1;

    return std::make_shared<std::vector<DCRTPoly>>(std::initializer_list<DCRTPoly>({std::move(b), std::move(a)}));
}

Ciphertext<DCRTPoly> PKEZImpl::Encrypt(ZEncoding ptxt) {
    if (!pk) {
        OPENFHE_THROW("Public key is not set for encryption");
    }
    auto zEncDCRTPoly = ptxt->GetElement<DCRTPoly>();

    auto ba = EncryptZeroCore(pk);
    (*ba)[0] += zEncDCRTPoly;

    auto ctxt = std::make_shared<CiphertextImpl<DCRTPoly>>(pk);
    ctxt->SetElements(std::move(*ba));
    ctxt->SetNoiseScaleDeg(1);
    ctxt->SetScalingFactorBFP(ptxt->GetScalingFactorBFP());
    ctxt->SetZEncodingParams(ptxt->GetZEncodingParams());
    return ctxt;
}

CiphertextGroup PKEZImpl::Encrypt(std::vector<ZEncoding> ptxts) {
    std::vector<Ciphertext<DCRTPoly>> ctxts;
    ctxts.reserve(ptxts.size());
    for (const auto& ptxt : ptxts) {
        ctxts.push_back(Encrypt(ptxt));
    }
    return CiphertextGroup(ctxts);
}

std::vector<uint64_t> PKEZImpl::DecryptSmallFast(CiphertextGroup cts) {
    if (!sk) {
        OPENFHE_THROW("Secret key is not set for decryption");
    }
    auto ct            = cts[0];
    auto b             = DecryptCore(ct->GetElements(), sk);
    auto sfBigFP       = ct->GetScalingFactorBFP();
    auto zEncodeParams = ct->GetZEncodingParams();
    auto sfBFP         = ct->GetScalingFactorBFP();
    auto zEncode       = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, zEncodeParams);

    if (zEncodeParams.isZMode()) {
        return ZEncodingImpl::decodeArithSmall(zEncode);
    }
    else if (zEncodeParams.isBModeFull()) {
        OPENFHE_THROW("Unknown encoding mode");
    }
    else {
        OPENFHE_THROW("Unknown encoding mode");
    }
}

ZDecryptResult PKEZImpl::Decrypt(CiphertextGroup cts) {
    if (!sk) {
        OPENFHE_THROW("Secret key is not set for decryption");
    }
    auto ct            = cts[0];
    auto b             = DecryptCore(ct->GetElements(), sk);
    auto sfBigFP       = ct->GetScalingFactorBFP();
    auto zEncodeParams = ct->GetZEncodingParams();
    auto zN            = zEncodeParams.getZN();
    auto zSlots        = zEncodeParams.getZSlots();
    auto level         = ct->GetLevel();
    auto sfBFP         = ct->GetScalingFactorBFP();
    auto encodingType  = zEncodeParams.getEncodingType();
    auto zEncode       = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, zEncodeParams);
    RPolynomial rPoly  = ZEncodingImpl::decodeR(zEncode);

    std::vector<BigInteger> values(zSlots);
    double log2MaxNoise;
    double log2MaxI;

    if (zEncodeParams.isZMode()) {
        auto cSlots = rPoly.toCSlots();
        std::vector<ZPolynomial> zPolys(zSlots, ZPolynomial(zN));
        std::vector<double> log2Noises(zSlots);
        std::vector<double> log2I(zSlots);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(zSlots))
        for (size_t i = 0; i != zSlots; ++i) {
            zPolys[i]     = cSlots.getZPolynomial(i);
            values[i]     = ZPolynomial::decode(zPolys[i]);
            log2Noises[i] = ZPolynomial::extractError(zPolys[i]).getLog2Norm();
            log2I[i]      = ZPolynomial::getMaxLogI(zPolys[i]);
        }
        log2MaxNoise = *std::max_element(log2Noises.begin(), log2Noises.end());
        log2MaxI     = *std::max_element(log2I.begin(), log2I.end());
    }
    else if (zEncodeParams.isBModeFull()) {
        auto b2            = DecryptCore(cts[1]->GetElements(), sk);
        auto zEncode2      = std::make_shared<ZEncodingImpl>(b2.GetParams(), b2, sfBigFP, zEncodeParams);
        RPolynomial rPoly2 = ZEncodingImpl::decodeR(zEncode2);
        auto cSlots        = rPoly.toCSlots();
        auto cSlots2       = rPoly2.toCSlots();

        std::vector<double> log2Noises(zSlots);

#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(zSlots))
        for (size_t i = 0; i != zSlots; ++i) {
            BigInteger value = 0;
            std::vector<double> noises(zN);
            for (size_t j = 0; j != zN / 2; ++j) {
                auto lowJ          = cSlots[i * (zN / 2) + j].getReal();
                auto highJ         = cSlots2[i * (zN / 2) + j].getReal();
                auto lowJInt       = lowJ.round();
                auto highJInt      = highJ.round();
                auto lowJBigInt    = lowJInt.getRoundedBigInteger();
                auto highJBigInt   = highJInt.getRoundedBigInteger();
                auto lowJNoise     = (lowJ - lowJInt).log2Norm();
                auto highJNoise    = (highJ - highJInt).log2Norm();
                noises[j]          = lowJNoise;
                noises[zN / 2 + j] = highJNoise;
                // reconstruct value
                value += lowJBigInt << j;
                value += highJBigInt << (j + zN / 2);
            }
            values[i]     = value;
            log2Noises[i] = *std::max_element(noises.begin(), noises.end());
        }
        log2MaxNoise = *std::max_element(log2Noises.begin(), log2Noises.end());
        log2MaxI     = 0;
    }
    else {
        OPENFHE_THROW("Unknown encoding mode");
    }

    return ZDecryptResult(values, log2MaxNoise, log2MaxI, level, sfBFP, encodingType);
}

static std::string bigIntegerToHexString(BigInteger value) {
    std::stringstream ss;
    NTL::ZZ v = value;

    if (v == 0) {
        return "0x00";
    }

    long nbytes = NumBytes(v);
    std::vector<unsigned char> buf(nbytes);
    BytesFromZZ(buf.data(), v, nbytes);
    std::reverse(buf.begin(), buf.end());

    std::ostringstream oss;
    oss << "0x";

    for (unsigned char b : buf) {
        oss << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(b);
    }
    return oss.str();
};

void ZDecryptResult::print(std::string msg, size_t maxSlotsToPrint) const {
    std::cout << msg << " DecryptResult: " << std::endl;
    for (size_t i = 0; i != std::min(values.size(), maxSlotsToPrint); ++i) {
        std::cout << "  Value[" << i << "]: " << bigIntegerToHexString(values[i]) << std::endl;
    }
    std::cout << "  ... (total " << values.size() << " slots)" << std::endl;
    std::cout << "  logMaxNoise: " << logMaxNoise << std::endl;
    if (logMaxI != 0)
        std::cout << "  logMaxI: " << logMaxI << std::endl;
    std::cout << "  level: " << level << std::endl;
    std::cout << "  Scaling factor log2: " << sfBFP.log2Norm() << std::endl;
    std::cout << "  Encoding type: " << encodingType << std::endl;
}

void PKEZImpl::debug(Ciphertext<DCRTPoly> ct, std::string msg) {
    auto b = DecryptCore(ct->GetElements(), sk);

    auto sfBigFP = ct->GetScalingFactorBFP();
    auto log2sf  = std::log2(sfBigFP.convertToDouble());
    std::cout << msg << "  Scaling factor log2: " << std::setprecision(20) << log2sf << std::endl;
    auto level = ct->GetLevel();
    std::cout << msg << "  level: " << level << std::endl;
    auto zEncodeParams = ct->GetZEncodingParams();
    auto zN            = zEncodeParams.getZN();
    auto zSlots        = zEncodeParams.getZSlots();
    auto zEncode       = std::make_shared<ZEncodingImpl>(b.GetParams(), b, sfBigFP, zEncodeParams);

    RPolynomial values = ZEncodingImpl::decodeR(zEncode);

    enum class DecodeMode { RDecode, CSlotsDecode, ZDecode };
    std::map<std::string, DecodeMode> decodeMap = {
        // RDecode
        {"Z2R", DecodeMode::RDecode},
        // CDecode
        {"Scaled", DecodeMode::CSlotsDecode},
        {"Z2C0", DecodeMode::CSlotsDecode},
        {"R2C0", DecodeMode::CSlotsDecode},
    };

    auto decodeModeIt = decodeMap.find(msg);
    if (decodeModeIt == decodeMap.end()) {
        return;
    }
    auto decodeMode = decodeModeIt->second;

    if (decodeMode == DecodeMode::RDecode) {
        // Threshold print to avoid too much output
        for (size_t i = 0; i != std::min(values.getCoefficients().size(), 8ul); ++i) {
            auto valueI = values[i];
            std::cout << msg << "  values [" << i << "]: " << values[i].toHexString(ceil(log2sf / 4.0)) << std::endl;
        }
        if (values.getCoefficients().size() > 4) {
            std::cout << msg << "  ... (total " << values.getCoefficients().size() << " coefficients)" << std::endl;
        }
    }
    if (decodeMode == DecodeMode::CSlotsDecode) {
        auto cSlots          = values.toCSlots();
        auto maxSlotsToPrint = std::min(zN / 2, uint(16));
        for (size_t i = 0; i != maxSlotsToPrint; ++i) {
            auto value    = cSlots[i].getReal();
            auto p        = BigFixedPoint::positive(16);
            auto lutPart  = (value * p).round() / p;
            auto fracPart = value - lutPart;
            //if (i % 4 == 3) {
            std::cout << msg << "  cSlots Slot " << i << " " << cSlots[i].getReal().toHexString(ceil(log2sf / 4.0))
                      << " (LUT part: " << lutPart.toHexString(ceil(log2sf / 4.0))
                      << ", frac part: " << fracPart.toHexString(ceil(log2sf / 4.0)) << ")" << std::endl;
        }
        if (zSlots > maxSlotsToPrint) {
            std::cout << msg << "  ... (total " << zSlots << " slots)" << std::endl;
        }
    }
    //if (zEncodeParams.isZMode()) {
    //    auto cSlots = values.toCSlots();
    //    std::vector<ZPolynomial> zPolys(zSlots, ZPolynomial(zN));
    //    auto maxSlotsToPrint = std::min(zSlots, uint(8));
    //#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(maxSlotsToPrint))
    //    for (size_t i = 0; i != maxSlotsToPrint; ++i) {
    //        zPolys[i] = cSlots.getZPolynomial(i);
    //    }
    //    for (size_t i = 0; i != maxSlotsToPrint; ++i) {
    //        //printZPoly(zPolys[i], msg);
    //    }
    //    if (zSlots > maxSlotsToPrint) {
    //        std::cout << msg << "  ... (total " << zSlots << " slots)" << std::endl;
    //    }
    //}
    return;
}

bool ZDecryptResult::valuesEqual(const ZDecryptResult& other) const {
    if (values.size() != other.values.size()) {
        return false;
    }
    for (size_t i = 0; i != values.size(); ++i) {
        if (values[i] != other.values[i]) {
            std::cout << "Value mismatch at index " << i << ": " << bigIntegerToHexString(values[i])
                      << " != " << bigIntegerToHexString(other.values[i]) << std::endl;
            return false;
        }
    }
    return true;
}

}  // namespace lbcrypto