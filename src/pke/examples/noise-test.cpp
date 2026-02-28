#include "openfhe.h"
#include "scheme/ckksrns/z-user.h"
#include "utils.h"
#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/z-pke.h"

using namespace lbcrypto;

void NoiseTestExample(int zN, std::string directive);

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <zN> <test|verify|bench>" << std::endl;
    }
    auto zN = 64;
    if (argc >= 2) {
        zN = std::stoi(argv[1]);
    }
    auto directive = "test";
    ;
    if (argc >= 3) {
        directive = argv[2];
    }

    NoiseTestExample(zN, directive);
}

void NoiseTestExample(int zN, std::string directive) {
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::SPARSE_ENCAPSULATED);

    uint32_t repeats = 1;
    if (directive == "bench" || directive == "verify") {
        parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetRingDim(1 << 16);
        if (directive == "bench")
            repeats = 20;
    }
    else {  // test
        parameters.SetSecurityLevel(lbcrypto::HEStd_NotSet);
        parameters.SetRingDim(1 << 10);
    }

    ScalingTechnique rescaleTech = FLEXIBLEMANUAL;
    uint32_t dcrtBits            = 43;
    // bit size for aux moduli in P
    // for HEXL acceleration. Extra 6 bit for SPARSE_ENCAPSULATED
    AUXMODSIZE_FLEXIBLEMANUAL = 50;
    std::cout << "Scaling Factor: " << dcrtBits << " bits\n";
    std::cout << "Auxiliary Prime Size: " << AUXMODSIZE_FLEXIBLEMANUAL << " bits\n";

    parameters.SetScalingModSize(dcrtBits);
    parameters.SetFirstModSize(dcrtBits);
    parameters.SetScalingTechnique(rescaleTech);
    parameters.SetNumLargeDigits(3);

    std::vector<uint32_t> levelBudget = {3, 2};
    std::cout << "Level Budget: ";
    for (auto lb : levelBudget) {
        std::cout << lb << " ";
    }
    std::cout << "\n";

    uint32_t mulDepth = 20;
    parameters.SetMultiplicativeDepth(mulDepth);
    std::cout << "Multiplicative Depth: " << mulDepth << "\n";

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);

    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);

    uint32_t ringDim = cc->GetRingDimension();
    std::cout << "CKKS scheme ring dimension: " << ringDim << "\n";

    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);

    const auto cryptoParams = std::dynamic_pointer_cast<CryptoParametersCKKSRNS>(cc->GetCryptoParameters());
    //std::cout << *cryptoParams << std::endl;
    //std::cout << *(cryptoParams->GetParamsP()) << " primes in the special prime modulus." << std::endl;
    double logQ = 0;
    double logP = 0;
    {
        auto moduliQ = cc->GetCryptoParameters()->GetElementParams()->GetModulus();
        auto moduliP = cryptoParams->GetParamsP()->GetModulus();
        logQ         = moduliQ.GetMSB();
        logP         = moduliP.GetMSB();
    }
    std::cout << "log2(Q) = " << logQ << " log2(P) = " << logP << " log2(QP) = " << logQ + logP << std::endl;

    LeveledZ z     = std::make_shared<LeveledZImpl>();
    UserZ u        = std::make_shared<UserZImpl>(z);
    AdvancedZ advZ = std::make_shared<AdvancedZImpl>(z);
    FHEZ fheZ      = std::make_shared<FHEZImpl>(z, advZ);
    PKEZ pkeZ      = std::make_shared<PKEZImpl>(keyPair.publicKey, keyPair.secretKey);
    pkeZ_global    = pkeZ;  // Temporary

    uint32_t zSlots = cc->GetRingDimension() / zN;  // Full packing
    std::cout << "Bootstrapping parameters: zN = " << zN << ", zSlots = " << zSlots << std::endl << std::endl;

    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -16, 1);
    fheZ->EvalBootstrapKeyGen(keyPair.secretKey, zN, zSlots);

    //------
    // User program
    //------

    auto elemParam = cc->GetCryptoParameters()->GetElementParams();
    auto sfq0      = cryptoParams->GetScalingFactorBFP(0);

    auto oneTest = [&]() -> std::array<double, 6> {
        std::vector<BigInteger> vec(zSlots, 0);
        std::vector<BigInteger> vec2(zSlots, 0);
        for (size_t i = 0; i != zSlots; ++i) {
            vec[i]  = i + rand();  // between 0 and RAND_MAX = 2147483647
            vec2[i] = 2 * i + rand();
        }
        auto ptxt1 = ZEncodingImpl::encodeArith(vec, zN, zSlots, elemParam, sfq0);
        auto ptxt2 = ZEncodingImpl::encodeArith(vec2, zN, zSlots, elemParam, sfq0);

        /// TEST ENCODE
        auto ct  = pkeZ->Encrypt(ptxt1);
        auto ct2 = pkeZ->Encrypt(ptxt2);

        auto ctDec = pkeZ->Decrypt(ct);
        //ctDec.print("Arith ct");
        auto ct2Dec = pkeZ->Decrypt(ct2);
        //std::cout << std::endl;

        // Test ct-ct for freshly A2A ct

        auto ctA2A  = fheZ->EvalArithToArith(ct);
        auto ct2A2A = fheZ->EvalArithToArith(ct2);
        //
        auto ctA2ADec = pkeZ->Decrypt(ctA2A);
        //ctDec.print("Arith ct");
        auto ct2A2ADec = pkeZ->Decrypt(ct2A2A);

        auto ctMult    = u->EvalMultFullInZ(ctA2A, ct2A2A);
        auto ctMultDec = pkeZ->Decrypt(ctMult);

        auto valueMult0         = ctMultDec[0];
        auto valueMultExpected0 = (ctA2ADec[0] * ct2A2ADec[0]) % (BigInteger(1) << zN);
        if (valueMult0 != valueMultExpected0) {
            std::cout << "Error in Mult!" << std::endl;
        }

        auto noiseGrowth = ctMultDec.getLogMaxNoise() - std::max(ctA2ADec.getLogMaxNoise(), ct2A2ADec.getLogMaxNoise());
        auto overflowGrowth = ctMultDec.getLogMaxI() - std::max(ctA2ADec.getLogMaxI(), ct2A2ADec.getLogMaxI());
        auto overflowA      = ctA2ADec.getLogMaxI();
        auto overflowB      = ct2A2ADec.getLogMaxI();
        //std::cout << "Noise growth after Mult on A2A ct: " << noiseGrowth << " bits." << std::endl;

        // Test ct-pt mult, i.e. multshort

        // vec2 is 32bit
        auto pt = BigInteger(vec2[0]);
        if (zN > 32) {
            for (int i = 1; i < zN / 32; i++) {
                pt = (pt << 32) + (BigInteger(vec2[0]));
            }
        }

        auto ctMultShort = u->EvalMultPtInZ(ctA2A, pt);

        auto ctMultShortDec          = pkeZ->Decrypt(ctMultShort);
        auto valueMultShort0         = ctMultShortDec[0];
        auto valueMultShortExpected0 = (ctA2ADec[0] * pt) % (BigInteger(1) << zN);
        if (valueMultShort0 != valueMultShortExpected0) {
            std::cout << "Error in MultShort!" << std::endl;
        }

        auto noiseGrowthShort    = ctMultShortDec.getLogMaxNoise() - ctA2ADec.getLogMaxNoise();
        auto overflowGrowthShort = ctMultShortDec.getLogMaxI() - ctA2ADec.getLogMaxI();
        //std::cout << "Noise growth after MultShort on A2A ct: " << noiseGrowthShort << " bits." << std::endl;
        return {noiseGrowth, noiseGrowthShort, overflowGrowth, overflowGrowthShort, overflowA, overflowB};
    };

    std::vector<double> noiseGrowths;
    std::vector<double> noiseGrowthShorts;
    std::vector<double> overflowGrowths;
    std::vector<double> overflowGrowthShorts;
    std::vector<double> overflows;
    for (uint32_t i = 0; i < repeats; i++) {
        auto [a, b, c, d, e, f] = oneTest();
        noiseGrowths.push_back(a);
        noiseGrowthShorts.push_back(b);
        overflowGrowths.push_back(c);
        overflowGrowthShorts.push_back(d);
        overflows.push_back(e);
        overflows.push_back(f);
    }
    if (repeats >= 1) {
        double avgNoiseGrowth = std::accumulate(noiseGrowths.begin(), noiseGrowths.end(), 0.0) / noiseGrowths.size();
        double avgNoiseGrowthShort =
            std::accumulate(noiseGrowthShorts.begin(), noiseGrowthShorts.end(), 0.0) / noiseGrowthShorts.size();
        double maxNoiseGrowth      = *std::max_element(noiseGrowths.begin(), noiseGrowths.end());
        double maxNoiseGrowthShort = *std::max_element(noiseGrowthShorts.begin(), noiseGrowthShorts.end());

        std::cout << "\nMax noise growth after Mult on A2A ct over " << repeats << " runs: " << maxNoiseGrowth
                  << " bits." << std::endl;
        std::cout << "Max noise growth after MultShort on A2A ct over " << repeats << " runs: " << maxNoiseGrowthShort
                  << " bits." << std::endl;
        std::cout << "\nAverage noise growth after Mult on A2A ct over " << repeats << " runs: " << avgNoiseGrowth
                  << " bits." << std::endl;
        std::cout << "Average noise growth after MultShort on A2A ct over " << repeats
                  << " runs: " << avgNoiseGrowthShort << " bits." << std::endl;

        double maxOverflowGrowth      = *std::max_element(overflowGrowths.begin(), overflowGrowths.end());
        double maxOverflowGrowthShort = *std::max_element(overflowGrowthShorts.begin(), overflowGrowthShorts.end());
        double maxOverflow            = *std::max_element(overflows.begin(), overflows.end());

        double avgOverflowGrowth =
            std::accumulate(overflowGrowths.begin(), overflowGrowths.end(), 0.0) / overflowGrowths.size();
        double avgOverflowGrowthShort = std::accumulate(overflowGrowthShorts.begin(), overflowGrowthShorts.end(), 0.0) /
                                        overflowGrowthShorts.size();
        double avgOverflow = std::accumulate(overflows.begin(), overflows.end(), 0.0) / overflows.size();
        std::cout << "\nMax overflow growth after Mult on A2A ct over " << repeats << " runs: " << maxOverflowGrowth
                  << " bits." << std::endl;
        std::cout << "Max overflow growth after MultShort on A2A ct over " << repeats
                  << " runs: " << maxOverflowGrowthShort << " bits." << std::endl;
        std::cout << "Max overflow before Mult on A2A ct over " << repeats << " runs: " << maxOverflow << " bits."
                  << std::endl;

        std::cout << "\nAverage overflow growth after Mult on A2A ct over " << repeats << " runs: " << avgOverflowGrowth
                  << " bits." << std::endl;
        std::cout << "Average overflow growth after MultShort on A2A ct over " << repeats
                  << " runs: " << avgOverflowGrowthShort << " bits." << std::endl;
        std::cout << "Average overflow before Mult on A2A ct over " << repeats << " runs: " << avgOverflow << " bits."
                  << std::endl;
    }
}
