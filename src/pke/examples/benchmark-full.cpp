#include "openfhe.h"
#include "scheme/ckksrns/z-user.h"
#include "utils.h"
#include "scheme/ckksrns/z-fhe.h"
#include "scheme/ckksrns/z-pke.h"

using namespace lbcrypto;

void SimpleBootstrapExample(int zN, std::string directive);

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

    SimpleBootstrapExample(zN, directive);
}

void SimpleBootstrapExample(int zN, std::string directive) {
    CCParams<CryptoContextCKKSRNS> parameters;

    parameters.SetSecretKeyDist(lbcrypto::SPARSE_ENCAPSULATED);

    uint32_t repeats = 1;
    if (directive == "bench" || directive == "verify") {
        parameters.SetSecurityLevel(lbcrypto::HEStd_128_classic);
        parameters.SetRingDim(1 << 16);
        if (directive == "bench")
            repeats = 5;
        if (directive == "verify")
            repeats = 0;
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

    fheZ->EvalBootstrapSetup(*cc, zN, zSlots, levelBudget, {0, 0}, 4, -24, 1);
    fheZ->EvalBootstrapKeyGen(keyPair.secretKey, zN, zSlots);

    auto elemParam = cc->GetCryptoParameters()->GetElementParams();
    auto sfq0      = cryptoParams->GetScalingFactorBFP(0);

    std::vector<BigInteger> vec(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec[i] = i + 0xdeadbeaf;
    }
    auto ptxt1 = ZEncodingImpl::encodeArith(vec, zN, zSlots, elemParam, sfq0);

    std::vector<BigInteger> vec2(zSlots, 0);
    for (size_t i = 0; i != zSlots; ++i) {
        vec2[i] = i + 0xf0f0f0ff;
    }
    auto ptxt2 = ZEncodingImpl::encodeArith(vec2, zN, zSlots, elemParam, sfq0);

    /// TEST ENCODE
    auto ct  = pkeZ->Encrypt(ptxt1);
    auto ct2 = pkeZ->Encrypt(ptxt2);

    auto ctDec = pkeZ->Decrypt(ct);
    ctDec.print("Arith ct");
    auto ct2Dec = pkeZ->Decrypt(ct2);

    /// TEST Boolean Encode
    auto ptxt3          = ZEncodingImpl::encodeBooleanFull(vec, zN, zSlots, elemParam, sfq0);
    CiphertextGroup ct3 = pkeZ->Encrypt(ptxt3);
    auto ct3Dec         = pkeZ->Decrypt(ct3);
    ct3Dec.print("Bool ct");

    auto ptxt4          = ZEncodingImpl::encodeBooleanFull(vec2, zN, zSlots, elemParam, sfq0);
    CiphertextGroup ct4 = pkeZ->Encrypt(ptxt4);
    auto ct4Dec         = pkeZ->Decrypt(ct4);

    std::cout << std::endl;

    // Mult
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = u->EvalMultFullInZ(ct, ct);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (ctResDec[0] != (ctDec[0] * ctDec[0]) % (BigInteger(1) << zN)) {
                       std::cout << "Error in MultFull!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ctDec, "MultFull");
                   return;
               }()),
               "MultFull");
        BENCHMARK(u->EvalMultFullInZ(ct, ct), repeats, "MultFull");
    }

    // Bool
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = u->EvalBooleanOR(ct3, ct4);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (ctResDec[0].ConvertToInt() != (ct3Dec[0].ConvertToInt() | ct4Dec[0].ConvertToInt())) {
                       std::cout << "Error in BooleanOR!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ctDec, "BooleanOR");
                   return;
               }()),
               "BooleanOR");
        BENCHMARK(u->EvalBooleanOR(ct3, ct4), repeats, "BooleanOR");
    }

    // BoolToArith
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = fheZ->EvalBooleanToArith(ct3);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (!ctResDec.valuesEqual(ct3Dec)) {
                       std::cout << "Error in B2A!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ctDec, "B2A");
                   return;
               }()),
               "B2A");
        BENCHMARK(fheZ->EvalBooleanToArith(ct3), repeats, "B2A");
    }

    // BooleanToBoolean
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = fheZ->EvalBooleanToBooleanFull(ct3);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (!ctResDec.valuesEqual(ct3Dec)) {
                       std::cout << "Error in B2B!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ct3Dec, "B2B");
                   return;
               }()),
               "B2B");
        BENCHMARK((fheZ->EvalBooleanToBooleanFull(ct3)), repeats, "B2B");
    }

    // ArithToArithHigh
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = fheZ->EvalArithToArithHigh(ct2);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (!ctResDec.valuesEqual(ct2Dec)) {
                       std::cout << "Error in A2AI!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ct2Dec, "A2AI");
                   return;
               }()),
               "A2AI");
        BENCHMARK(fheZ->EvalArithToArithHigh(ct2), repeats, "A2AI");
    }

    // ArithToArithNoise
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = fheZ->EvalArithToArithNoise(ct2);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (!ctResDec.valuesEqual(ct2Dec)) {
                       std::cout << "Error in A2Ae!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ct2Dec, "A2Ae");
                   return;
               }()),
               "A2Ae");
        BENCHMARK(fheZ->EvalArithToArithNoise(ct2), repeats, "A2Ae");
    }

    // ArithToArith
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = fheZ->EvalArithToArith(ct2);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (!ctResDec.valuesEqual(ct2Dec)) {
                       std::cout << "Error in A2A!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ct2Dec, "A2A");
                   return;
               }()),
               "A2A");
        BENCHMARK(fheZ->EvalArithToArith(ct2), repeats, "A2A");
    }

    // ArithToBooleanBatched
    if (1) {
        std::vector<Ciphertext<DCRTPoly>> batchCts = {ct, ct2};
        size_t batchSize                             = zN / 4;  // zN / w
        while (batchCts.size() < batchSize) {
            batchCts.push_back(ct);
        }

        WARMUP(([&]() {
                   auto ctGroupBool = fheZ->EvalArithToBooleanBatched(batchCts);
                   auto ctResDec0   = pkeZ->Decrypt(CiphertextGroup({ctGroupBool[0], ctGroupBool[1]}));
                   auto ctResDec1   = pkeZ->Decrypt(CiphertextGroup({ctGroupBool[2], ctGroupBool[3]}));
                   if (!ctResDec0.valuesEqual(ctDec)) {
                       std::cout << "Error in B-A2B!" << std::endl;
                   }
                   if (!ctResDec1.valuesEqual(ct2Dec)) {
                       std::cout << "Error in B-A2B!" << std::endl;
                   }
                   ctResDec0.printNoiseComparison(ctDec, "B-A2B-0");
                   ctResDec1.printNoiseComparison(ct2Dec, "B-A2B-1");
                   return;
               }()),
               "B-A2B");
        BENCHMARK(fheZ->EvalArithToBooleanBatched(batchCts), repeats, "B-A2B");
    }

    //ArithToBooleanFull
    if (1) {
        WARMUP(([&]() {
                   auto ctRes    = fheZ->EvalArithToBooleanFull(ct2);
                   auto ctResDec = pkeZ->Decrypt(ctRes);
                   if (!ctResDec.valuesEqual(ct2Dec)) {
                       std::cout << "Error in A2B!" << std::endl;
                   }
                   ctResDec.printNoiseComparison(ct2Dec, "A2B");
                   return;
               }()),
               "A2B");
        BENCHMARK(fheZ->EvalArithToBooleanFull(ct2), repeats, "A2B");
    }
}
