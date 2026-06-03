// Final compact version: MPI CRT + point-to-point collection + intra-rank threaded CRT NTT

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <pthread.h>
#include <string>
#include <vector>

#ifndef THREAD_COUNT
#define THREAD_COUNT 8
#endif

#ifndef ORDINARY_MPI_TASK
#define ORDINARY_MPI_TASK 0
#endif

#ifndef CRT_COLLECT_METHOD
#define CRT_COLLECT_METHOD 0
#endif

#ifndef CRT_INTRA_THREAD
#define CRT_INTRA_THREAD 0
#endif

#if ORDINARY_MPI_TASK != 0 && ORDINARY_MPI_TASK != 1
#error "ORDINARY_MPI_TASK must be 0 or 1"
#endif

#if CRT_COLLECT_METHOD != 0 && CRT_COLLECT_METHOD != 1
#error "CRT_COLLECT_METHOD must be 0 or 1"
#endif

#if CRT_INTRA_THREAD != 0 && CRT_INTRA_THREAD != 1
#error "CRT_INTRA_THREAD must be 0 or 1"
#endif

// 鍙互鑷娣诲姞闇€瑕佺殑澶存枃浠?
// 鏈増鏈彧淇濈暀鏈€缁堢湡姝ｄ娇鐢ㄧ殑涓ゆ潯涓荤嚎锛?
// 1. input 0~3锛歅thread whole-pipeline NTT
// 2. input 4锛?-mod CRT + Montgomery NTT + Fast Garner

//======================================================================
// 鏁版嵁杈撳叆銆佹鏌ュ拰璋冭瘯杈撳嚭鍑芥暟
//======================================================================

void fRead(int *a, int *b, int *n, int *p, int input_id) {
    // 鏁版嵁杈撳叆鍑芥暟
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";

    std::ifstream fin;
    fin.open(strin, std::ios::in);

    fin >> *n >> *p;

    for (int i = 0; i < *n; i++) {
        fin >> a[i];
    }

    for (int i = 0; i < *n; i++) {
        fin >> b[i];
    }
}

void fCheck(int *ab, int n, int input_id) {
    // 鍒ゆ柇澶氶」寮忎箻娉曠粨鏋滄槸鍚︽纭?
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";

    std::ifstream fin;
    fin.open(strout, std::ios::in);

    for (int i = 0; i < n * 2 - 1; i++) {
        int x;
        fin >> x;

        if (x != ab[i]) {
            std::cout << "澶氶」寮忎箻娉曠粨鏋滈敊璇? << std::endl;
            return;
        }
    }

    std::cout << "澶氶」寮忎箻娉曠粨鏋滄纭? << std::endl;
}

void fWrite(int *ab, int n, int input_id) {
    // 鏁版嵁杈撳嚭鍑芥暟, 鍙互鐢ㄦ潵杈撳嚭鏈€缁堢粨鏋? 涔熷彲鐢ㄤ簬璋冭瘯鏃惰緭鍑轰腑闂存暟缁?
    // 浣跨敤鍓嶈纭褰撳墠鐩綍涓嬪凡缁忓瓨鍦?files 鏂囦欢澶?
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";

    std::ofstream fout;
    fout.open(strout, std::ios::out);

    for (int i = 0; i < n * 2 - 1; i++) {
        fout << ab[i] << '\n';
    }
}

void fReadLarge(long long *a, long long *b, int *n, long long *p, int input_id) {
    // 澶фā鏁拌緭鍏ュ嚱鏁?
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";

    std::ifstream fin;
    fin.open(strin, std::ios::in);

    fin >> *n >> *p;

    for (int i = 0; i < *n; i++) {
        fin >> a[i];
    }

    for (int i = 0; i < *n; i++) {
        fin >> b[i];
    }
}

void fCheckLarge(long long *ab, int n, int input_id) {
    // 鍒ゆ柇澶фā鏁板椤瑰紡涔樻硶缁撴灉鏄惁姝ｇ‘
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";

    std::ifstream fin;
    fin.open(strout, std::ios::in);

    for (int i = 0; i < n * 2 - 1; i++) {
        long long x;
        fin >> x;

        if (x != ab[i]) {
            std::cout << "澶фā鏁板椤瑰紡涔樻硶缁撴灉閿欒" << std::endl;
            return;
        }
    }

    std::cout << "澶фā鏁板椤瑰紡涔樻硶缁撴灉姝ｇ‘" << std::endl;
}

void fWriteLarge(long long *ab, int n, int input_id) {
    // 澶фā鏁版暟鎹緭鍑哄嚱鏁? 璋冭瘯 input 4 鏃跺彲浠ヤ娇鐢?
    // 浣跨敤鍓嶈纭褰撳墠鐩綍涓嬪凡缁忓瓨鍦?files 鏂囦欢澶?
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";

    std::ofstream fout;
    fout.open(strout, std::ios::out);

    for (int i = 0; i < n * 2 - 1; i++) {
        fout << ab[i] << '\n';
    }
}

//======================================================================
// 閫氱敤蹇€熷箓鍜?Montgomery 妯′箻
//======================================================================

long long quickPowerMod(long long base, long long exponent, long long modulus) {
    // 蹇€熷箓锛氱敤浜庤绠?base^exponent mod modulus
    // 鍦?NTT 涓敤浜庤绠楀崟浣嶆牴銆侀€嗗崟浣嶆牴鍜岄暱搴︾殑妯￠€嗗厓
    long long result = 1;

    while (exponent > 0) {
        if (exponent & 1LL) {
            result = result * base % modulus;
        }

        base = base * base % modulus;
        exponent >>= 1LL;
    }

    return result;
}

static inline uint32_t computeMontgomeryNegInverse32(uint32_t modulus) {
    // 璁＄畻 -modulus^{-1} mod 2^32
    // Newton 杩唬姣忚疆浼氳姝ｇ‘浣嶆暟缈诲€?
    uint32_t inverse = modulus;

    for (int i = 0; i < 5; ++i) {
        inverse *= 2u - modulus * inverse;
    }

    return 0u - inverse;
}

static inline uint32_t montgomeryReduce32(
    uint64_t value,
    uint32_t modulus,
    uint32_t modulusNegInverse
) {
    // Montgomery REDC锛孯 = 2^32
    // 杩斿洖 value * R^{-1} mod modulus
    uint32_t q = static_cast<uint32_t>(value) * modulusNegInverse;
    uint64_t reduced = (value + static_cast<uint64_t>(q) * modulus) >> 32;

    uint32_t result = static_cast<uint32_t>(reduced);

    if (result >= modulus) {
        result -= modulus;
    }

    return result;
}

//======================================================================
// 32 浣嶆ā鏁帮細Pthread whole-pipeline NTT
//======================================================================

struct PthreadWholeNttContext {
    int *firstInput;
    int *secondInput;
    int *resultOutput;

    int polynomialLength;
    int transformLength;
    int modulus;
    int threadCount;

    int *firstTransformed;
    int *secondTransformed;

    uint32_t modulusUnsigned;
    uint32_t montgomeryNegInverse;
    uint32_t montgomeryR;

    int *stageRoots;
    int *stageRootsMontgomery;

    long long inverseLength;

    pthread_barrier_t barrier;
};

struct PthreadWholeWorkerParam {
    int threadId;
    PthreadWholeNttContext *context;
};

void bitReverseSerial(int *data, int transformLength) {
    // bit-reversal permutation
    // 杩欓噷淇濇寔涓茶锛岀敱 thread 0 鎵ц锛岄伩鍏嶅苟琛?swap 甯︽潵鐨勫啓鍐茬獊
    for (int currentIndex = 1, reversedIndex = 0;
         currentIndex < transformLength;
         ++currentIndex) {

        int highestBit = transformLength >> 1;

        while (reversedIndex & highestBit) {
            reversedIndex ^= highestBit;
            highestBit >>= 1;
        }

        reversedIndex ^= highestBit;

        if (currentIndex < reversedIndex) {
            std::swap(data[currentIndex], data[reversedIndex]);
        }
    }
}

void pthreadWholeRunNttStages(
    PthreadWholeNttContext *ctx,
    int *data,
    bool inverseTransform,
    int threadId
) {
    const int primitiveRoot = 3;

    for (int butterflySize = 1;
         butterflySize < ctx->transformLength;
         butterflySize <<= 1) {

        int blockSize = butterflySize << 1;

        if (threadId == 0) {
            long long rootStep = quickPowerMod(
                primitiveRoot,
                (ctx->modulus - 1) / blockSize,
                ctx->modulus
            );

            if (inverseTransform) {
                rootStep = quickPowerMod(rootStep, ctx->modulus - 2, ctx->modulus);
            }

            // stageRoots 棰勮绠楋細姣忎竴灞傜殑鍗曚綅鏍瑰彧绠椾竴娆?
            ctx->stageRoots[0] = 1;

            for (int offset = 1; offset < butterflySize; ++offset) {
                ctx->stageRoots[offset] = static_cast<int>(
                    1LL * ctx->stageRoots[offset - 1] * rootStep % ctx->modulus
                );
            }

            // 杞垚 Montgomery 褰㈠紡锛屽噺灏?butterfly 涓殑鍙栨ā寮€閿€
            for (int offset = 0; offset < butterflySize; ++offset) {
                ctx->stageRootsMontgomery[offset] = static_cast<int>(
                    static_cast<uint64_t>(ctx->stageRoots[offset]) *
                    ctx->montgomeryR %
                    ctx->modulusUnsigned
                );
            }
        }

        // 绛夊緟 thread 0 瀹屾垚褰撳墠 stage 鐨勫崟浣嶆牴棰勮绠?
        pthread_barrier_wait(&ctx->barrier);

        int blockCount = ctx->transformLength / blockSize;

        // 鍜?OpenMP schedule(static) 绫讳技锛氭瘡涓嚎绋嬪鐞嗕竴娈佃繛缁?block
        // 杩欐牱瀹炵幇绠€鍗曪紝鏁版嵁璁块棶涔熸瘮杈冭繛缁?
        int beginBlock = static_cast<int>(
            1LL * blockCount * threadId / ctx->threadCount
        );

        int endBlock = static_cast<int>(
            1LL * blockCount * (threadId + 1) / ctx->threadCount
        );

        for (int blockIndex = beginBlock; blockIndex < endBlock; ++blockIndex) {
            int blockStart = blockIndex * blockSize;

            for (int offset = 0; offset < butterflySize; ++offset) {
                int upperIndex = blockStart + offset;
                int lowerIndex = upperIndex + butterflySize;

                int upperValue = data[upperIndex];

                uint64_t product =
                    static_cast<uint64_t>(
                        static_cast<uint32_t>(data[lowerIndex])
                    ) *
                    static_cast<uint32_t>(ctx->stageRootsMontgomery[offset]);

                int lowerValue = static_cast<int>(
                    montgomeryReduce32(
                        product,
                        ctx->modulusUnsigned,
                        ctx->montgomeryNegInverse
                    )
                );

                int addedValue = upperValue + lowerValue;
                if (addedValue >= ctx->modulus) {
                    addedValue -= ctx->modulus;
                }

                int subtractedValue = upperValue - lowerValue;
                if (subtractedValue < 0) {
                    subtractedValue += ctx->modulus;
                }

                data[upperIndex] = addedValue;
                data[lowerIndex] = subtractedValue;
            }
        }

        // 绛夊緟褰撳墠 stage 鍏ㄩ儴瀹屾垚锛屾墠鑳借繘鍏ヤ笅涓€ stage
        pthread_barrier_wait(&ctx->barrier);
    }

    if (inverseTransform) {
        if (threadId == 0) {
            ctx->inverseLength = quickPowerMod(
                ctx->transformLength,
                ctx->modulus - 2,
                ctx->modulus
            );
        }

        pthread_barrier_wait(&ctx->barrier);

        int beginIndex = static_cast<int>(
            1LL * ctx->transformLength * threadId / ctx->threadCount
        );

        int endIndex = static_cast<int>(
            1LL * ctx->transformLength * (threadId + 1) / ctx->threadCount
        );

        for (int index = beginIndex; index < endIndex; ++index) {
            data[index] = static_cast<int>(
                1LL * data[index] * ctx->inverseLength % ctx->modulus
            );
        }

        pthread_barrier_wait(&ctx->barrier);
    }
}

void *pthreadWholeNttWorker(void *arg) {
    PthreadWholeWorkerParam *param = static_cast<PthreadWholeWorkerParam *>(arg);
    int threadId = param->threadId;
    PthreadWholeNttContext *ctx = param->context;

    int beginIndex = static_cast<int>(
        1LL * ctx->transformLength * threadId / ctx->threadCount
    );

    int endIndex = static_cast<int>(
        1LL * ctx->transformLength * (threadId + 1) / ctx->threadCount
    );

    // 1. 娓呴浂鍜岃緭鍏ユ嫹璐濆苟琛屽寲
    for (int index = beginIndex; index < endIndex; ++index) {
        ctx->firstTransformed[index] = 0;
        ctx->secondTransformed[index] = 0;

        if (index < ctx->polynomialLength) {
            ctx->firstTransformed[index] = ctx->firstInput[index] % ctx->modulus;
            ctx->secondTransformed[index] = ctx->secondInput[index] % ctx->modulus;
        }
    }

    pthread_barrier_wait(&ctx->barrier);

    // 2. first 姝?NTT
    if (threadId == 0) {
        bitReverseSerial(ctx->firstTransformed, ctx->transformLength);
    }

    pthread_barrier_wait(&ctx->barrier);
    pthreadWholeRunNttStages(ctx, ctx->firstTransformed, false, threadId);

    // 3. second 姝?NTT
    if (threadId == 0) {
        bitReverseSerial(ctx->secondTransformed, ctx->transformLength);
    }

    pthread_barrier_wait(&ctx->barrier);
    pthreadWholeRunNttStages(ctx, ctx->secondTransformed, false, threadId);

    // 4. 鐐瑰€间箻娉曞苟琛屽寲
    beginIndex = static_cast<int>(
        1LL * ctx->transformLength * threadId / ctx->threadCount
    );

    endIndex = static_cast<int>(
        1LL * ctx->transformLength * (threadId + 1) / ctx->threadCount
    );

    for (int index = beginIndex; index < endIndex; ++index) {
        ctx->firstTransformed[index] = static_cast<int>(
            1LL * ctx->firstTransformed[index] *
            ctx->secondTransformed[index] %
            ctx->modulus
        );
    }

    pthread_barrier_wait(&ctx->barrier);

    // 5. 閫?NTT
    if (threadId == 0) {
        bitReverseSerial(ctx->firstTransformed, ctx->transformLength);
    }

    pthread_barrier_wait(&ctx->barrier);
    pthreadWholeRunNttStages(ctx, ctx->firstTransformed, true, threadId);

    // 6. 鎷疯礉鐪熷疄缁撴灉骞惰鍖?
    int resultLength = 2 * ctx->polynomialLength - 1;

    int beginResult = static_cast<int>(
        1LL * resultLength * threadId / ctx->threadCount
    );

    int endResult = static_cast<int>(
        1LL * resultLength * (threadId + 1) / ctx->threadCount
    );

    for (int index = beginResult; index < endResult; ++index) {
        ctx->resultOutput[index] = ctx->firstTransformed[index];
    }

    pthread_barrier_wait(&ctx->barrier);

    return nullptr;
}

void multiplyByNttPthreadWhole(
    int *firstPolynomial,
    int *secondPolynomial,
    int *resultPolynomial,
    int polynomialLength,
    int modulus,
    int threadCount
) {
    // 浣跨敤 static 鏁扮粍閬垮厤澶ф暟缁勫弽澶嶅湪鏍堜笂鐢宠锛屽噺灏戞爤婧㈠嚭椋庨櫓
    static int firstTransformed[300000];
    static int secondTransformed[300000];
    static int stageRoots[300000];
    static int stageRootsMontgomery[300000];

    int transformLength = 1;

    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    PthreadWholeNttContext context;

    context.firstInput = firstPolynomial;
    context.secondInput = secondPolynomial;
    context.resultOutput = resultPolynomial;

    context.polynomialLength = polynomialLength;
    context.transformLength = transformLength;
    context.modulus = modulus;
    context.threadCount = threadCount;

    context.firstTransformed = firstTransformed;
    context.secondTransformed = secondTransformed;

    context.modulusUnsigned = static_cast<uint32_t>(modulus);
    context.montgomeryNegInverse = computeMontgomeryNegInverse32(context.modulusUnsigned);
    context.montgomeryR = static_cast<uint32_t>((1ULL << 32) % context.modulusUnsigned);

    context.stageRoots = stageRoots;
    context.stageRootsMontgomery = stageRootsMontgomery;
    context.inverseLength = 1;

    pthread_barrier_init(&context.barrier, nullptr, threadCount);

    std::vector<pthread_t> handles(threadCount);
    std::vector<PthreadWholeWorkerParam> params(threadCount);

    // 鏁翠釜澶氶」寮忎箻娉曞彧鍒涘缓涓€娆＄嚎绋嬬粍
    // 娓呴浂銆佷袱娆℃ NTT銆佺偣鍊间箻娉曘€侀€?NTT銆佺粨鏋滃洖鎷烽兘鍦ㄥ悓涓€鎵圭嚎绋嬮噷瀹屾垚
    for (int threadId = 0; threadId < threadCount; ++threadId) {
        params[threadId].threadId = threadId;
        params[threadId].context = &context;

        pthread_create(
            &handles[threadId],
            nullptr,
            pthreadWholeNttWorker,
            &params[threadId]
        );
    }

    for (int threadId = 0; threadId < threadCount; ++threadId) {
        pthread_join(handles[threadId], nullptr);
    }

    pthread_barrier_destroy(&context.barrier);
}


//======================================================================
// 32-bit modulus: task-level MPI NTT
//======================================================================

void numberTheoreticTransformOrdinarySerial(
    int *data,
    int transformLength,
    int modulus,
    bool inverseTransform
) {
    const int primitiveRoot = 3;

    bitReverseSerial(data, transformLength);

    for (int butterflySize = 1;
         butterflySize < transformLength;
         butterflySize <<= 1) {

        int blockSize = butterflySize << 1;
        long long rootStep = quickPowerMod(
            primitiveRoot,
            (modulus - 1) / blockSize,
            modulus
        );

        if (inverseTransform) {
            rootStep = quickPowerMod(rootStep, modulus - 2, modulus);
        }

        for (int blockStart = 0;
             blockStart < transformLength;
             blockStart += blockSize) {

            long long currentRoot = 1;

            for (int offset = 0; offset < butterflySize; ++offset) {
                int upperIndex = blockStart + offset;
                int lowerIndex = upperIndex + butterflySize;

                int upperValue = data[upperIndex];
                int lowerValue = static_cast<int>(
                    currentRoot * data[lowerIndex] % modulus
                );

                int addedValue = upperValue + lowerValue;
                if (addedValue >= modulus) {
                    addedValue -= modulus;
                }

                int subtractedValue = upperValue - lowerValue;
                if (subtractedValue < 0) {
                    subtractedValue += modulus;
                }

                data[upperIndex] = addedValue;
                data[lowerIndex] = subtractedValue;

                currentRoot = currentRoot * rootStep % modulus;
            }
        }
    }

    if (inverseTransform) {
        long long inverseLength = quickPowerMod(
            transformLength,
            modulus - 2,
            modulus
        );

        for (int index = 0; index < transformLength; ++index) {
            data[index] = static_cast<int>(
                data[index] * inverseLength % modulus
            );
        }
    }
}

void multiplyByNttTaskMPI(
    int *firstPolynomial,
    int *secondPolynomial,
    int *resultPolynomial,
    int polynomialLength,
    int modulus,
    int threadCount,
    int rank,
    int size
) {
    if (size < 2) {
        if (rank == 0) {
            multiplyByNttPthreadWhole(
                firstPolynomial,
                secondPolynomial,
                resultPolynomial,
                polynomialLength,
                modulus,
                threadCount
            );
        }
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }

    int transformLength = 1;
    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    int resultLength = 2 * polynomialLength - 1;

    std::vector<int> transformed(transformLength, 0);

    if (rank == 0) {
        for (int index = 0; index < polynomialLength; ++index) {
            transformed[index] = firstPolynomial[index] % modulus;
        }

        numberTheoreticTransformOrdinarySerial(
            transformed.data(),
            transformLength,
            modulus,
            false
        );

        std::vector<int> secondTransformed(transformLength, 0);
        MPI_Recv(
            secondTransformed.data(),
            transformLength,
            MPI_INT,
            1,
            1001,
            MPI_COMM_WORLD,
            MPI_STATUS_IGNORE
        );

        for (int index = 0; index < transformLength; ++index) {
            transformed[index] = static_cast<int>(
                1LL * transformed[index] *
                secondTransformed[index] %
                modulus
            );
        }

        numberTheoreticTransformOrdinarySerial(
            transformed.data(),
            transformLength,
            modulus,
            true
        );

        for (int index = 0; index < resultLength; ++index) {
            resultPolynomial[index] = transformed[index];
        }
    } else if (rank == 1) {
        for (int index = 0; index < polynomialLength; ++index) {
            transformed[index] = secondPolynomial[index] % modulus;
        }

        numberTheoreticTransformOrdinarySerial(
            transformed.data(),
            transformLength,
            modulus,
            false
        );

        MPI_Send(
            transformed.data(),
            transformLength,
            MPI_INT,
            0,
            1001,
            MPI_COMM_WORLD
        );
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

//======================================================================
// 澶фā鏁?input 4锛?-mod CRT + Montgomery NTT + Fast Garner
//======================================================================

const int CRT_MOD_COUNT = 4;

const int CRT_MODS[CRT_MOD_COUNT] = {
    469762049,
    754974721,
    880803841,
    998244353
};

const int CRT_ROOTS[CRT_MOD_COUNT] = {
    3,
    11,
    26,
    3
};

static int crtWorkA[CRT_MOD_COUNT][300000];
static int crtWorkB[CRT_MOD_COUNT][300000];
static int crtStageRootsMontgomery[CRT_MOD_COUNT][300000];
static int crtInverseLengthMontgomery[CRT_MOD_COUNT];

static inline int multiplyModCrtInt(int x, int y, int modulus) {
    return static_cast<int>(1LL * x * y % modulus);
}

int quickPowerModCrtInt(int base, long long exponent, int modulus) {
    // CRT 灏忔ā鏁板揩閫熷箓
    long long result = 1;
    long long current = base % modulus;

    while (exponent > 0) {
        if (exponent & 1LL) {
            result = result * current % modulus;
        }

        current = current * current % modulus;
        exponent >>= 1LL;
    }

    return static_cast<int>(result);
}

static inline int positiveModLongLong(long long value, int modulus) {
    value %= modulus;

    if (value < 0) {
        value += modulus;
    }

    return static_cast<int>(value);
}

void numberTheoreticTransformCrtModMontgomery(
    int *data,
    int transformLength,
    int modulus,
    int primitiveRoot,
    bool inverseTransform,
    int *stageRootsMontgomery
) {
    // 鍗曚釜 CRT 灏忔ā鏁颁笅鐨?NTT锛屼緵 V10-1 fallback 浣跨敤
    for (int currentIndex = 1, reversedIndex = 0;
         currentIndex < transformLength;
         ++currentIndex) {

        int highestBit = transformLength >> 1;

        while (reversedIndex & highestBit) {
            reversedIndex ^= highestBit;
            highestBit >>= 1;
        }

        reversedIndex ^= highestBit;

        if (currentIndex < reversedIndex) {
            std::swap(data[currentIndex], data[reversedIndex]);
        }
    }

    const uint32_t modulusUnsigned = static_cast<uint32_t>(modulus);
    const uint32_t montgomeryNegInverse = computeMontgomeryNegInverse32(modulusUnsigned);
    const uint32_t montgomeryR = static_cast<uint32_t>((1ULL << 32) % modulusUnsigned);

    for (int butterflySize = 1;
         butterflySize < transformLength;
         butterflySize <<= 1) {

        int blockSize = butterflySize << 1;

        int rootStep = quickPowerModCrtInt(
            primitiveRoot,
            (modulus - 1) / blockSize,
            modulus
        );

        if (inverseTransform) {
            rootStep = quickPowerModCrtInt(rootStep, modulus - 2, modulus);
        }

        int currentRoot = 1;

        for (int offset = 0; offset < butterflySize; ++offset) {
            stageRootsMontgomery[offset] = static_cast<int>(
                static_cast<uint64_t>(currentRoot) *
                montgomeryR %
                modulusUnsigned
            );

            currentRoot = multiplyModCrtInt(currentRoot, rootStep, modulus);
        }

        for (int blockStart = 0;
             blockStart < transformLength;
             blockStart += blockSize) {

            for (int offset = 0; offset < butterflySize; ++offset) {
                int upperIndex = blockStart + offset;
                int lowerIndex = upperIndex + butterflySize;

                int upperValue = data[upperIndex];

                uint64_t product =
                    static_cast<uint64_t>(
                        static_cast<uint32_t>(data[lowerIndex])
                    ) *
                    static_cast<uint32_t>(stageRootsMontgomery[offset]);

                int lowerValue = static_cast<int>(
                    montgomeryReduce32(product, modulusUnsigned, montgomeryNegInverse)
                );

                int addedValue = upperValue + lowerValue;
                if (addedValue >= modulus) {
                    addedValue -= modulus;
                }

                int subtractedValue = upperValue - lowerValue;
                if (subtractedValue < 0) {
                    subtractedValue += modulus;
                }

                data[upperIndex] = addedValue;
                data[lowerIndex] = subtractedValue;
            }
        }
    }

    if (inverseTransform) {
        int inverseLength = quickPowerModCrtInt(
            transformLength,
            modulus - 2,
            modulus
        );

        int inverseLengthMontgomery = static_cast<int>(
            static_cast<uint64_t>(inverseLength) *
            montgomeryR %
            modulusUnsigned
        );

        for (int index = 0; index < transformLength; ++index) {
            uint64_t product =
                static_cast<uint64_t>(
                    static_cast<uint32_t>(data[index])
                ) *
                static_cast<uint32_t>(inverseLengthMontgomery);

            data[index] = static_cast<int>(
                montgomeryReduce32(product, modulusUnsigned, montgomeryNegInverse)
            );
        }
    }
}

void numberTheoreticTransformCrtModGroupedOpenMP(
    int *data,
    int transformLength,
    int modulusIndex,
    bool inverseTransform,
    int localThreadId,
    int localThreadCount
) {
    // 鍒嗙粍 CRT NTT锛歍HREAD_COUNT = 8 鏃讹紝姣忎釜 CRT 妯℃暟鍒嗗埌 2 涓?OpenMP 绾跨▼
    // 娉ㄦ剰杩欓噷鐨?barrier 鏄暣涓?OpenMP parallel 鍖哄煙鐨勫叏灞€ barrier
    // 鍥犱负 4 涓ā鏁扮殑 transformLength 鍜?stage 缁撴瀯瀹屽叏鐩稿悓锛屾墍浠ユ墍鏈夌嚎绋嬩細鎸夊悓鏍烽『搴忓埌杈?barrier
    int modulus = CRT_MODS[modulusIndex];
    int primitiveRoot = CRT_ROOTS[modulusIndex];

    const uint32_t modulusUnsigned = static_cast<uint32_t>(modulus);
    const uint32_t montgomeryNegInverse = computeMontgomeryNegInverse32(modulusUnsigned);
    const uint32_t montgomeryR = static_cast<uint32_t>((1ULL << 32) % modulusUnsigned);

    if (localThreadId == 0) {
        for (int currentIndex = 1, reversedIndex = 0;
             currentIndex < transformLength;
             ++currentIndex) {

            int highestBit = transformLength >> 1;

            while (reversedIndex & highestBit) {
                reversedIndex ^= highestBit;
                highestBit >>= 1;
            }

            reversedIndex ^= highestBit;

            if (currentIndex < reversedIndex) {
                std::swap(data[currentIndex], data[reversedIndex]);
            }
        }
    }

    #pragma omp barrier

    for (int butterflySize = 1;
         butterflySize < transformLength;
         butterflySize <<= 1) {

        int blockSize = butterflySize << 1;

        if (localThreadId == 0) {
            int rootStep = quickPowerModCrtInt(
                primitiveRoot,
                (modulus - 1) / blockSize,
                modulus
            );

            if (inverseTransform) {
                rootStep = quickPowerModCrtInt(rootStep, modulus - 2, modulus);
            }

            int currentRoot = 1;

            for (int offset = 0; offset < butterflySize; ++offset) {
                crtStageRootsMontgomery[modulusIndex][offset] = static_cast<int>(
                    static_cast<uint64_t>(currentRoot) *
                    montgomeryR %
                    modulusUnsigned
                );

                currentRoot = multiplyModCrtInt(currentRoot, rootStep, modulus);
            }
        }

        #pragma omp barrier

        int totalButterflies = transformLength >> 1;

        int beginButterfly = static_cast<int>(
            1LL * totalButterflies * localThreadId / localThreadCount
        );

        int endButterfly = static_cast<int>(
            1LL * totalButterflies * (localThreadId + 1) / localThreadCount
        );

        // flat-butterfly 鍒掑垎锛氭瘡灞?butterfly 鎬绘暟鎭掍负 transformLength / 2
        // 杩欐牱鍙互閬垮厤鍚庢湡 stage block 鏁板お灏戝鑷寸嚎绋嬬┖杞?
        int blockIndex = beginButterfly / butterflySize;
        int offset = beginButterfly - blockIndex * butterflySize;
        int blockStart = blockIndex * blockSize;

        for (int butterflyId = beginButterfly;
             butterflyId < endButterfly;
             ++butterflyId) {

            int upperIndex = blockStart + offset;
            int lowerIndex = upperIndex + butterflySize;

            int upperValue = data[upperIndex];

            uint64_t product =
                static_cast<uint64_t>(
                    static_cast<uint32_t>(data[lowerIndex])
                ) *
                static_cast<uint32_t>(
                    crtStageRootsMontgomery[modulusIndex][offset]
                );

            int lowerValue = static_cast<int>(
                montgomeryReduce32(product, modulusUnsigned, montgomeryNegInverse)
            );

            int addedValue = upperValue + lowerValue;
            if (addedValue >= modulus) {
                addedValue -= modulus;
            }

            int subtractedValue = upperValue - lowerValue;
            if (subtractedValue < 0) {
                subtractedValue += modulus;
            }

            data[upperIndex] = addedValue;
            data[lowerIndex] = subtractedValue;

            ++offset;

            if (offset == butterflySize) {
                offset = 0;
                blockStart += blockSize;
            }
        }

        #pragma omp barrier
    }

    if (inverseTransform) {
        if (localThreadId == 0) {
            int inverseLength = quickPowerModCrtInt(
                transformLength,
                modulus - 2,
                modulus
            );

            crtInverseLengthMontgomery[modulusIndex] = static_cast<int>(
                static_cast<uint64_t>(inverseLength) *
                montgomeryR %
                modulusUnsigned
            );
        }

        #pragma omp barrier

        int beginIndex = static_cast<int>(
            1LL * transformLength * localThreadId / localThreadCount
        );

        int endIndex = static_cast<int>(
            1LL * transformLength * (localThreadId + 1) / localThreadCount
        );

        for (int index = beginIndex; index < endIndex; ++index) {
            uint64_t product =
                static_cast<uint64_t>(
                    static_cast<uint32_t>(data[index])
                ) *
                static_cast<uint32_t>(
                    crtInverseLengthMontgomery[modulusIndex]
                );

            data[index] = static_cast<int>(
                montgomeryReduce32(product, modulusUnsigned, montgomeryNegInverse)
            );
        }

        #pragma omp barrier
    }
}

struct CrtGarnerPrecomputed {
    long long inv01;
    long long inv02;
    long long inv12;
    long long inv03;
    long long inv13;
    long long inv23;

    long long coeff0;
    long long coeff1;
    long long coeff2;
    long long coeff3;

    long long targetModulus;
};

CrtGarnerPrecomputed buildCrtGarnerPrecomputed(long long targetModulus) {
    // 棰勮绠?Garner 鍚堝苟鎵€闇€鐨勬ā閫嗗拰绯绘暟
    CrtGarnerPrecomputed precomputed;

    precomputed.inv01 = quickPowerModCrtInt(
        CRT_MODS[0] % CRT_MODS[1],
        CRT_MODS[1] - 2,
        CRT_MODS[1]
    );

    precomputed.inv02 = quickPowerModCrtInt(
        CRT_MODS[0] % CRT_MODS[2],
        CRT_MODS[2] - 2,
        CRT_MODS[2]
    );

    precomputed.inv12 = quickPowerModCrtInt(
        CRT_MODS[1] % CRT_MODS[2],
        CRT_MODS[2] - 2,
        CRT_MODS[2]
    );

    precomputed.inv03 = quickPowerModCrtInt(
        CRT_MODS[0] % CRT_MODS[3],
        CRT_MODS[3] - 2,
        CRT_MODS[3]
    );

    precomputed.inv13 = quickPowerModCrtInt(
        CRT_MODS[1] % CRT_MODS[3],
        CRT_MODS[3] - 2,
        CRT_MODS[3]
    );

    precomputed.inv23 = quickPowerModCrtInt(
        CRT_MODS[2] % CRT_MODS[3],
        CRT_MODS[3] - 2,
        CRT_MODS[3]
    );

    precomputed.coeff0 = 1 % targetModulus;
    precomputed.coeff1 = CRT_MODS[0] % targetModulus;

    precomputed.coeff2 = static_cast<long long>(
        static_cast<__int128>(CRT_MODS[0]) *
        CRT_MODS[1] %
        targetModulus
    );

    precomputed.coeff3 = static_cast<long long>(
        static_cast<__int128>(precomputed.coeff2) *
        CRT_MODS[2] %
        targetModulus
    );

    precomputed.targetModulus = targetModulus;

    return precomputed;
}

static inline long long garnerCombineToTargetModulusFast4(
    int residue0,
    int residue1,
    int residue2,
    int residue3,
    const CrtGarnerPrecomputed &precomputed
) {
    // 4 妯℃暟 Garner 鍚堝苟锛屾渶鍚庡彧淇濈暀鐩爣妯℃暟 targetModulus 涓嬬殑缁撴灉
    long long x0 = residue0;

    long long x1 = positiveModLongLong(residue1 - x0, CRT_MODS[1]);
    x1 = x1 * precomputed.inv01 % CRT_MODS[1];

    long long x2 = positiveModLongLong(residue2 - x0, CRT_MODS[2]);
    x2 = x2 * precomputed.inv02 % CRT_MODS[2];
    x2 = positiveModLongLong(x2 - x1, CRT_MODS[2]);
    x2 = x2 * precomputed.inv12 % CRT_MODS[2];

    long long x3 = positiveModLongLong(residue3 - x0, CRT_MODS[3]);
    x3 = x3 * precomputed.inv03 % CRT_MODS[3];
    x3 = positiveModLongLong(x3 - x1, CRT_MODS[3]);
    x3 = x3 * precomputed.inv13 % CRT_MODS[3];
    x3 = positiveModLongLong(x3 - x2, CRT_MODS[3]);
    x3 = x3 * precomputed.inv23 % CRT_MODS[3];

    __int128 result = 0;

    result += static_cast<__int128>(x0) * precomputed.coeff0;
    result += static_cast<__int128>(x1) * precomputed.coeff1;
    result += static_cast<__int128>(x2) * precomputed.coeff2;
    result += static_cast<__int128>(x3) * precomputed.coeff3;

    result %= precomputed.targetModulus;

    return static_cast<long long>(result);
}

void multiplySingleCrtModForMPI(
    long long *firstPolynomial,
    long long *secondPolynomial,
    int polynomialLength,
    int transformLength,
    int modulusIndex,
    int threadCount
) {
    int modulus = CRT_MODS[modulusIndex];

#if CRT_INTRA_THREAD == 0
    for (int index = 0; index < transformLength; ++index) {
        crtWorkA[modulusIndex][index] = 0;
        crtWorkB[modulusIndex][index] = 0;
    }

    for (int index = 0; index < polynomialLength; ++index) {
        crtWorkA[modulusIndex][index] = static_cast<int>(
            firstPolynomial[index] % modulus
        );
        crtWorkB[modulusIndex][index] = static_cast<int>(
            secondPolynomial[index] % modulus
        );
    }

    numberTheoreticTransformCrtModMontgomery(
        crtWorkA[modulusIndex],
        transformLength,
        modulus,
        CRT_ROOTS[modulusIndex],
        false,
        crtStageRootsMontgomery[modulusIndex]
    );

    numberTheoreticTransformCrtModMontgomery(
        crtWorkB[modulusIndex],
        transformLength,
        modulus,
        CRT_ROOTS[modulusIndex],
        false,
        crtStageRootsMontgomery[modulusIndex]
    );

    for (int index = 0; index < transformLength; ++index) {
        crtWorkA[modulusIndex][index] = multiplyModCrtInt(
            crtWorkA[modulusIndex][index],
            crtWorkB[modulusIndex][index],
            modulus
        );
    }

    numberTheoreticTransformCrtModMontgomery(
        crtWorkA[modulusIndex],
        transformLength,
        modulus,
        CRT_ROOTS[modulusIndex],
        true,
        crtStageRootsMontgomery[modulusIndex]
    );
#else
    int effectiveThreadCount = std::max(1, threadCount);

    #pragma omp parallel num_threads(effectiveThreadCount)
    {
        int localThreadId = omp_get_thread_num();
        int localThreadCount = omp_get_num_threads();

        int beginIndex = static_cast<int>(
            1LL * transformLength * localThreadId / localThreadCount
        );

        int endIndex = static_cast<int>(
            1LL * transformLength * (localThreadId + 1) / localThreadCount
        );

        for (int index = beginIndex; index < endIndex; ++index) {
            crtWorkA[modulusIndex][index] = 0;
            crtWorkB[modulusIndex][index] = 0;

            if (index < polynomialLength) {
                crtWorkA[modulusIndex][index] = static_cast<int>(
                    firstPolynomial[index] % modulus
                );
                crtWorkB[modulusIndex][index] = static_cast<int>(
                    secondPolynomial[index] % modulus
                );
            }
        }

        #pragma omp barrier

        numberTheoreticTransformCrtModGroupedOpenMP(
            crtWorkA[modulusIndex],
            transformLength,
            modulusIndex,
            false,
            localThreadId,
            localThreadCount
        );

        numberTheoreticTransformCrtModGroupedOpenMP(
            crtWorkB[modulusIndex],
            transformLength,
            modulusIndex,
            false,
            localThreadId,
            localThreadCount
        );

        beginIndex = static_cast<int>(
            1LL * transformLength * localThreadId / localThreadCount
        );

        endIndex = static_cast<int>(
            1LL * transformLength * (localThreadId + 1) / localThreadCount
        );

        for (int index = beginIndex; index < endIndex; ++index) {
            crtWorkA[modulusIndex][index] = multiplyModCrtInt(
                crtWorkA[modulusIndex][index],
                crtWorkB[modulusIndex][index],
                modulus
            );
        }

        #pragma omp barrier

        numberTheoreticTransformCrtModGroupedOpenMP(
            crtWorkA[modulusIndex],
            transformLength,
            modulusIndex,
            true,
            localThreadId,
            localThreadCount
        );
    }
#endif
}


void multiplyByNttLargeCRTMPI(
    long long *firstPolynomial,
    long long *secondPolynomial,
    long long *resultPolynomial,
    int polynomialLength,
    long long targetModulus,
    int threadCount,
    int rank,
    int size
) {
    // 绗竴鐗堢湡姝?MPI 骞惰锛氭妸 4 涓?CRT 灏忔ā鏁版寜 rank 鍒嗛厤銆?
    // rank = 0,1,2,3 鏃跺垎鍒鐞嗕竴涓ā鏁帮紱np=2 鏃舵瘡涓?rank 澶勭悊涓や釜妯℃暟銆?
    // 姣忎釜 rank 鍙妸鑷繁璐熻矗鐨?residue 鏀惧叆 localResidues锛屾渶鍚?Reduce 鍒?rank 0銆?
    int transformLength = 1;
    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    int resultLength = 2 * polynomialLength - 1;
#if CRT_COLLECT_METHOD == 0
    std::vector<int> localResidues(CRT_MOD_COUNT * resultLength, 0);
#endif
    std::vector<int> gatheredResidues;
    if (rank == 0) {
        gatheredResidues.assign(CRT_MOD_COUNT * resultLength, 0);
    }

    for (int modulusIndex = rank;
         modulusIndex < CRT_MOD_COUNT;
         modulusIndex += size) {

        multiplySingleCrtModForMPI(
            firstPolynomial,
            secondPolynomial,
            polynomialLength,
            transformLength,
            modulusIndex,
            threadCount
        );

        int base = modulusIndex * resultLength;
#if CRT_COLLECT_METHOD == 0
        #if CRT_INTRA_THREAD == 1
        #pragma omp parallel for num_threads(threadCount) schedule(static)
        #endif
        for (int index = 0; index < resultLength; ++index) {
            localResidues[base + index] = crtWorkA[modulusIndex][index];
        }
#else
        if (rank == 0) {
            #if CRT_INTRA_THREAD == 1
            #pragma omp parallel for num_threads(threadCount) schedule(static)
            #endif
            for (int index = 0; index < resultLength; ++index) {
                gatheredResidues[base + index] = crtWorkA[modulusIndex][index];
            }
        } else {
            MPI_Send(
                &modulusIndex,
                1,
                MPI_INT,
                0,
                2001,
                MPI_COMM_WORLD
            );
            MPI_Send(
                crtWorkA[modulusIndex],
                resultLength,
                MPI_INT,
                0,
                2002,
                MPI_COMM_WORLD
            );
        }
#endif
    }

#if CRT_COLLECT_METHOD == 0
    MPI_Reduce(
        localResidues.data(),
        rank == 0 ? gatheredResidues.data() : localResidues.data(),
        CRT_MOD_COUNT * resultLength,
        MPI_INT,
        MPI_SUM,
        0,
        MPI_COMM_WORLD
    );
#else
    if (rank == 0) {
        for (int source = 1; source < size; ++source) {
            for (int expectedModulusIndex = source;
                 expectedModulusIndex < CRT_MOD_COUNT;
                 expectedModulusIndex += size) {

                int receivedModulusIndex = -1;
                MPI_Recv(
                    &receivedModulusIndex,
                    1,
                    MPI_INT,
                    source,
                    2001,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );

                MPI_Recv(
                    &gatheredResidues[receivedModulusIndex * resultLength],
                    resultLength,
                    MPI_INT,
                    source,
                    2002,
                    MPI_COMM_WORLD,
                    MPI_STATUS_IGNORE
                );
            }
        }
    }
#endif

    if (rank == 0) {
        CrtGarnerPrecomputed precomputed = buildCrtGarnerPrecomputed(targetModulus);

        #pragma omp parallel for num_threads(threadCount) schedule(static)
        for (int index = 0; index < resultLength; ++index) {
            resultPolynomial[index] = garnerCombineToTargetModulusFast4(
                gatheredResidues[0 * resultLength + index],
                gatheredResidues[1 * resultLength + index],
                gatheredResidues[2 * resultLength + index],
                gatheredResidues[3 * resultLength + index],
                precomputed
            );
        }
    }
}

//======================================================================
// 涓诲嚱鏁?
//======================================================================

int a[300000], b[300000], ab[300000];
long long aLarge[300000], bLarge[300000], abLarge[300000];

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 淇濊瘉杈撳叆鐨勬墍鏈夋ā鏁扮殑鍘熸牴鍧囦负 3, 涓旀ā鏁伴兘鑳借〃绀轰负 a 脳 4^k + 1 鐨勫舰寮?
    // 杈撳叆妯℃暟鍒嗗埆涓?7340033 104857601 469762049 浠ュ強涓€涓ぇ妯℃暟
    // 绗洓涓ぇ妯℃暟瓒呰繃浜嗘暣鍨嬭〃绀鸿寖鍥? 杩欓噷浣跨敤 CRT 鏂规澶勭悊
    // 杈撳叆鏂囦欢鍏变簲涓? 绗竴涓緭鍏ユ枃浠?n = 4, 鍏朵綑鍥涗釜鏂囦欢 n = 131072

    const int threadCount = THREAD_COUNT;
    omp_set_num_threads(threadCount);

    if (rank == 0) {
        std::cout << "MPI size = " << size
                  << ", per-process thread count = " << threadCount
                  << std::endl;
#if ORDINARY_MPI_TASK == 0
        std::cout << "input 0~3: rank 0 Pthread whole-pipeline baseline" << std::endl;
#else
        std::cout << "input 0~3: task-level MPI NTT, rank0 forward A, rank1 forward B" << std::endl;
#endif
#if CRT_COLLECT_METHOD == 0
        #if CRT_INTRA_THREAD == 0
        std::cout << "input 4: MPI CRT parallel, collect method = reduce-sum" << std::endl;
        #else
        std::cout << "input 4: MPI CRT parallel, collect method = reduce-sum, intra-rank threads = "
                  << threadCount
                  << std::endl;
        #endif
#else
        #if CRT_INTRA_THREAD == 0
        std::cout << "input 4: MPI CRT parallel, collect method = point-to-point" << std::endl;
        #else
        std::cout << "input 4: MPI CRT parallel, collect method = point-to-point, intra-rank threads = "
                  << threadCount
                  << std::endl;
        #endif
#endif
    }

    int test_begin = 0;
    int test_end = 4;

    for (int i = test_begin; i <= test_end; ++i) {
        if (i <= 3) {
#if ORDINARY_MPI_TASK == 0
            // 绗竴闃舵鍏堟妸 32-bit 杈撳叆浣滀负 rank 0 鐨?Pthread baseline 淇濈暀锛?
            // 鍚庣画鍐嶅崌绾т负 stage-wise butterfly MPI銆?
            if (rank == 0) {
                int n_, p_;

                fRead(a, b, &n_, &p_, i);
                memset(ab, 0, sizeof(ab));

                double startTime = MPI_Wtime();
                multiplyByNttPthreadWhole(a, b, ab, n_, p_, threadCount);
                double elapsedMs = (MPI_Wtime() - startTime) * 1000.0;

                fCheck(ab, n_, i);

                std::cout << "rank0 pthread baseline latency for input "
                          << i
                          << " n = "
                          << n_
                          << " p = "
                          << p_
                          << " : "
                          << elapsedMs
                          << " (ms)"
                          << std::endl;
            }
#else
            int n_ = 0;
            int p_ = 0;

            if (rank == 0) {
                fRead(a, b, &n_, &p_, i);
                memset(ab, 0, sizeof(ab));
            }

            MPI_Bcast(&n_, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&p_, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(a, n_, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(b, n_, MPI_INT, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);
            double startTime = MPI_Wtime();

            multiplyByNttTaskMPI(
                a,
                b,
                ab,
                n_,
                p_,
                threadCount,
                rank,
                size
            );

            double localElapsedMs = (MPI_Wtime() - startTime) * 1000.0;
            double globalElapsedMs = 0.0;

            MPI_Reduce(
                &localElapsedMs,
                &globalElapsedMs,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                0,
                MPI_COMM_WORLD
            );

            if (rank == 0) {
                fCheck(ab, n_, i);

                std::cout << "task-level MPI NTT latency for input "
                          << i
                          << " n = "
                          << n_
                          << " p = "
                          << p_
                          << " np = "
                          << size
                          << " : "
                          << globalElapsedMs
                          << " (ms)"
                          << std::endl;
            }
#endif
        } else {
            int n_ = 0;
            long long p_ = 0;

            if (rank == 0) {
                fReadLarge(aLarge, bLarge, &n_, &p_, i);
                memset(abLarge, 0, sizeof(abLarge));
            }

            MPI_Bcast(&n_, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&p_, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(aLarge, n_, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
            MPI_Bcast(bLarge, n_, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

            MPI_Barrier(MPI_COMM_WORLD);
            double startTime = MPI_Wtime();

            multiplyByNttLargeCRTMPI(
                aLarge,
                bLarge,
                abLarge,
                n_,
                p_,
                threadCount,
                rank,
                size
            );

            double localElapsedMs = (MPI_Wtime() - startTime) * 1000.0;
            double globalElapsedMs = 0.0;

            MPI_Reduce(
                &localElapsedMs,
                &globalElapsedMs,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                0,
                MPI_COMM_WORLD
            );

            if (rank == 0) {
                fCheckLarge(abLarge, n_, i);

                std::cout << "MPI CRT latency for input "
                          << i
                          << " n = "
                          << n_
                          << " p = "
                          << p_
                          << " np = "
                          << size
                          << " : "
                          << globalElapsedMs
                          << " (ms)"
                          << std::endl;
            }
        }
    }

    MPI_Finalize();
    return 0;
}

