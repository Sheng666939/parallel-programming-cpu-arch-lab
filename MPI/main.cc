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

#if ORDINARY_MPI_TASK != 0 && ORDINARY_MPI_TASK != 1
#error "ORDINARY_MPI_TASK must be 0 or 1"
#endif

#if CRT_COLLECT_METHOD != 0 && CRT_COLLECT_METHOD != 1
#error "CRT_COLLECT_METHOD must be 0 or 1"
#endif

// 可以自行添加需要的头文件
// 本版本只保留最终真正使用的两条主线：
// 1. input 0~3：Pthread whole-pipeline NTT
// 2. input 4：4-mod CRT + Montgomery NTT + Fast Garner

//======================================================================
// 数据输入、检查和调试输出函数
//======================================================================

void fRead(int *a, int *b, int *n, int *p, int input_id) {
    // 数据输入函数
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
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";

    std::ifstream fin;
    fin.open(strout, std::ios::in);

    for (int i = 0; i < n * 2 - 1; i++) {
        int x;
        fin >> x;

        if (x != ab[i]) {
            std::cout << "多项式乘法结果错误" << std::endl;
            return;
        }
    }

    std::cout << "多项式乘法结果正确" << std::endl;
}

void fWrite(int *ab, int n, int input_id) {
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    // 使用前请确认当前目录下已经存在 files 文件夹
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
    // 大模数输入函数
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
    // 判断大模数多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";

    std::ifstream fin;
    fin.open(strout, std::ios::in);

    for (int i = 0; i < n * 2 - 1; i++) {
        long long x;
        fin >> x;

        if (x != ab[i]) {
            std::cout << "大模数多项式乘法结果错误" << std::endl;
            return;
        }
    }

    std::cout << "大模数多项式乘法结果正确" << std::endl;
}

void fWriteLarge(long long *ab, int n, int input_id) {
    // 大模数数据输出函数, 调试 input 4 时可以使用
    // 使用前请确认当前目录下已经存在 files 文件夹
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
// 通用快速幂和 Montgomery 模乘
//======================================================================

long long quickPowerMod(long long base, long long exponent, long long modulus) {
    // 快速幂：用于计算 base^exponent mod modulus
    // 在 NTT 中用于计算单位根、逆单位根和长度的模逆元
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
    // 计算 -modulus^{-1} mod 2^32
    // Newton 迭代每轮会让正确位数翻倍
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
    // Montgomery REDC，R = 2^32
    // 返回 value * R^{-1} mod modulus
    uint32_t q = static_cast<uint32_t>(value) * modulusNegInverse;
    uint64_t reduced = (value + static_cast<uint64_t>(q) * modulus) >> 32;

    uint32_t result = static_cast<uint32_t>(reduced);

    if (result >= modulus) {
        result -= modulus;
    }

    return result;
}

//======================================================================
// 32 位模数：Pthread whole-pipeline NTT
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
    // 这里保持串行，由 thread 0 执行，避免并行 swap 带来的写冲突
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

            // stageRoots 预计算：每一层的单位根只算一次
            ctx->stageRoots[0] = 1;

            for (int offset = 1; offset < butterflySize; ++offset) {
                ctx->stageRoots[offset] = static_cast<int>(
                    1LL * ctx->stageRoots[offset - 1] * rootStep % ctx->modulus
                );
            }

            // 转成 Montgomery 形式，减少 butterfly 中的取模开销
            for (int offset = 0; offset < butterflySize; ++offset) {
                ctx->stageRootsMontgomery[offset] = static_cast<int>(
                    static_cast<uint64_t>(ctx->stageRoots[offset]) *
                    ctx->montgomeryR %
                    ctx->modulusUnsigned
                );
            }
        }

        // 等待 thread 0 完成当前 stage 的单位根预计算
        pthread_barrier_wait(&ctx->barrier);

        int blockCount = ctx->transformLength / blockSize;

        // 和 OpenMP schedule(static) 类似：每个线程处理一段连续 block
        // 这样实现简单，数据访问也比较连续
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

        // 等待当前 stage 全部完成，才能进入下一 stage
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

    // 1. 清零和输入拷贝并行化
    for (int index = beginIndex; index < endIndex; ++index) {
        ctx->firstTransformed[index] = 0;
        ctx->secondTransformed[index] = 0;

        if (index < ctx->polynomialLength) {
            ctx->firstTransformed[index] = ctx->firstInput[index] % ctx->modulus;
            ctx->secondTransformed[index] = ctx->secondInput[index] % ctx->modulus;
        }
    }

    pthread_barrier_wait(&ctx->barrier);

    // 2. first 正 NTT
    if (threadId == 0) {
        bitReverseSerial(ctx->firstTransformed, ctx->transformLength);
    }

    pthread_barrier_wait(&ctx->barrier);
    pthreadWholeRunNttStages(ctx, ctx->firstTransformed, false, threadId);

    // 3. second 正 NTT
    if (threadId == 0) {
        bitReverseSerial(ctx->secondTransformed, ctx->transformLength);
    }

    pthread_barrier_wait(&ctx->barrier);
    pthreadWholeRunNttStages(ctx, ctx->secondTransformed, false, threadId);

    // 4. 点值乘法并行化
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

    // 5. 逆 NTT
    if (threadId == 0) {
        bitReverseSerial(ctx->firstTransformed, ctx->transformLength);
    }

    pthread_barrier_wait(&ctx->barrier);
    pthreadWholeRunNttStages(ctx, ctx->firstTransformed, true, threadId);

    // 6. 拷贝真实结果并行化
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
    // 使用 static 数组避免大数组反复在栈上申请，减少栈溢出风险
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

    // 整个多项式乘法只创建一次线程组
    // 清零、两次正 NTT、点值乘法、逆 NTT、结果回拷都在同一批线程里完成
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
// 大模数 input 4：4-mod CRT + Montgomery NTT + Fast Garner
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

static int crtResidues[CRT_MOD_COUNT][300000];
static int crtWorkA[CRT_MOD_COUNT][300000];
static int crtWorkB[CRT_MOD_COUNT][300000];
static int crtStageRootsMontgomery[CRT_MOD_COUNT][300000];
static int crtInverseLengthMontgomery[CRT_MOD_COUNT];

static inline int multiplyModCrtInt(int x, int y, int modulus) {
    return static_cast<int>(1LL * x * y % modulus);
}

int quickPowerModCrtInt(int base, long long exponent, int modulus) {
    // CRT 小模数快速幂
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
    // 单个 CRT 小模数下的 NTT，供 V10-1 fallback 使用
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
    // 分组 CRT NTT：THREAD_COUNT = 8 时，每个 CRT 模数分到 2 个 OpenMP 线程
    // 注意这里的 barrier 是整个 OpenMP parallel 区域的全局 barrier
    // 因为 4 个模数的 transformLength 和 stage 结构完全相同，所以所有线程会按同样顺序到达 barrier
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

        // flat-butterfly 划分：每层 butterfly 总数恒为 transformLength / 2
        // 这样可以避免后期 stage block 数太少导致线程空转
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
    // 预计算 Garner 合并所需的模逆和系数
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
    // 4 模数 Garner 合并，最后只保留目标模数 targetModulus 下的结果
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

void multiplyByNttLargeCRTModParallel(
    long long *firstPolynomial,
    long long *secondPolynomial,
    long long *resultPolynomial,
    int polynomialLength,
    long long targetModulus,
    int threadCount
) {
    // V10-1 fallback：不同 CRT 模数之间并行
    int transformLength = 1;

    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int modulusIndex = 0;
         modulusIndex < CRT_MOD_COUNT;
         ++modulusIndex) {

        int modulus = CRT_MODS[modulusIndex];

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

        int resultLength = 2 * polynomialLength - 1;

        for (int index = 0; index < resultLength; ++index) {
            crtResidues[modulusIndex][index] = crtWorkA[modulusIndex][index];
        }
    }

    int resultLength = 2 * polynomialLength - 1;
    CrtGarnerPrecomputed precomputed = buildCrtGarnerPrecomputed(targetModulus);

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < resultLength; ++index) {
        resultPolynomial[index] = garnerCombineToTargetModulusFast4(
            crtResidues[0][index],
            crtResidues[1][index],
            crtResidues[2][index],
            crtResidues[3][index],
            precomputed
        );
    }
}

void multiplyByNttLargeCRTGrouped(
    long long *firstPolynomial,
    long long *secondPolynomial,
    long long *resultPolynomial,
    int polynomialLength,
    long long targetModulus,
    int threadCount
) {
    // V10-2：4-mod grouped CRT
    // 当 THREAD_COUNT = 8 时，每个 CRT 模数分配 2 个线程
    int transformLength = 1;

    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    int resultLength = 2 * polynomialLength - 1;
    int localThreadCount = threadCount / CRT_MOD_COUNT;

    if (localThreadCount < 2) {
        // 线程数太少时，不强行分组，直接退回到模数级并行版本
        multiplyByNttLargeCRTModParallel(
            firstPolynomial,
            secondPolynomial,
            resultPolynomial,
            polynomialLength,
            targetModulus,
            threadCount
        );

        return;
    }

    int activeThreadCount = localThreadCount * CRT_MOD_COUNT;

    #pragma omp parallel num_threads(activeThreadCount)
    {
        int globalThreadId = omp_get_thread_num();
        int modulusIndex = globalThreadId / localThreadCount;
        int localThreadId = globalThreadId - modulusIndex * localThreadCount;

        int modulus = CRT_MODS[modulusIndex];

        int beginIndex = static_cast<int>(
            1LL * transformLength * localThreadId / localThreadCount
        );

        int endIndex = static_cast<int>(
            1LL * transformLength * (localThreadId + 1) / localThreadCount
        );

        // 清零和输入映射到当前 CRT 模数
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

        // 两个输入多项式分别做正 NTT
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

        // 点值表示下逐点相乘
        for (int index = beginIndex; index < endIndex; ++index) {
            crtWorkA[modulusIndex][index] = multiplyModCrtInt(
                crtWorkA[modulusIndex][index],
                crtWorkB[modulusIndex][index],
                modulus
            );
        }

        #pragma omp barrier

        // 逆 NTT 回到系数表示
        numberTheoreticTransformCrtModGroupedOpenMP(
            crtWorkA[modulusIndex],
            transformLength,
            modulusIndex,
            true,
            localThreadId,
            localThreadCount
        );

        int beginResult = static_cast<int>(
            1LL * resultLength * localThreadId / localThreadCount
        );

        int endResult = static_cast<int>(
            1LL * resultLength * (localThreadId + 1) / localThreadCount
        );

        for (int index = beginResult; index < endResult; ++index) {
            crtResidues[modulusIndex][index] = crtWorkA[modulusIndex][index];
        }

        #pragma omp barrier
    }

    CrtGarnerPrecomputed precomputed = buildCrtGarnerPrecomputed(targetModulus);

    // Garner 合并阶段继续用 OpenMP 并行
    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < resultLength; ++index) {
        resultPolynomial[index] = garnerCombineToTargetModulusFast4(
            crtResidues[0][index],
            crtResidues[1][index],
            crtResidues[2][index],
            crtResidues[3][index],
            precomputed
        );
    }
}

void multiplyByNttLargeCRT(
    long long *firstPolynomial,
    long long *secondPolynomial,
    long long *resultPolynomial,
    int polynomialLength,
    long long targetModulus,
    int threadCount
) {
    // 大模数最终版本入口
    multiplyByNttLargeCRTGrouped(
        firstPolynomial,
        secondPolynomial,
        resultPolynomial,
        polynomialLength,
        targetModulus,
        threadCount
    );
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
    // 第一版真正 MPI 并行：把 4 个 CRT 小模数按 rank 分配。
    // rank = 0,1,2,3 时分别处理一个模数；np=2 时每个 rank 处理两个模数。
    // 每个 rank 只把自己负责的 residue 放入 localResidues，最后 Reduce 到 rank 0。
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

        int modulus = CRT_MODS[modulusIndex];

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

        int base = modulusIndex * resultLength;
#if CRT_COLLECT_METHOD == 0
        for (int index = 0; index < resultLength; ++index) {
            localResidues[base + index] = crtWorkA[modulusIndex][index];
        }
#else
        if (rank == 0) {
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
// 主函数
//======================================================================

int a[300000], b[300000], ab[300000];
long long aLarge[300000], bLarge[300000], abLarge[300000];

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank = 0;
    int size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a × 4^k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 以及一个大模数
    // 第四个大模数超过了整型表示范围, 这里使用 CRT 方案处理
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件 n = 131072

    const int threadCount = THREAD_COUNT;

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
        std::cout << "input 4: MPI CRT parallel, collect method = reduce-sum" << std::endl;
#else
        std::cout << "input 4: MPI CRT parallel, collect method = point-to-point" << std::endl;
#endif
    }

    int test_begin = 0;
    int test_end = 4;

    for (int i = test_begin; i <= test_end; ++i) {
        if (i <= 3) {
#if ORDINARY_MPI_TASK == 0
            // 第一阶段先把 32-bit 输入作为 rank 0 的 Pthread baseline 保留，
            // 后续再升级为 stage-wise butterfly MPI。
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
