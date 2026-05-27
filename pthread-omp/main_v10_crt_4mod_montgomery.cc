#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sys/time.h>
#include <omp.h>
#include <algorithm>
#include <cstdint>
#include <pthread.h>
#include <vector>
#ifndef THREAD_COUNT
#define THREAD_COUNT 8
#endif
// 可以自行添加需要的头文件

void fRead(int *a, int *b, int *n, int *p, int input_id){
    // 数据输入函数
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";
    char data_path[strin.size() + 1];
    std::copy(strin.begin(), strin.end(), data_path);
    data_path[strin.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    fin>>*n>>*p;
    for (int i = 0; i < *n; i++){
        fin>>a[i];
    }
    for (int i = 0; i < *n; i++){   
        fin>>b[i];
    }
}

void fCheck(int *ab, int n, int input_id){
    // 判断多项式乘法结果是否正确
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char data_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), data_path);
    data_path[strout.size()] = '\0';
    std::ifstream fin;
    fin.open(data_path, std::ios::in);
    for (int i = 0; i < n * 2 - 1; i++){
        int x;
        fin>>x;
        if(x != ab[i]){
            std::cout<<"多项式乘法结果错误"<<std::endl;
            return;
        }
    }
    std::cout<<"多项式乘法结果正确"<<std::endl;
    return;
}

void fWrite(int *ab, int n, int input_id){
    // 数据输出函数, 可以用来输出最终结果, 也可用于调试时输出中间数组
    std::string str1 = "files/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";
    char output_path[strout.size() + 1];
    std::copy(strout.begin(), strout.end(), output_path);
    output_path[strout.size()] = '\0';
    std::ofstream fout;
    fout.open(output_path, std::ios::out);
    for (int i = 0; i < n * 2 - 1; i++){
        fout<<ab[i]<<'\n';
    }
}

/*void poly_multiply(int *a, int *b, int *ab, int n, int p){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j){
            ab[i+j]=(1LL * a[i] * b[j] % p + ab[i+j]) % p;
        }
    }
}
*///原来的朴素算法
long long quickPowerMod(long long base, long long exponent, long long modulus) {
    // 快速幂：用于计算 base^exponent mod modulus。
    // 在 NTT 中，我们会用它计算单位根、逆单位根，以及 len 的模逆元。
    long long result = 1;

    while (exponent > 0) {
        if (exponent & 1) {
            result = result * base % modulus;
        }

        base = base * base % modulus;
        exponent >>= 1;
    }

    return result;
}
static inline uint32_t computeMontgomeryNegInverse32(uint32_t modulus) {
    // Compute -modulus^{-1} modulo 2^32.
    // Newton iteration doubles the number of correct bits each round.
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
    // Standard Montgomery REDC with R = 2^32:
    // result = value * R^{-1} mod modulus.
    uint32_t q = static_cast<uint32_t>(value) * modulusNegInverse;
    uint64_t reduced = (value + static_cast<uint64_t>(q) * modulus) >> 32;

    uint32_t result = static_cast<uint32_t>(reduced);
    if (result >= modulus) {
        result -= modulus;
    }
    return result;
}

void numberTheoreticTransformScalar(
    int *data,
    int transformLength,
    int modulus,
    bool inverseTransform
) {
    // bit-reversal permutation
    for (int currentIndex = 1, reversedIndex = 0; currentIndex < transformLength; ++currentIndex) {
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

    const int primitiveRoot = 3;

    const uint32_t modulusUnsigned = static_cast<uint32_t>(modulus);
    const uint32_t montgomeryNegInverse = computeMontgomeryNegInverse32(modulusUnsigned);
    const uint32_t montgomeryR = static_cast<uint32_t>((1ULL << 32) % modulusUnsigned);

    static int stageRoots[300000];
    static int stageRootsMontgomery[300000];

    for (int butterflySize = 1; butterflySize < transformLength; butterflySize <<= 1) {
        int blockSize = butterflySize << 1;

        long long rootStep = quickPowerMod(
            primitiveRoot,
            (modulus - 1) / blockSize,
            modulus
        );

        if (inverseTransform) {
            rootStep = quickPowerMod(rootStep, modulus - 2, modulus);
        }

        // stageRoots 预计算：标量优化，不是 SIMD
        stageRoots[0] = 1;
        for (int offset = 1; offset < butterflySize; ++offset) {
            stageRoots[offset] = static_cast<int>(
                1LL * stageRoots[offset - 1] * rootStep % modulus
            );
        }

        for (int offset = 0; offset < butterflySize; ++offset) {
            stageRootsMontgomery[offset] = static_cast<int>(
                static_cast<uint64_t>(stageRoots[offset]) * montgomeryR % modulusUnsigned
            );
        }

        // 纯标量蝶形运算
        for (int blockStart = 0; blockStart < transformLength; blockStart += blockSize) {
            for (int offset = 0; offset < butterflySize; ++offset) {
                int upperValue = data[blockStart + offset];

                uint64_t product =
                    static_cast<uint64_t>(
                        static_cast<uint32_t>(data[blockStart + offset + butterflySize])
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

                data[blockStart + offset] = addedValue;
                data[blockStart + offset + butterflySize] = subtractedValue;
            }
        }
    }

    if (inverseTransform) {
        long long inverseLength = quickPowerMod(transformLength, modulus - 2, modulus);

        for (int index = 0; index < transformLength; ++index) {
            data[index] = static_cast<int>(
                1LL * data[index] * inverseLength % modulus
            );
        }
    }
}
void numberTheoreticTransformOpenMP(
    int *data,
    int transformLength,
    int modulus,
    bool inverseTransform,
    int threadCount
) {
    // bit-reversal permutation 暂时保持串行。
    // 该部分不是本版本优化重点，避免并行 swap 产生复杂一致性问题。
    for (int currentIndex = 1, reversedIndex = 0; currentIndex < transformLength; ++currentIndex) {
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

    const int primitiveRoot = 3;

    const uint32_t modulusUnsigned = static_cast<uint32_t>(modulus);
    const uint32_t montgomeryNegInverse = computeMontgomeryNegInverse32(modulusUnsigned);
    const uint32_t montgomeryR = static_cast<uint32_t>((1ULL << 32) % modulusUnsigned);

    static int stageRoots[300000];
    static int stageRootsMontgomery[300000];

    // 在整个 NTT stage 循环外创建一次 OpenMP 线程组。
    // 这样避免每一层 butterflySize 都重复 fork-join。
    #pragma omp parallel num_threads(threadCount)
    {
        for (int butterflySize = 1; butterflySize < transformLength; butterflySize <<= 1) {
            int blockSize = butterflySize << 1;

            // 每一层的旋转因子只由一个线程预计算。
            #pragma omp single
            {
                long long rootStep = quickPowerMod(
                    primitiveRoot,
                    (modulus - 1) / blockSize,
                    modulus
                );

                if (inverseTransform) {
                    rootStep = quickPowerMod(rootStep, modulus - 2, modulus);
                }

                stageRoots[0] = 1;
                for (int offset = 1; offset < butterflySize; ++offset) {
                    stageRoots[offset] = static_cast<int>(
                        1LL * stageRoots[offset - 1] * rootStep % modulus
                    );
                }

                for (int offset = 0; offset < butterflySize; ++offset) {
                    stageRootsMontgomery[offset] = static_cast<int>(
                        static_cast<uint64_t>(stageRoots[offset]) * montgomeryR % modulusUnsigned
                    );
                }
            }
            // single 末尾默认有隐式 barrier。
            // 因此所有线程都会等 stageRootsMontgomery 准备好后再进入 butterfly。

            int blockCount = transformLength / blockSize;

            // 每个 block 的数据区间互不重叠，所以 block 级并行是安全的。
            #pragma omp for schedule(static)
            for (int blockIndex = 0; blockIndex < blockCount; ++blockIndex) {
                int blockStart = blockIndex * blockSize;

                for (int offset = 0; offset < butterflySize; ++offset) {
                    int upperValue = data[blockStart + offset];

                    uint64_t product =
                        static_cast<uint64_t>(
                            static_cast<uint32_t>(data[blockStart + offset + butterflySize])
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

                    data[blockStart + offset] = addedValue;
                    data[blockStart + offset + butterflySize] = subtractedValue;
                }
            }
            // omp for 末尾默认有隐式 barrier。
            // 这正好保证当前 stage 完成后才能进入下一 stage。
        }
    }

    if (inverseTransform) {
        long long inverseLength = quickPowerMod(transformLength, modulus - 2, modulus);

        #pragma omp parallel for num_threads(threadCount) schedule(static)
        for (int index = 0; index < transformLength; ++index) {
            data[index] = static_cast<int>(
                1LL * data[index] * inverseLength % modulus
            );
        }
    }
}
struct PthreadNttContext {
    int *data;
    int transformLength;
    int modulus;
    bool inverseTransform;
    int threadCount;

    uint32_t modulusUnsigned;
    uint32_t montgomeryNegInverse;
    uint32_t montgomeryR;

    int *stageRoots;
    int *stageRootsMontgomery;

    long long inverseLength;

    pthread_barrier_t barrier;
};

struct PthreadNttWorkerParam {
    int threadId;
    PthreadNttContext *context;
};

void *pthreadNttWorker(void *arg) {
    PthreadNttWorkerParam *param = static_cast<PthreadNttWorkerParam *>(arg);
    int threadId = param->threadId;
    PthreadNttContext *ctx = param->context;

    const int primitiveRoot = 3;

    for (int butterflySize = 1; butterflySize < ctx->transformLength; butterflySize <<= 1) {
        int blockSize = butterflySize << 1;

        if (threadId == 0) {
            long long rootStep = quickPowerMod(
                primitiveRoot,
                (ctx->modulus - 1) / blockSize,
                ctx->modulus
            );

            if (ctx->inverseTransform) {
                rootStep = quickPowerMod(rootStep, ctx->modulus - 2, ctx->modulus);
            }

            ctx->stageRoots[0] = 1;
            for (int offset = 1; offset < butterflySize; ++offset) {
                ctx->stageRoots[offset] = static_cast<int>(
                    1LL * ctx->stageRoots[offset - 1] * rootStep % ctx->modulus
                );
            }

            for (int offset = 0; offset < butterflySize; ++offset) {
                ctx->stageRootsMontgomery[offset] = static_cast<int>(
                    static_cast<uint64_t>(ctx->stageRoots[offset]) *
                    ctx->montgomeryR %
                    ctx->modulusUnsigned
                );
            }
        }

        // 等待 thread 0 完成当前 stage 的旋转因子预计算。
        pthread_barrier_wait(&ctx->barrier);

        int totalButterflies = ctx->transformLength >> 1;

        // contiguous flat-butterfly 任务划分。
        // 每个线程处理一段连续 butterfly 区间，减少最内层的整数除法开销，
        // 同时保持每层总任务数为 transformLength / 2，从而避免后期 stage 负载不均。
        int beginButterfly = static_cast<int>(
            1LL * totalButterflies * threadId / ctx->threadCount
        );
        int endButterfly = static_cast<int>(
            1LL * totalButterflies * (threadId + 1) / ctx->threadCount
        );

        int blockIndex = beginButterfly / butterflySize;
        int offset = beginButterfly - blockIndex * butterflySize;
        int blockStart = blockIndex * blockSize;

        for (int butterflyId = beginButterfly;
            butterflyId < endButterfly;
            ++butterflyId) {

            int upperIndex = blockStart + offset;
            int lowerIndex = upperIndex + butterflySize;

            int upperValue = ctx->data[upperIndex];

            uint64_t product =
                static_cast<uint64_t>(
                    static_cast<uint32_t>(ctx->data[lowerIndex])
                ) *
                static_cast<uint32_t>(ctx->stageRootsMontgomery[offset]);

            int lowerValue = static_cast<int>(
                montgomeryReduce32(product, ctx->modulusUnsigned, ctx->montgomeryNegInverse)
            );

            int addedValue = upperValue + lowerValue;
            if (addedValue >= ctx->modulus) {
                addedValue -= ctx->modulus;
            }

            int subtractedValue = upperValue - lowerValue;
            if (subtractedValue < 0) {
                subtractedValue += ctx->modulus;
            }

            ctx->data[upperIndex] = addedValue;
            ctx->data[lowerIndex] = subtractedValue;

            ++offset;
            if (offset == butterflySize) {
                offset = 0;
                blockStart += blockSize;
            }
        }

        // 等待当前 stage 所有 block 完成，才能进入下一 stage。
        pthread_barrier_wait(&ctx->barrier);
    }

    if (ctx->inverseTransform) {
        if (threadId == 0) {
            ctx->inverseLength = quickPowerMod(
                ctx->transformLength,
                ctx->modulus - 2,
                ctx->modulus
            );
        }

        pthread_barrier_wait(&ctx->barrier);

        for (int index = threadId; index < ctx->transformLength; index += ctx->threadCount) {
            ctx->data[index] = static_cast<int>(
                1LL * ctx->data[index] * ctx->inverseLength % ctx->modulus
            );
        }

        pthread_barrier_wait(&ctx->barrier);
    }

    return nullptr;
}

void numberTheoreticTransformPthread(
    int *data,
    int transformLength,
    int modulus,
    bool inverseTransform,
    int threadCount
) {
    // bit-reversal permutation 先保持串行。
    for (int currentIndex = 1, reversedIndex = 0; currentIndex < transformLength; ++currentIndex) {
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

    static int stageRoots[300000];
    static int stageRootsMontgomery[300000];

    PthreadNttContext context;
    context.data = data;
    context.transformLength = transformLength;
    context.modulus = modulus;
    context.inverseTransform = inverseTransform;
    context.threadCount = threadCount;

    context.modulusUnsigned = static_cast<uint32_t>(modulus);
    context.montgomeryNegInverse = computeMontgomeryNegInverse32(context.modulusUnsigned);
    context.montgomeryR = static_cast<uint32_t>((1ULL << 32) % context.modulusUnsigned);

    context.stageRoots = stageRoots;
    context.stageRootsMontgomery = stageRootsMontgomery;
    context.inverseLength = 1;

    pthread_barrier_init(&context.barrier, nullptr, threadCount);

    std::vector<pthread_t> handles(threadCount);
    std::vector<PthreadNttWorkerParam> params(threadCount);

    for (int threadId = 0; threadId < threadCount; ++threadId) {
        params[threadId].threadId = threadId;
        params[threadId].context = &context;

        pthread_create(
            &handles[threadId],
            nullptr,
            pthreadNttWorker,
            &params[threadId]
        );
    }

    for (int threadId = 0; threadId < threadCount; ++threadId) {
        pthread_join(handles[threadId], nullptr);
    }

    pthread_barrier_destroy(&context.barrier);
}

void multiplyByNtt(int *firstPolynomial, int *secondPolynomial, int *resultPolynomial, int polynomialLength, int modulus) {
    // 使用 static 数组避免大数组反复在栈上申请，减少栈溢出风险。
    // transformLength 最大需要覆盖 2 * polynomialLength - 1 项。
    static int firstTransformed[300000];
    static int secondTransformed[300000];
    int transformLength = 1;
    // 多项式乘法结果长度为 2n - 1。
    // NTT 长度必须是 2 的幂，所以这里向上补到最近的 2 的幂。
    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }
    // 每次测试前都要清空工作数组。
    // 否则上一次测试残留的数据会污染本次 NTT 结果。
    for (int index = 0; index < transformLength; ++index) {
        firstTransformed[index] = 0;
        secondTransformed[index] = 0;
    }
    // 把两个输入多项式复制到工作数组中。
    // 后面的高位保持为 0，相当于补零。
    for (int index = 0; index < polynomialLength; ++index) {
        firstTransformed[index] = firstPolynomial[index] % modulus;
        secondTransformed[index] = secondPolynomial[index] % modulus;
    }
    // 先把两个多项式从系数表示转换到点值表示。
    numberTheoreticTransformScalar(firstTransformed, transformLength, modulus, false);
    numberTheoreticTransformScalar(secondTransformed, transformLength, modulus, false);
    // 点值表示下，多项式乘法变成逐点相乘。
    for (int index = 0; index < transformLength; ++index) {
        firstTransformed[index] =
            static_cast<int>(1LL * firstTransformed[index] * secondTransformed[index] % modulus);
    }
    // 逆 NTT 把点值表示转换回系数表示。
    numberTheoreticTransformScalar(firstTransformed, transformLength, modulus, true);
    // 只拷贝真实结果长度 2n - 1。
    // transformLength 中多出来的部分只是补零空间，不属于最终答案。
    for (int index = 0; index < 2 * polynomialLength - 1; ++index) {
        resultPolynomial[index] = firstTransformed[index];
    }
}
void multiplyByNttOpenMP(
    int *firstPolynomial,
    int *secondPolynomial,
    int *resultPolynomial,
    int polynomialLength,
    int modulus,
    int threadCount
) {
    static int firstTransformed[300000];
    static int secondTransformed[300000];

    int transformLength = 1;
    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < transformLength; ++index) {
        firstTransformed[index] = 0;
        secondTransformed[index] = 0;
    }

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < polynomialLength; ++index) {
        firstTransformed[index] = firstPolynomial[index] % modulus;
        secondTransformed[index] = secondPolynomial[index] % modulus;
    }

    numberTheoreticTransformOpenMP(firstTransformed, transformLength, modulus, false, threadCount);
    numberTheoreticTransformOpenMP(secondTransformed, transformLength, modulus, false, threadCount);

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < transformLength; ++index) {
        firstTransformed[index] =
            static_cast<int>(1LL * firstTransformed[index] * secondTransformed[index] % modulus);
    }

    numberTheoreticTransformOpenMP(firstTransformed, transformLength, modulus, true, threadCount);

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < 2 * polynomialLength - 1; ++index) {
        resultPolynomial[index] = firstTransformed[index];
    }
}

void multiplyByNttPthread(
    int *firstPolynomial,
    int *secondPolynomial,
    int *resultPolynomial,
    int polynomialLength,
    int modulus,
    int threadCount
) {
    static int firstTransformed[300000];
    static int secondTransformed[300000];

    int transformLength = 1;
    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    // 这里先保持串行，保证 Pthread 版本不混入 OpenMP。
    for (int index = 0; index < transformLength; ++index) {
        firstTransformed[index] = 0;
        secondTransformed[index] = 0;
    }

    for (int index = 0; index < polynomialLength; ++index) {
        firstTransformed[index] = firstPolynomial[index] % modulus;
        secondTransformed[index] = secondPolynomial[index] % modulus;
    }

    numberTheoreticTransformPthread(firstTransformed, transformLength, modulus, false, threadCount);
    numberTheoreticTransformPthread(secondTransformed, transformLength, modulus, false, threadCount);

    // 点值乘法先保持串行，后面可以再做 Pthread 化作为优化版本。
    for (int index = 0; index < transformLength; ++index) {
        firstTransformed[index] =
            static_cast<int>(1LL * firstTransformed[index] * secondTransformed[index] % modulus);
    }

    numberTheoreticTransformPthread(firstTransformed, transformLength, modulus, true, threadCount);

    for (int index = 0; index < 2 * polynomialLength - 1; ++index) {
        resultPolynomial[index] = firstTransformed[index];
    }
}
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
    for (int currentIndex = 1, reversedIndex = 0; currentIndex < transformLength; ++currentIndex) {
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

    for (int butterflySize = 1; butterflySize < ctx->transformLength; butterflySize <<= 1) {
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

            ctx->stageRoots[0] = 1;
            for (int offset = 1; offset < butterflySize; ++offset) {
                ctx->stageRoots[offset] = static_cast<int>(
                    1LL * ctx->stageRoots[offset - 1] * rootStep % ctx->modulus
                );
            }

            for (int offset = 0; offset < butterflySize; ++offset) {
                ctx->stageRootsMontgomery[offset] = static_cast<int>(
                    static_cast<uint64_t>(ctx->stageRoots[offset]) *
                    ctx->montgomeryR %
                    ctx->modulusUnsigned
                );
            }
        }

        pthread_barrier_wait(&ctx->barrier);

        int blockCount = ctx->transformLength / blockSize;

        // 和 OpenMP schedule(static) 类似：每个线程拿一段连续 block。
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
                    montgomeryReduce32(product, ctx->modulusUnsigned, ctx->montgomeryNegInverse)
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
        ctx->firstTransformed[index] =
            static_cast<int>(
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
void fReadLarge(long long *a, long long *b, int *n, long long *p, int input_id) {
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strin = str1 + str2 + ".in";

    std::ifstream fin;
    fin.open(strin, std::ios::in);

    fin >> *n >> *p;

    for (int i = 0; i < *n; ++i) {
        fin >> a[i];
    }

    for (int i = 0; i < *n; ++i) {
        fin >> b[i];
    }
}

void fCheckLarge(long long *ab, int n, int input_id) {
    std::string str1 = "/nttdata/";
    std::string str2 = std::to_string(input_id);
    std::string strout = str1 + str2 + ".out";

    std::ifstream fin;
    fin.open(strout, std::ios::in);

    for (int i = 0; i < n * 2 - 1; ++i) {
        long long x;
        fin >> x;

        if (x != ab[i]) {
            std::cout << "大模数多项式乘法结果错误" << std::endl;
            return;
        }
    }

    std::cout << "大模数多项式乘法结果正确" << std::endl;
}

static inline long long multiplyModLarge(
    long long x,
    long long y,
    long long modulus
) {
    return static_cast<long long>(
        static_cast<__int128>(x) * y % modulus
    );
}

long long quickPowerModLarge(
    long long base,
    long long exponent,
    long long modulus
) {
    long long result = 1 % modulus;
    base %= modulus;

    while (exponent > 0) {
        if (exponent & 1LL) {
            result = multiplyModLarge(result, base, modulus);
        }

        base = multiplyModLarge(base, base, modulus);
        exponent >>= 1LL;
    }

    return result;
}

void numberTheoreticTransformLargeScalar(
    long long *data,
    int transformLength,
    long long modulus,
    bool inverseTransform
) {
    // bit-reversal permutation
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

    const long long primitiveRoot = 3;

    static long long stageRootsLarge[300000];

    for (int butterflySize = 1;
         butterflySize < transformLength;
         butterflySize <<= 1) {

        int blockSize = butterflySize << 1;

        long long rootStep = quickPowerModLarge(
            primitiveRoot,
            (modulus - 1) / blockSize,
            modulus
        );

        if (inverseTransform) {
            rootStep = quickPowerModLarge(rootStep, modulus - 2, modulus);
        }

        stageRootsLarge[0] = 1;

        for (int offset = 1; offset < butterflySize; ++offset) {
            stageRootsLarge[offset] = multiplyModLarge(
                stageRootsLarge[offset - 1],
                rootStep,
                modulus
            );
        }

        for (int blockStart = 0;
             blockStart < transformLength;
             blockStart += blockSize) {

            for (int offset = 0; offset < butterflySize; ++offset) {
                int upperIndex = blockStart + offset;
                int lowerIndex = upperIndex + butterflySize;

                long long upperValue = data[upperIndex];

                long long lowerValue = multiplyModLarge(
                    data[lowerIndex],
                    stageRootsLarge[offset],
                    modulus
                );

                long long addedValue = upperValue + lowerValue;
                if (addedValue >= modulus) {
                    addedValue -= modulus;
                }

                long long subtractedValue = upperValue - lowerValue;
                if (subtractedValue < 0) {
                    subtractedValue += modulus;
                }

                data[upperIndex] = addedValue;
                data[lowerIndex] = subtractedValue;
            }
        }
    }

    if (inverseTransform) {
        long long inverseLength = quickPowerModLarge(
            transformLength,
            modulus - 2,
            modulus
        );

        for (int index = 0; index < transformLength; ++index) {
            data[index] = multiplyModLarge(
                data[index],
                inverseLength,
                modulus
            );
        }
    }
}

void multiplyByNttLargeScalar(
    long long *firstPolynomial,
    long long *secondPolynomial,
    long long *resultPolynomial,
    int polynomialLength,
    long long modulus
) {
    static long long firstTransformedLarge[300000];
    static long long secondTransformedLarge[300000];

    int transformLength = 1;

    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    for (int index = 0; index < transformLength; ++index) {
        firstTransformedLarge[index] = 0;
        secondTransformedLarge[index] = 0;
    }

    for (int index = 0; index < polynomialLength; ++index) {
        firstTransformedLarge[index] = firstPolynomial[index] % modulus;
        secondTransformedLarge[index] = secondPolynomial[index] % modulus;
    }

    numberTheoreticTransformLargeScalar(
        firstTransformedLarge,
        transformLength,
        modulus,
        false
    );

    numberTheoreticTransformLargeScalar(
        secondTransformedLarge,
        transformLength,
        modulus,
        false
    );

    for (int index = 0; index < transformLength; ++index) {
        firstTransformedLarge[index] = multiplyModLarge(
            firstTransformedLarge[index],
            secondTransformedLarge[index],
            modulus
        );
    }

    numberTheoreticTransformLargeScalar(
        firstTransformedLarge,
        transformLength,
        modulus,
        true
    );

    for (int index = 0; index < 2 * polynomialLength - 1; ++index) {
        resultPolynomial[index] = firstTransformedLarge[index];
    }
}
void numberTheoreticTransformLargeOpenMP(
    long long *data,
    int transformLength,
    long long modulus,
    bool inverseTransform,
    int threadCount
) {
    // bit-reversal 暂时保持串行，避免并行 swap 的一致性问题。
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

    const long long primitiveRoot = 3;

    static long long stageRootsLargeOpenMP[300000];

    #pragma omp parallel num_threads(threadCount)
    {
        for (int butterflySize = 1;
             butterflySize < transformLength;
             butterflySize <<= 1) {

            int blockSize = butterflySize << 1;

            #pragma omp single
            {
                long long rootStep = quickPowerModLarge(
                    primitiveRoot,
                    (modulus - 1) / blockSize,
                    modulus
                );

                if (inverseTransform) {
                    rootStep = quickPowerModLarge(rootStep, modulus - 2, modulus);
                }

                stageRootsLargeOpenMP[0] = 1;

                for (int offset = 1; offset < butterflySize; ++offset) {
                    stageRootsLargeOpenMP[offset] = multiplyModLarge(
                        stageRootsLargeOpenMP[offset - 1],
                        rootStep,
                        modulus
                    );
                }
            }

            // single 结束有隐式 barrier，保证 roots 已经准备好。
            int blockCount = transformLength / blockSize;

            #pragma omp for schedule(static)
            for (int blockIndex = 0; blockIndex < blockCount; ++blockIndex) {
                int blockStart = blockIndex * blockSize;

                for (int offset = 0; offset < butterflySize; ++offset) {
                    int upperIndex = blockStart + offset;
                    int lowerIndex = upperIndex + butterflySize;

                    long long upperValue = data[upperIndex];

                    long long lowerValue = multiplyModLarge(
                        data[lowerIndex],
                        stageRootsLargeOpenMP[offset],
                        modulus
                    );

                    long long addedValue = upperValue + lowerValue;
                    if (addedValue >= modulus) {
                        addedValue -= modulus;
                    }

                    long long subtractedValue = upperValue - lowerValue;
                    if (subtractedValue < 0) {
                        subtractedValue += modulus;
                    }

                    data[upperIndex] = addedValue;
                    data[lowerIndex] = subtractedValue;
                }
            }

            // omp for 结束有隐式 barrier，保证当前 stage 完成。
        }
    }

    if (inverseTransform) {
        long long inverseLength = quickPowerModLarge(
            transformLength,
            modulus - 2,
            modulus
        );

        #pragma omp parallel for num_threads(threadCount) schedule(static)
        for (int index = 0; index < transformLength; ++index) {
            data[index] = multiplyModLarge(
                data[index],
                inverseLength,
                modulus
            );
        }
    }
}

void multiplyByNttLargeOpenMP(
    long long *firstPolynomial,
    long long *secondPolynomial,
    long long *resultPolynomial,
    int polynomialLength,
    long long modulus,
    int threadCount
) {
    static long long firstTransformedLarge[300000];
    static long long secondTransformedLarge[300000];

    int transformLength = 1;

    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < transformLength; ++index) {
        firstTransformedLarge[index] = 0;
        secondTransformedLarge[index] = 0;
    }

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < polynomialLength; ++index) {
        firstTransformedLarge[index] = firstPolynomial[index] % modulus;
        secondTransformedLarge[index] = secondPolynomial[index] % modulus;
    }

    numberTheoreticTransformLargeOpenMP(
        firstTransformedLarge,
        transformLength,
        modulus,
        false,
        threadCount
    );

    numberTheoreticTransformLargeOpenMP(
        secondTransformedLarge,
        transformLength,
        modulus,
        false,
        threadCount
    );

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < transformLength; ++index) {
        firstTransformedLarge[index] = multiplyModLarge(
            firstTransformedLarge[index],
            secondTransformedLarge[index],
            modulus
        );
    }

    numberTheoreticTransformLargeOpenMP(
        firstTransformedLarge,
        transformLength,
        modulus,
        true,
        threadCount
    );

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < 2 * polynomialLength - 1; ++index) {
        resultPolynomial[index] = firstTransformedLarge[index];
    }
}

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

static inline int multiplyModCrtInt(int x, int y, int modulus) {
    return static_cast<int>(1LL * x * y % modulus);
}

int quickPowerModCrtInt(int base, long long exponent, int modulus) {
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

void numberTheoreticTransformCrtModMontgomery(
    int *data,
    int transformLength,
    int modulus,
    int primitiveRoot,
    bool inverseTransform,
    int *stageRootsMontgomery
) {
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
    const uint32_t montgomeryNegInverse =
        computeMontgomeryNegInverse32(modulusUnsigned);
    const uint32_t montgomeryR =
        static_cast<uint32_t>((1ULL << 32) % modulusUnsigned);

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

            currentRoot = multiplyModCrtInt(
                currentRoot,
                rootStep,
                modulus
            );
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
                    montgomeryReduce32(
                        product,
                        modulusUnsigned,
                        montgomeryNegInverse
                    )
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
                montgomeryReduce32(
                    product,
                    modulusUnsigned,
                    montgomeryNegInverse
                )
            );
        }
    }
}

long long garnerCombineToTargetModulus(
    const int *residues,
    const long long inverseTable[CRT_MOD_COUNT][CRT_MOD_COUNT],
    long long targetModulus
) {
    long long mixedRadix[CRT_MOD_COUNT];

    for (int i = 0; i < CRT_MOD_COUNT; ++i) {
        long long currentValue = residues[i];

        for (int j = 0; j < i; ++j) {
            currentValue -= mixedRadix[j];
            currentValue %= CRT_MODS[i];

            if (currentValue < 0) {
                currentValue += CRT_MODS[i];
            }

            currentValue =
                currentValue *
                inverseTable[j][i] %
                CRT_MODS[i];
        }

        mixedRadix[i] = currentValue;
    }

    long long resultModuloTarget = 0;
    long long multiplierModuloTarget = 1 % targetModulus;

    for (int i = 0; i < CRT_MOD_COUNT; ++i) {
        resultModuloTarget = static_cast<long long>(
            (
                static_cast<__int128>(mixedRadix[i]) *
                multiplierModuloTarget +
                resultModuloTarget
            ) % targetModulus
        );

        multiplierModuloTarget = static_cast<long long>(
            static_cast<__int128>(multiplierModuloTarget) *
            CRT_MODS[i] %
            targetModulus
        );
    }

    return resultModuloTarget;
}

void multiplyByNttLargeCRT(
    long long *firstPolynomial,
    long long *secondPolynomial,
    long long *resultPolynomial,
    int polynomialLength,
    long long targetModulus,
    int threadCount
) {
    int transformLength = 1;

    while (transformLength < 2 * polynomialLength - 1) {
        transformLength <<= 1;
    }

    // V10-1:
    // 1. Remove the smallest redundant CRT modulus from V9.
    // 2. Keep four NTT-friendly moduli whose product still covers the
    //    worst-case integer convolution bound for input 4.
    // 3. Use Montgomery multiplication inside each CRT NTT butterfly.
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
            crtWorkA[modulusIndex][index] =
                static_cast<int>(firstPolynomial[index] % modulus);

            crtWorkB[modulusIndex][index] =
                static_cast<int>(secondPolynomial[index] % modulus);
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
            crtResidues[modulusIndex][index] =
                crtWorkA[modulusIndex][index];
        }
    }

    long long inverseTable[CRT_MOD_COUNT][CRT_MOD_COUNT];

    for (int i = 0; i < CRT_MOD_COUNT; ++i) {
        for (int j = 0; j < CRT_MOD_COUNT; ++j) {
            inverseTable[i][j] = 0;
        }
    }

    for (int i = 0; i < CRT_MOD_COUNT; ++i) {
        for (int j = i + 1; j < CRT_MOD_COUNT; ++j) {
            inverseTable[i][j] = quickPowerModCrtInt(
                CRT_MODS[i] % CRT_MODS[j],
                CRT_MODS[j] - 2,
                CRT_MODS[j]
            );
        }
    }

    int resultLength = 2 * polynomialLength - 1;

    #pragma omp parallel for num_threads(threadCount) schedule(static)
    for (int index = 0; index < resultLength; ++index) {
        int residuesForIndex[CRT_MOD_COUNT];

        for (int modulusIndex = 0;
             modulusIndex < CRT_MOD_COUNT;
             ++modulusIndex) {

            residuesForIndex[modulusIndex] =
                crtResidues[modulusIndex][index];
        }

        resultPolynomial[index] = garnerCombineToTargetModulus(
            residuesForIndex,
            inverseTable,
            targetModulus
        );
    }
}


int a[300000], b[300000], ab[300000];
long long aLarge[300000];
long long bLarge[300000];
long long abLarge[300000];
int main(int argc, char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    const int threadCount = THREAD_COUNT;
    std::cout << "32-bit Pthread whole-pipeline + large-modulus CRT thread count = "<< threadCount << std::endl;
    int test_begin = 0;
    int test_end = 4;
    for(int i = test_begin; i <= test_end; ++i){
        if (i <= 3) {
            long double ans = 0;
            int n_, p_;
            fRead(a, b, &n_, &p_, i);
            memset(ab, 0, sizeof(ab));
            auto Start = std::chrono::high_resolution_clock::now();
            // TODO : 将 poly_multiply 函数替换成你写的 ntt
            multiplyByNttPthreadWhole(a, b, ab, n_, p_, threadCount);
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
            ans += elapsed.count();
            fCheck(ab, n_, i);
            std::cout<<"average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (ms) "<<std::endl;
            // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
            // 禁止使用 cout 一次性输出大量文件内容
            //fWrite(ab, n_, i);
        }else {
            long double ans = 0;
            int n_;
            long long p_;
            fReadLarge(aLarge, bLarge, &n_, &p_, i);
            memset(abLarge, 0, sizeof(abLarge));
            auto Start = std::chrono::high_resolution_clock::now();
            multiplyByNttLargeCRT(aLarge, bLarge, abLarge, n_, p_, threadCount);
            auto End = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double, std::ratio<1,1000>> elapsed =End - Start;
            ans += elapsed.count();
            fCheckLarge(abLarge, n_, i);
            std::cout << "large modulus CRT latency for n = "<< n_<< " p = "<< p_<< " : "<< ans<< " (ms) "<< std::endl;
        }
    }
    return 0;
}


