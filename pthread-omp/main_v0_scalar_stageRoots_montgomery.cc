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
int a[300000], b[300000], ab[300000];
int main(int argc, char *argv[])
{
    
    // 保证输入的所有模数的原根均为 3, 且模数都能表示为 a \times 4 ^ k + 1 的形式
    // 输入模数分别为 7340033 104857601 469762049 263882790666241
    // 第四个模数超过了整型表示范围, 如果实现此模数意义下的多项式乘法需要修改框架
    // 对第四个模数的输入数据不做必要要求, 如果要自行探索大模数 NTT, 请在完成前三个模数的基础代码及优化后实现大模数 NTT
    // 输入文件共五个, 第一个输入文件 n = 4, 其余四个文件分别对应四个模数, n = 131072
    // 在实现快速数论变化前, 后四个测试样例运行时间较久, 推荐调试正确性时只使用输入文件 1
    int test_begin = 0;
    int test_end = 3;
    for(int i = test_begin; i <= test_end; ++i){
        long double ans = 0;
        int n_, p_;
        fRead(a, b, &n_, &p_, i);
        memset(ab, 0, sizeof(ab));
        auto Start = std::chrono::high_resolution_clock::now();
        // TODO : 将 poly_multiply 函数替换成你写的 ntt
        multiplyByNtt(a, b, ab, n_, p_);
        auto End = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::ratio<1,1000>>elapsed = End - Start;
        ans += elapsed.count();
        fCheck(ab, n_, i);
        std::cout<<"average latency for n = "<<n_<<" p = "<<p_<<" : "<<ans<<" (ms) "<<std::endl;
        // 可以使用 fWrite 函数将 ab 的输出结果打印到 files 文件夹下
        // 禁止使用 cout 一次性输出大量文件内容
        //fWrite(ab, n_, i);
    }
    return 0;
}
