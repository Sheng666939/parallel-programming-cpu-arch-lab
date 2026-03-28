#include <iostream>
#include <vector>
#include <iomanip>
#include <sys/time.h>

using namespace std;

double nowMs() {
    timeval timeValue;
    gettimeofday(&timeValue, nullptr);
    return timeValue.tv_sec * 1000.0 + timeValue.tv_usec / 1000.0;
}

double sumNaive(const vector<double>& dataVec) {
    double totalSum = 0.0;
    for (size_t idx = 0; idx < dataVec.size(); ++idx) {
        totalSum += dataVec[idx];
    }
    return totalSum;
}

double sumSuperscalarUnroll4(const vector<double>& dataVec) {
    double sum0 = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;

    size_t idx = 0;
    size_t limit = dataVec.size() / 4 * 4;

    for (; idx < limit; idx += 4) {
        sum0 += dataVec[idx];
        sum1 += dataVec[idx + 1];
        sum2 += dataVec[idx + 2];
        sum3 += dataVec[idx + 3];
    }

    for (; idx < dataVec.size(); ++idx) {
        sum0 += dataVec[idx];
    }

    return sum0 + sum1 + sum2 + sum3;
}

int main() {
    vector<int> sizeList = {100000, 500000, 1000000, 2000000, 5000000, 10000000};

    cout << "n,algo,repeats,total_ms,avg_ms,result\n";
    cout << fixed << setprecision(6);

    for (int n : sizeList) {
        vector<double> dataVec(n);

        for (int idx = 0; idx < n; ++idx) {
            dataVec[idx] = (idx % 100) * 0.01 + 1.0;
        }

        long long repeatCount = 0;
        double resultValue = 0.0;
        double startMs = nowMs();

        do {
            resultValue = sumNaive(dataVec);
            ++repeatCount;
        } while (nowMs() - startMs < 300.0);

        double totalMs = nowMs() - startMs;
        cout << n << ",naive," << repeatCount << "," << totalMs << ","
             << totalMs / repeatCount << "," << resultValue << "\n";

        repeatCount = 0;
        startMs = nowMs();

        do {
            resultValue = sumSuperscalarUnroll4(dataVec);
            ++repeatCount;
        } while (nowMs() - startMs < 300.0);

        totalMs = nowMs() - startMs;
        cout << n << ",superscalar_unroll4," << repeatCount << "," << totalMs << ","
             << totalMs / repeatCount << "," << resultValue << "\n";
    }

    return 0;
}
