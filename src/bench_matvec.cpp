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

void matVecNaiveColDot(const vector<double>& matrixData,
                       const vector<double>& inputVec,
                       vector<double>& resultVec,
                       int n) {
    fill(resultVec.begin(), resultVec.end(), 0.0);

    for (int col = 0; col < n; ++col) {
        for (int row = 0; row < n; ++row) {
            resultVec[col] += matrixData[row * n + col] * inputVec[row];
        }
    }
}

void matVecCacheFriendly(const vector<double>& matrixData,
                         const vector<double>& inputVec,
                         vector<double>& resultVec,
                         int n) {
    fill(resultVec.begin(), resultVec.end(), 0.0);

    for (int row = 0; row < n; ++row) {
        double currentX = inputVec[row];
        int rowBase = row * n;
        for (int col = 0; col < n; ++col) {
            resultVec[col] += matrixData[rowBase + col] * currentX;
        }
    }
}

double calcChecksum(const vector<double>& dataVec) {
    double totalSum = 0.0;
    for (double value : dataVec) {
        totalSum += value;
    }
    return totalSum;
}

int main() {
    vector<int> sizeList = {128, 256, 512, 1024, 1536, 2048};

    cout << "n,algo,repeats,total_ms,avg_ms,checksum\n";
    cout << fixed << setprecision(6);

    for (int n : sizeList) {
        vector<double> matrixData(static_cast<size_t>(n) * n);
        vector<double> inputVec(n);
        vector<double> resultVec(n, 0.0);

        for (int row = 0; row < n; ++row) {
            inputVec[row] = (row % 100) * 0.01 + 1.0;
            for (int col = 0; col < n; ++col) {
                matrixData[static_cast<size_t>(row) * n + col] = ((row + col) % 200) * 0.01;
            }
        }

        long long repeatCount = 0;
        double startMs = nowMs();

        do {
            matVecNaiveColDot(matrixData, inputVec, resultVec, n);
            ++repeatCount;
        } while (nowMs() - startMs < 200.0);

        double totalMs = nowMs() - startMs;
        cout << n << ",naive," << repeatCount << "," << totalMs << ","
             << totalMs / repeatCount << "," << calcChecksum(resultVec) << "\n";

        repeatCount = 0;
        startMs = nowMs();

        do {
            matVecCacheFriendly(matrixData, inputVec, resultVec, n);
            ++repeatCount;
        } while (nowMs() - startMs < 200.0);

        totalMs = nowMs() - startMs;
        cout << n << ",cache," << repeatCount << "," << totalMs << ","
             << totalMs / repeatCount << "," << calcChecksum(resultVec) << "\n";
    }

    return 0;
}
