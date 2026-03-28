#include <iostream>
#include <vector>

using namespace std;

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
    vector<double> dataVec = {1, 2, 3, 4, 5, 6, 7, 8};

    double resultValue = sumSuperscalarUnroll4(dataVec);

    cout << "Result:\n";
    cout << resultValue << "\n";

    cout << "Expected:\n";
    cout << "36\n";

    return 0;
}
