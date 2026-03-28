#include <iostream>
#include <vector>

using namespace std;

double sumNaive(const vector<double>& dataVec) {
    double totalSum = 0.0;
    for (size_t idx = 0; idx < dataVec.size(); ++idx) {
        totalSum += dataVec[idx];
    }
    return totalSum;
}

int main() {
    vector<double> dataVec = {1, 2, 3, 4, 5, 6, 7, 8};

    double resultValue = sumNaive(dataVec);

    cout << "Result:\n";
    cout << resultValue << "\n";

    cout << "Expected:\n";
    cout << "36\n";

    return 0;
}
