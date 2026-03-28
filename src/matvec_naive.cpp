#include <iostream>
#include <vector>

using namespace std;

vector<double> matVecNaiveColDot(const vector<double>& matrixData,
                                 const vector<double>& inputVec,
                                 int n) {
    vector<double> resultVec(n, 0.0);

    for (int col = 0; col < n; ++col) {
        for (int row = 0; row < n; ++row) {
            resultVec[col] += matrixData[row * n + col] * inputVec[row];
        }
    }

    return resultVec;
}

int main() {
    const int n = 4;

    vector<double> matrixData = {
        1,  2,  3,  4,
        5,  6,  7,  8,
        9, 10, 11, 12,
        13, 14, 15, 16
    };

    vector<double> inputVec = {1, 2, 3, 4};

    vector<double> resultVec = matVecNaiveColDot(matrixData, inputVec, n);

    cout << "Result:\n";
    for (int col = 0; col < n; ++col) {
        cout << resultVec[col] << (col + 1 == n ? '\n' : ' ');
    }

    cout << "Expected:\n";
    cout << "90 100 110 120\n";

    return 0;
}
