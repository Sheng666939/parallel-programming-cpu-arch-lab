#include <iostream>
#include <vector>

using namespace std;

vector<double> matVecCacheFriendly(const vector<double>& matrixData,
                                   const vector<double>& inputVec,
                                   int n) {
    vector<double> resultVec(n, 0.0);

    for (int row = 0; row < n; ++row) {
        double currentX = inputVec[row];
        int rowBase = row * n;
        for (int col = 0; col < n; ++col) {
            resultVec[col] += matrixData[rowBase + col] * currentX;
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

    vector<double> resultVec = matVecCacheFriendly(matrixData, inputVec, n);

    cout << "Result:\n";
    for (int col = 0; col < n; ++col) {
        cout << resultVec[col] << (col + 1 == n ? '\n' : ' ');
    }

    cout << "Expected:\n";
    cout << "90 100 110 120\n";

    return 0;
}
