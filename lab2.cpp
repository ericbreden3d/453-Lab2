#include <iostream>
#include <stdlib.h>
#include <time.h>
using namespace std;

void make_matrix(int n) {
    srand(time(NULL));
    int m[n] = {};
    for (int i = 0; i < n; i++) {
        cout << m[i] << endl;
    }
}

int main() {
    make_matrix();
}