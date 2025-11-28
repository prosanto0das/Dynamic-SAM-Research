#include <bits/stdc++.h>
using namespace std;

string randomDNA(int len) {
    static const char bases[] = {'A', 'T', 'G', 'C'};
    string s;
    s.reserve(len);
    for (int i = 0; i < len; i++)
        s += bases[rand() % 4];
    return s;
}

int main() {
    srand(time(NULL));

    // You can modify these test cases as needed
    vector<pair<int,int>> testCases = {
        // {50, 20},      // N=50,  M=20
		// {50, 40},      // N=50,  M=40
		// {50, 50},      // N=50,  M=50
		// {50, 200},     // N=50,  M=200
		// {50, 500},     // N=50,  M=500
		// {50, 1000}   // N=50,  M=1000

		// {100, 20},     // N=100, M=20
		// {100, 40},     // N=100, M=40
		// {100, 50},     // N=100, M=50
		// {100, 200},    // N=100, M=200
		// {100, 500},    // N=100, M=500
		// {100, 1000}   // N=100, M=1000

        // {300, 20},     // N=100, M=20
        // {300, 40},     // N=100, M=40
        // {300, 50},     // N=100, M=50
        // {300, 200},    // N=100, M=200
        // {300, 500},    // N=100, M=500
        // {300, 1000}   // N=100, M=1000

        {1000, 20},     // N=100, M=20
        {1000, 40},     // N=100, M=40
        {1000, 50},     // N=100, M=50
        {1000, 200},    // N=100, M=200
        {1000, 500},    // N=100, M=500
        {1000, 1000}    // N=100, M=1000


    };

    ofstream out("dna_dataset.txt");

    if (!out.is_open()) {
        cerr << "Error opening output file!" << endl;
        return 1;
    }

    for (auto &tc : testCases) {
        int N = tc.first;
        int M = tc.second;

        out << N << " " << M << "\n";

        for (int i = 0; i < N; i++) {
            out << randomDNA(M) << "\n";
        }

        out << "---\n";   // separator for different test cases
    }

    out.close();
    cout << "DNA dataset generated successfully in dna_dataset.txt\n";

    return 0;
}
