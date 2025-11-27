#include <bits/stdc++.h>
using namespace std;

// Global storage
vector<string> string_list;
vector<vector<int>> ans_grid;

// Your original brute force logic
void find_solution() {
    int n = string_list.size();
    ans_grid.assign(n, vector<int>(n, 0));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {

            int a = min(string_list[i].size(), string_list[j].size());
            int b = string_list[j].size();

            for (int k = a; k >= 1; k--) {
                if (string_list[i].substr(0, k) == string_list[j].substr(b - k, k)) {
                    ans_grid[i][j] = k;
                    break;
                }
            }
        }
    }
}

int main() {

    ifstream in("dna_dataset.txt");
    if (!in.is_open()) {
        cerr << "Error: Could not open dna_dataset.txt\n";
        return 1;
    }

    cout << "Running Incremental Brute Force APSP...\n\n";

    while (true) {
        int N, M;
        if (!(in >> N >> M)) break;  // EOF

        string_list.clear();
        string_list.reserve(N);

        // Read the N DNA strings
        for (int i = 0; i < N; i++) {
            string s;
            in >> s;
            string_list.push_back(s);
        }

        string sep;
        in >> sep; // read "---"

        // Time measurement
        auto start = chrono::steady_clock::now();

        // Your incremental re-execution:
        // For every addition, recompute everything
        for (int i = 0; i < N; i++) {
            find_solution();
        }

        auto end = chrono::steady_clock::now();
        long long ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();

        cout << "N = " << N << ", M = " << M
             << "  →  Time: " << ms << " ms\n";
    }

    return 0;
}


