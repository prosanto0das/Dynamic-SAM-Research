// // ==========================================================
// // STEP 4 — SEMI-DYNAMIC APSPM (Incremental insertion + rolling hash)
// // ==========================================================
// void apspm_semi_dynamic(const vector<string>& strings) {
//     // This semi-dynamic approach inserts strings one-by-one.
//     // For each newly inserted string s_k, it computes:
//     //   - overlaps result[i][k] for all previous i (suffix of i -> prefix of k)
//     //   - overlaps result[k][i] for all previous i (suffix of k -> prefix of i)
//     //
//     // Complexity: for inserting k-th string, it does O(k * log m) hash checks
//     // (binary search per pair). Overall O(n^2 log m) in worst-case but avoids
//     // global rebuilds and is incremental / online-friendly.
//     //
//     // Note: if you want stronger collision guarantees, switch HashedString to
//     // a pair of moduli.

//     auto start = chrono::high_resolution_clock::now();

//     int n = strings.size();
//     vector<vector<int>> result(n, vector<int>(n, 0));

//     // Reuse HashedString (defined earlier in the template).
//     // Build an incremental vector of hashed strings as we insert:
//     vector<HashedString> hashed; hashed.reserve(n);

//     // Helper lambda: compute overlap length of suffix(A) with prefix(B)
//     auto match_suffix_prefix = [&](const HashedString &HA, const string &A,
//                                    const HashedString &HB, const string &B) -> int {
//         int nA = A.size(), nB = B.size();
//         int low = 0, high = min(nA, nB), best = 0;
//         while (low <= high) {
//             int mid = (low + high) >> 1;
//             // hash of suffix A[nA-mid .. nA)
//             long long h1 = (mid == 0 ? 0 : HA.getHash(nA - mid, nA));
//             // hash of prefix B[0 .. mid)
//             long long h2 = (mid == 0 ? 0 : HB.getHash(0, mid));
//             if (h1 == h2) {
//                 best = mid;
//                 low = mid + 1;
//             } else high = mid - 1;
//         }
//         return best;
//     };

//     // Insert strings incrementally
//     for (int k = 0; k < n; ++k) {
//         // prepare hashed for the new string
//         HashedString Hk(strings[k]);
//         // compute overlaps between new string k and all previous i
//         for (int i = 0; i < (int)hashed.size(); ++i) {
//             // hashed[i] corresponds to strings[i]
//             // 1) suffix(strings[i]) -> prefix(strings[k])
//             result[i][k] = match_suffix_prefix(hashed[i], strings[i], Hk, strings[k]);

//             // 2) suffix(strings[k]) -> prefix(strings[i])
//             result[k][i] = match_suffix_prefix(Hk, strings[k], hashed[i], strings[i]);
//         }

//         // Also, for completeness, set result[k][k] = 0 (no self-overlap)
//         result[k][k] = 0;

//         // push the new hashed string into vector for future insertions
//         hashed.push_back(std::move(Hk));
//     }

//     auto end = chrono::high_resolution_clock::now();
//     double ms = chrono::duration<double, milli>(end - start).count();

//     // optional: you may want to print or save the result here.
//     // For compatibility with the rest of the template, call print_table:
//     print_table("Step 4: Semi-Dynamic APSPM (incremental)", ms, result);

//     // If you want to debug the result matrix, uncomment below (careful for large n):
//     /*
//     for (int i = 0; i < n; ++i) {
//         for (int j = 0; j < n; ++j) {
//             cerr << result[i][j] << " ";
//         }
//         cerr << "\n";
//     }
//     */
//     // Visual reference image (if you want to use it externally):
//     // /mnt/data/4b235fc6-28cc-4933-bdf7-e0277d3bb5ba.png
// }


#include <bits/stdc++.h>
using namespace std;
using namespace std::chrono; // for chrono::steady_clock

// --- Dynamic Suffix Tree Node ---
struct Node {
    int start, *end;
    unordered_map<char, Node*> children;
    Node* suffixLink;
    vector<int> string_ids_here; // track which strings reach this node
    Node(int s, int *e) : start(s), end(e), suffixLink(nullptr) {}
};

// --- Dynamic Suffix Tree (Ukkonen’s Algorithm) ---
struct SuffixTree {
    string S;
    Node* root;
    Node* activeNode;
    int activeEdge, activeLength;
    int remainingSuffix;
    int leafEnd;
    int size;

    SuffixTree() {
        root = new Node(-1, new int(-1));
        root->suffixLink = root;
        activeNode = root;
        activeEdge = -1;
        activeLength = 0;
        remainingSuffix = 0;
        leafEnd = -1;
        size = 0;
    }

    int edgeLength(Node* n) { return *(n->end) - n->start + 1; }

    void extend(char c, int pos, int cur_str_id, vector<vector<int>>& ans_grid) {
        S += c;
        leafEnd = pos;
        remainingSuffix++;
        Node* lastNewNode = nullptr;

        while (remainingSuffix > 0) {
            if (activeLength == 0) activeEdge = pos;
            char edgeChar = S[activeEdge];

            if (!activeNode->children.count(edgeChar)) {
                Node* leaf = new Node(pos, &leafEnd);
                leaf->string_ids_here.push_back(cur_str_id);
                activeNode->children[edgeChar] = leaf;

                Node* temp = activeNode;
                int length = 1;
                while (temp && !temp->string_ids_here.empty()) {
                    for (int prev_id : temp->string_ids_here) {
                        if (prev_id != cur_str_id)
                            ans_grid[prev_id][cur_str_id] = max(ans_grid[prev_id][cur_str_id], length);
                    }
                    temp = temp->suffixLink;
                    length++;
                }

                if (lastNewNode) lastNewNode->suffixLink = activeNode;
                lastNewNode = nullptr;
            } else {
                Node* next = activeNode->children[edgeChar];
                int len = edgeLength(next);
                if (activeLength >= len) {
                    activeEdge += len;
                    activeLength -= len;
                    activeNode = next;
                    continue;
                }
                if (S[next->start + activeLength] == c) {
                    activeLength++;
                    if (lastNewNode) lastNewNode->suffixLink = activeNode;
                    break;
                }

                int* splitEnd = new int(next->start + activeLength - 1);
                Node* split = new Node(next->start, splitEnd);
                split->string_ids_here.insert(split->string_ids_here.end(),
                                             next->string_ids_here.begin(), next->string_ids_here.end());
                activeNode->children[edgeChar] = split;
                split->children[c] = new Node(pos, &leafEnd);
                split->children[c]->string_ids_here.push_back(cur_str_id);
                next->start += activeLength;
                split->children[S[next->start]] = next;

                if (lastNewNode) lastNewNode->suffixLink = split;
                lastNewNode = split;
            }

            remainingSuffix--;
            if (activeNode == root && activeLength > 0) {
                activeLength--;
                activeEdge = pos - remainingSuffix + 1;
            } else if (activeNode != root) {
                activeNode = activeNode->suffixLink;
            }
        }
    }
};

// --- APSPM Using Dynamic Suffix Tree ---
vector<string> string_list;
vector<vector<int>> ans_grid;

void find_solution_dynamic_tree() {
    int n = string_list.size();
    ans_grid.assign(n, vector<int>(n, 0));

    SuffixTree tree;
    int tot_len = 0;

    for (int s_idx = 0; s_idx < n; s_idx++) {
        string s = string_list[s_idx] + '@'; // separator
        for (int j = 0; j < (int)s.size(); j++) {
            tree.extend(s[j], tot_len++, s_idx, ans_grid);
        }
    }
}

int main() {
    ifstream in("dna_dataset.txt");
    if (!in.is_open()) {
        cerr << "Error: Could not open dna_dataset.txt\n";
        return 1;
    }

    cout << "Running Dynamic Suffix Tree APSPM...\n\n";

    while (true) {
        int N, M;
        if (!(in >> N >> M)) break;  // EOF

        string_list.clear();
        string_list.reserve(N);

        for (int i = 0; i < N; i++) {
            string s;
            in >> s;
            string_list.push_back(s);
        }

        string sep;
        in >> sep; // read "---"

        auto start = steady_clock::now();
        find_solution_dynamic_tree();
        auto end = steady_clock::now();
        auto ms = duration_cast<milliseconds>(end - start).count();

        cout << "N = " << N << ", M = " << M
             << "  →  Time: " << ms << " ms\n";
    }
}
