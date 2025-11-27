// treap_apspm.cpp
// Treap (implemented via std::set with comparator) + rolling double-hash
// Incremental APSP: for each newly appended string, insert its suffixes into an ordered set
// and scan neighbors to compute overlaps. Input/Output format matches your previous program:
// reads dna_dataset.txt with blocks: N M  [N strings] ---
// prints: N = X, M = Y  ->  Time: Z ms

#include <bits/stdc++.h>
using namespace std;

/* ------------------ Rolling double-hash that supports append ------------------ */
struct RollingHasher {
    static const long long MOD1 = 1000000007LL;
    static const long long MOD2 = 1000000009LL;
    static const long long BASE1 = 91138233LL % MOD1;
    static const long long BASE2 = 97266353LL % MOD2;

    vector<long long> pref1, pref2, pow1, pow2;
    int n;
    RollingHasher() : n(0) {
        pref1 = {0};
        pref2 = {0};
        pow1  = {1};
        pow2  = {1};
    }

    void append_char(char c) {
        int v = (unsigned char)c;
        long long p1 = (pref1.back() * BASE1 + v) % MOD1;
        long long p2 = (pref2.back() * BASE2 + v) % MOD2;
        pref1.push_back(p1);
        pref2.push_back(p2);
        pow1.push_back((pow1.back() * BASE1) % MOD1);
        pow2.push_back((pow2.back() * BASE2) % MOD2);
        n++;
    }

    // hash of substring [l, r) ; 0-based
    pair<long long,long long> getHash(int l, int r) const {
        if (l >= r) return {0,0};
        long long x1 = (pref1[r] - (pref1[l] * pow1[r - l]) % MOD1) % MOD1;
        if (x1 < 0) x1 += MOD1;
        long long x2 = (pref2[r] - (pref2[l] * pow2[r - l]) % MOD2) % MOD2;
        if (x2 < 0) x2 += MOD2;
        return {x1, x2};
    }

    // LCP length of suffixes starting at a and b (positions in current S)
    int lcp_between(const string &S, int a, int b) const {
        if (a == b) return (int)S.size() - a;
        int maxLen = min((int)S.size() - a, (int)S.size() - b);
        int lo = 0, hi = maxLen, best = 0;
        while (lo <= hi) {
            int mid = (lo + hi) >> 1;
            auto ha = getHash(a, a + mid);
            auto hb = getHash(b, b + mid);
            if (ha == hb) {
                best = mid;
                lo = mid + 1;
            } else hi = mid - 1;
        }
        return best;
    }
};

/* ------------------ Global state used by comparator and algorithm ------------------ */
static string G_S;
static RollingHasher G_hash;
static vector<int> G_str_pos; // which string each position belongs to
static vector<int> G_suf_len; // distance from pos to end-of-string (like saa)
static vector<string> G_string_list;

/* comparator for suffixes: lexicographic order of suffix starting at positions a and b */
struct SuffixComp {
    bool operator()(int a, int b) const {
        if (a == b) return false;
        int l = G_hash.lcp_between(G_S, a, b);
        int na = (int)G_S.size() - a;
        int nb = (int)G_S.size() - b;
        if (l >= na || l >= nb) {
            // one suffix is prefix of the other: shorter suffix is smaller
            return na < nb;
        }
        // compare next character
        return G_S[a + l] < G_S[b + l];
    }
};

/* ------------------ APSP storage ------------------ */
static vector<vector<int>> ans_grid; // ans_grid[i][j] = overlap len from i -> j

/* Insert many positions (suffix starts) into set and update structures */
void insert_positions_into_set(set<int, SuffixComp> &st, const vector<int> &pos_list) {
    for (int p : pos_list) {
        st.insert(p);
    }
}

/* For one start position (the first char of the newly added string), scan neighbors and update ans_grid */
void fill_row_using_set(set<int, SuffixComp> &st, int start_pos, int row_id) {
    auto it = st.find(start_pos);
    if (it == st.end()) return;

    // scan left (predecessors)
    auto itl = it;
    while (itl != st.begin()) {
        --itl;
        int q = *itl;
        int lcp = G_hash.lcp_between(G_S, start_pos, q);
        int avail = G_suf_len[q];
        if (lcp >= avail) {
            int col = G_str_pos[q];
            ans_grid[row_id][col] = max(ans_grid[row_id][col], avail);
        } else {
            break; // left scan stops when condition fails
        }
    }

    // scan right (successors)
    auto itr = it;
    ++itr;
    int need_full = (int)G_string_list[row_id].size();
    while (itr != st.end()) {
        int q = *itr;
        int lcp = G_hash.lcp_between(G_S, start_pos, q);
        if (lcp >= need_full) {
            int col = G_str_pos[q];
            ans_grid[row_id][col] = max(ans_grid[row_id][col], G_suf_len[q]);
            ++itr;
        } else {
            break;
        }
    }
}

/* run incremental treap-based APSP for a single test-case (N strings already read in G_string_list) */
long long run_incremental_with_set(int N) {
    // prepare global arrays
    ans_grid.assign(N, vector<int>(N, 0));
    G_str_pos.clear();
    G_suf_len.clear();
    G_S.clear();
    // Reset hasher
    G_hash = RollingHasher();

    set<int, SuffixComp> st; // order suffixes lexicographically

    int cur_pos = 0;
    // incremental insertion for each string
    for (int id = 0; id < N; ++id) {
        const string &s = G_string_list[id];
        int L = (int)s.size();
        // collect the positions within S that will be inserted
        vector<int> pos_list;
        pos_list.reserve(L + 1);

        // append s and '@' to G_S and update hasher and metadata
        for (int i = 0; i < L; ++i) {
            G_S.push_back(s[i]);
            G_hash.append_char(s[i]);
            // for the newly appended character at position cur_pos:
            G_str_pos.push_back(id);
            G_suf_len.push_back(L - i); // distance to string end (excluding separator)
            pos_list.push_back(cur_pos);
            cur_pos++;
        }
        // separator '@'
        G_S.push_back('@');
        G_hash.append_char('@');
        G_str_pos.push_back(id);
        G_suf_len.push_back(0); // separator position has suf_len 0 (not used)
        pos_list.push_back(cur_pos);
        cur_pos++;

        // insert all suffix-start positions that belong to this appended block
        insert_positions_into_set(st, pos_list);

        // call fill_row on the start position of this string (first char pos)
        int start_pos = pos_list[0];
        fill_row_using_set(st, start_pos, id);
    }

    // This implementation follows the "incremental updating after each inserted string"
    // We return no matrix but we measured running time externally.
    return 0;
}

/* ------------------ Main: read dna_dataset.txt and run ------------------ */
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream in("dna_dataset.txt");
    if (!in.is_open()) {
        cerr << "Error: Could not open dna_dataset.txt\n";
        return 1;
    }

    cout << "Running Treap-based (ordered-set) Incremental APSP...\n\n";

    while (true) {
        int N, M;
        if (!(in >> N >> M)) break; // EOF
        G_string_list.clear();
        G_string_list.reserve(N);
        for (int i = 0; i < N; ++i) {
            string s;
            in >> s; // read string of length M
            G_string_list.push_back(s);
        }
        string sep;
        in >> sep; // read '---'

        // measure time
        auto start = chrono::steady_clock::now();

        // run the treap-based insertion & update algorithm
        // For timing we will implement the same outer behavior: call incremental insertion,
        // but function run_incremental_with_set executes the incremental logic.
        run_incremental_with_set(N);

        auto end = chrono::steady_clock::now();
        long long ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();

        cout << "N = " << N << ", M = " << M
             << "  ->  Time: " << ms << " ms\n";
    }

    return 0;
}
