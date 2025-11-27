// suffix array solution for incremental APSP (ASAP)
// Algorithm logic fully preserved — only input/output handling modified

#include <bits/stdc++.h>
using namespace std;

// ---------- SA-IS SUFFIX ARRAY IMPLEMENTATION (unchanged) ----------
void induced_sort(const vector<int> &vec, int val_range, vector<int> &SA,
                  const vector<bool> &sl, const vector<int> &lms_idx) {
    vector<int> l(val_range, 0), r(val_range, 0);
    for (int c : vec) {
        if (c + 1 < val_range) ++l[c + 1];
        ++r[c];
    }
    partial_sum(l.begin(), l.end(), l.begin());
    partial_sum(r.begin(), r.end(), r.begin());
    fill(SA.begin(), SA.end(), -1);

    for (int i = lms_idx.size() - 1; i >= 0; --i)
        SA[--r[vec[lms_idx[i]]]] = lms_idx[i];

    for (int i : SA)
        if (i >= 1 && sl[i - 1])
            SA[l[vec[i - 1]]++] = i - 1;

    fill(r.begin(), r.end(), 0);
    for (int c : vec) ++r[c];
    partial_sum(r.begin(), r.end(), r.begin());

    for (int k = SA.size() - 1, i = SA[k]; k >= 1; --k, i = SA[k])
        if (i >= 1 && !sl[i - 1])
            SA[--r[vec[i - 1]]] = i - 1;
}

vector<int> SA_IS(const vector<int> &vec, int val_range) {
    int n = vec.size();
    vector<int> SA(n), lms_idx;
    vector<bool> sl(n);

    sl[n - 1] = false;
    for (int i = n - 2; i >= 0; --i) {
        sl[i] = (vec[i] > vec[i + 1] || (vec[i] == vec[i + 1] && sl[i + 1]));
        if (sl[i] && !sl[i + 1]) lms_idx.push_back(i + 1);
    }
    reverse(lms_idx.begin(), lms_idx.end());
    induced_sort(vec, val_range, SA, sl, lms_idx);

    vector<int> new_lms_idx(lms_idx.size()), lms_vec(lms_idx.size());
    for (int i = 0, k = 0; i < n; ++i)
        if (!sl[SA[i]] && SA[i] >= 1 && sl[SA[i] - 1])
            new_lms_idx[k++] = SA[i];

    int cur = 0;
    SA[n - 1] = cur;

    for (size_t k = 1; k < new_lms_idx.size(); ++k) {
        int i = new_lms_idx[k - 1], j = new_lms_idx[k];
        if (vec[i] != vec[j]) {
            SA[j] = ++cur;
            continue;
        }
        bool flag = false;
        for (int a = i + 1, b = j + 1;; ++a, ++b) {
            if (vec[a] != vec[b]) {
                flag = true;
                break;
            }
            if ((!sl[a] && sl[a - 1]) || (!sl[b] && sl[b - 1])) {
                flag = !((!sl[a] && sl[a - 1]) && (!sl[b] && sl[b - 1]));
                break;
            }
        }
        SA[j] = (flag ? ++cur : cur);
    }

    for (size_t i = 0; i < lms_idx.size(); ++i)
        lms_vec[i] = SA[lms_idx[i]];

    if (cur + 1 < (int)lms_idx.size()) {
        auto lms_SA = SA_IS(lms_vec, cur + 1);
        for (size_t i = 0; i < lms_idx.size(); ++i)
            new_lms_idx[i] = lms_idx[lms_SA[i]];
    }

    induced_sort(vec, val_range, SA, sl, new_lms_idx);
    return SA;
}

vector<int> suffix_array(const string &s, int LIM = 128) {
    vector<int> vec(s.size() + 1);
    copy(s.begin(), s.end(), vec.begin());
    vec.back() = '$';
    auto ret = SA_IS(vec, LIM);
    ret.erase(ret.begin());
    return ret;
}

// ---------- SuffixArray Class (unchanged) ----------
struct SuffixArray {
    int n;
    string s;
    vector<int> sa, rank, lcp;
    static const int LG = 22;
    vector<vector<int>> t;
    vector<int> lg;

    SuffixArray() {}
    SuffixArray(string _s) {
        n = _s.size();
        s = _s;
        sa = suffix_array(s);
        rank.resize(n);
        for (int i = 0; i < n; i++) rank[sa[i]] = i;
        build_lcp();
        preprocess();
        build_rmq();
    }

    void build_lcp() {
        int k = 0;
        lcp.resize(n - 1);
        for (int i = 0; i < n; i++) {
            if (rank[i] == n - 1) { k = 0; continue; }
            int j = sa[rank[i] + 1];
            while (i + k < n && j + k < n && s[i + k] == s[j + k]) k++;
            lcp[rank[i]] = k;
            if (k) k--;
        }
    }

    void preprocess() {
        lg.resize(n, 0);
        for (int i = 2; i < n; i++)
            lg[i] = lg[i / 2] + 1;
    }

    void build_rmq() {
        int sz = n - 1;
        t.resize(sz, vector<int>(LG));
        for (int i = 0; i < sz; i++)
            t[i][0] = lcp[i];

        for (int k = 1; k < LG; k++)
            for (int i = 0; i + (1 << k) - 1 < sz; i++)
                t[i][k] = min(t[i][k - 1], t[i + (1 << (k - 1))][k - 1]);
    }

    int query(int l, int r) {
        int k = lg[r - l + 1];
        return min(t[l][k], t[r - (1 << k) + 1][k]);
    }
};

// ---------- GLOBALS USED BY YOUR ORIGINAL LOGIC ----------
vector<string> string_list;
vector<vector<int>> ans_grid;
vector<int> str_pos, suf_len;
string S;
int total_len = 0;

// ---------- Your find_solution() Logic (unchanged) ----------
void find_solution() {
    int n = string_list.size();
    ans_grid.assign(n, vector<int>(n, 0));

    SuffixArray SA(S);

    auto fill_row = [&](int row, int ind) {
        for (int j = ind - 1; j >= 0; j--) {
            int rmq = SA.query(j, ind - 1);
            if (rmq >= suf_len[SA.sa[j]]) {
                int col = str_pos[SA.sa[j]];
                ans_grid[row][col] = max(ans_grid[row][col], suf_len[SA.sa[j]]);
            }
        }
        for (int j = ind; j < total_len - 1; j++) {
            int rmq = SA.query(ind, j);
            if (rmq < string_list[row].size()) break;
            int col = str_pos[SA.sa[j + 1]];
            ans_grid[row][col] = max(ans_grid[row][col], suf_len[SA.sa[j + 1]]);
        }
    };

    for (int i = 0; i < total_len; i++) {
        int a = SA.sa[i];
        int id = str_pos[a];
        if (a == 0 || S[a - 1] == '@') {
            fill_row(id, i);
        }
    }
}

// ---------- Modified solve() for Your Required I/O ----------
void solve() {
    ifstream in("dna_dataset.txt");
    if (!in.is_open()) {
        cerr << "Error opening dna_dataset.txt\n";
        exit(1);
    }

    cout << "Running Incremental Suffix Array APSP...\n\n";

    while (true) {
        int N, M;
        if (!(in >> N >> M)) break;  // EOF

        string_list.clear();
        S.clear();
        str_pos.clear();
        suf_len.clear();
        total_len = 0;

        string_list.reserve(N);

        for (int i = 0; i < N; i++) {
            string s;
            in >> s;
            string_list.push_back(s + "@");

            int id = i;
            int L = s.size() + 1;

            for (int j = 0; j < L; j++) {
                str_pos.push_back(id);
                suf_len.push_back(L - 1 - j);
            }

            S += s + "@";
            total_len += L;
        }

        string sep;
        in >> sep;   // read "---"

        auto start = chrono::steady_clock::now();

        // incremental recomputation
        for (int i = 0; i < N; i++)
            find_solution();

        auto end = chrono::steady_clock::now();
        long long ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();

        cout << "N = " << N << ", M = " << M
             << "  →  Time: " << ms << " ms\n";
    }
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    solve();
    return 0;
}
