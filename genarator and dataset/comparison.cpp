#include <bits/stdc++.h>
using namespace std;

// ==========================================================
//                 PRINT RESULT TABLE (simple)
// ==========================================================
void print_table(const string& name, double ms, const vector<vector<int>>& /*mat*/) {
    // keep the same behaviour as your template: only report time
    cout << name << " completed in " << fixed << setprecision(3) << ms << " ms\n";
}

// ==========================================================
//          RANDOM DNA SEQUENCE GENERATOR
// ==========================================================
vector<string> generate_random_dna(int n, int m) {
    static const char dna[4] = {'A', 'C', 'G', 'T'};
    vector<string> result;
    result.reserve(n);

    for (int i = 0; i < n; i++) {
        string s;
        s.reserve(m);
        for (int j = 0; j < m; j++) {
            s.push_back(dna[rand() % 4]);
        }
        result.push_back(s);
    }
    return result;
}

// ==========================================================
// STEP 1: Naive Brute-Force APSPM
// ==========================================================
vector<vector<int>> compute_naive_APSPM(const vector<string>& arr) {
    int n = arr.size();
    vector<vector<int>> result(n, vector<int>(n, 0));

    auto suffix_prefix_match = [&](const string& a, const string& b) {
        int maxLen = min<int>(a.size(), b.size());
        for (int len = maxLen; len > 0; len--) {
            if (a.substr(a.size() - len) == b.substr(0, len)) {
                return len;
            }
        }
        return 0;
    };

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
                result[i][j] = suffix_prefix_match(arr[i], arr[j]);

    return result;
}

void apspm_bruteforce(const vector<string>& strings) {
    auto start = chrono::high_resolution_clock::now();
    auto r = compute_naive_APSPM(strings);
    auto end = chrono::high_resolution_clock::now();
    double ms = chrono::duration<double, milli>(end - start).count();
    print_table("Step 1: Naive Brute Force", ms, r);
}

// ==========================================================
// STEP 2: Optimized Brute-Force (DOUBLE HASHING + Binary Search)
// ==========================================================
struct HashedString {
    // double hashing
    static const long long MOD1 = 1000000007LL;
    static const long long MOD2 = 1000000009LL;
    static const long long BASE1 = 9113823LL;   // relatively small base
    static const long long BASE2 = 97266353LL;  // different base

    vector<long long> pref1, pref2, power1, power2;
    int n;

    HashedString() : n(0) {}
    HashedString(const string& s) {
        n = (int)s.size();
        pref1.assign(n + 1, 0);
        pref2.assign(n + 1, 0);
        power1.assign(n + 1, 0);
        power2.assign(n + 1, 0);

        power1[0] = 1;
        power2[0] = 1;
        for (int i = 0; i < n; i++) {
            int v = (unsigned char)s[i]; // use the character value
            pref1[i + 1] = ( (pref1[i] * BASE1) % MOD1 + v ) % MOD1;
            pref2[i + 1] = ( (pref2[i] * BASE2) % MOD2 + v ) % MOD2;
            power1[i + 1] = (power1[i] * BASE1) % MOD1;
            power2[i + 1] = (power2[i] * BASE2) % MOD2;
        }
    }

    // return pair hash for s[l..r-1]
    pair<long long,long long> getHash(int l, int r) const {
        if (l >= r) return {0,0};
        long long x1 = (pref1[r] - (pref1[l] * power1[r - l]) % MOD1) % MOD1;
        if (x1 < 0) x1 += MOD1;
        long long x2 = (pref2[r] - (pref2[l] * power2[r - l]) % MOD2) % MOD2;
        if (x2 < 0) x2 += MOD2;
        return {x1, x2};
    }
};

vector<vector<int>> compute_hashed_APSPM(const vector<string>& arr) {
    int n = arr.size();
    vector<vector<int>> result(n, vector<int>(n, 0));
    vector<HashedString> H; H.reserve(n);
    for (const auto &s : arr) H.emplace_back(s);

    auto match_len = [&](int i, int j) {
        const string &A = arr[i], &B = arr[j];
        int nA = A.size(), nB = B.size();
        int low = 0, high = min(nA, nB), best = 0;

        while (low <= high) {
                    int mid = (low + high) / 2;
                    if (mid > nA || mid > nB) {  // out of bounds
            high = mid - 1;
            continue;
        }
        auto h1 = H[i].getHash(nA - mid, nA);
        auto h2 = H[j].getHash(0, mid);

            if (h1 == h2) {
                best = mid;
                low = mid + 1;
            } else high = mid - 1;
        }
        return best;
    };

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (i != j)
                result[i][j] = match_len(i, j);

    return result;
}

void apspm_bruteforce_hash(const vector<string>& strings) {
    auto start = chrono::high_resolution_clock::now();
    auto r = compute_hashed_APSPM(strings);
    auto end = chrono::high_resolution_clock::now();
    double ms = chrono::duration<double, milli>(end - start).count();
    print_table("Step 2: Brute + Double-Hash", ms, r);
}

// ==========================================================
// STEP 3 — STATIC SUFFIX ARRAY + LCP + RMQ  (FULL CODE)
// ==========================================================
#include <climits>

void apspm_suffix_array(const vector<string>& strings) {

    auto start = chrono::high_resolution_clock::now();

    int nStr = strings.size();

    unordered_map<char,int> cmap = {{'A',1000},{'C',1001},{'G',1002},{'T',1003}};

    vector<int> concat;
    vector<int> start_pos(nStr), length(nStr);

    for(int i=0; i<nStr; i++){
        start_pos[i] = concat.size();
        length[i] = strings[i].size();
        for(char c : strings[i]) concat.push_back(cmap[c]);
        concat.push_back(1000 + 4 + i + 1); // ensure separators > alphabet
    }

    int N = concat.size();

    // ----- build SA (O(N log^2 N) simple stable approach) -----
    vector<int> sa(N), rank_(N), tmp(N);
    for(int i=0;i<N;i++){ sa[i]=i; rank_[i]=concat[i]; }

    for(int k=1;;k<<=1) {
        auto cmp = [&](int a, int b){
            if(rank_[a] != rank_[b]) return rank_[a] < rank_[b];
            int ra = (a+k < N)? rank_[a+k] : -1;
            int rb = (b+k < N)? rank_[b+k] : -1;
            return ra < rb;
        };
        sort(sa.begin(), sa.end(), cmp);

        tmp[sa[0]] = 0;
        for(int i=1;i<N;i++){
            tmp[sa[i]] = tmp[sa[i-1]] + ( cmp(sa[i-1], sa[i]) ? 1 : 0 );
        }
        for(int i=0;i<N;i++) rank_[i] = tmp[i];
        if(rank_[sa[N-1]] == N-1) break;
    }

    // ----- build LCP -----
    vector<int> rank_sa(N);
    for(int i=0;i<N;i++) rank_sa[sa[i]] = i;

    vector<int> lcp(max(0,N-1));
    int h = 0;
    for(int i=0;i<N;i++){
        int r = rank_sa[i];
        if(r == 0) continue;
        int j = sa[r-1];
        while(i+h<N && j+h<N && concat[i+h] == concat[j+h]) h++;
        lcp[r-1] = h;
        if(h>0) h--;
    }

    // ----- RMQ -----
    struct RMQ {
        int n, LOG;
        vector<vector<int>> st;
        vector<int> lg;

        void build(const vector<int>& arr){
            n = arr.size();
            if(n == 0) return;
            LOG = 32 - __builtin_clz(n);
            st.assign(LOG, vector<int>(n));
            st[0] = arr;

            for(int k=1;k<LOG;k++)
                for(int i=0;i+(1<<k)<=n;i++)
                    st[k][i] = min(st[k-1][i], st[k-1][i+(1<<(k-1))]);

            lg.assign(n+1,0);
            for(int i=2;i<=n;i++) lg[i] = lg[i/2]+1;
        }

        // returns min over [l, r-1] (assumes l < r)
        int query(int l,int r){
            if(l>r) swap(l,r);
            if(l==r) return INT_MAX;
            int len = r - l;
            int k = lg[len];
            return min(st[k][l], st[k][r - (1<<k)]);
        }
    } rmq;

    if(!lcp.empty()) rmq.build(lcp);

    auto lcp_query = [&](int x, int y){
        if(x==y) return INT_MAX;
        int rx = rank_sa[x], ry = rank_sa[y];
        if(rx > ry) swap(rx, ry);
        // LCP of suffixes x and y is min over lcp[rx .. ry-1] which matches rmq.query(rx, ry)
        return rmq.query(rx, ry);
    };

    vector<vector<int>> result(nStr, vector<int>(nStr,0));

    for(int i=0;i<nStr;i++){
        for(int j=0;j<nStr;j++){
            if(i==j) continue;
            int best = 0;
            int posJ = start_pos[j];
            int maxPossible = min(length[i], length[j]);

            int begin = start_pos[i];
            int end   = start_pos[i] + length[i] - 1;

            for(int p = begin; p <= end; p++){
                int val = lcp_query(p, posJ);
                if(val == INT_MAX) continue;
                int avail = min(start_pos[i] + length[i] - p, length[j]);
                best = max(best, min(val, avail));
                if(best == maxPossible) break;
            }

            result[i][j] = best;
        }
    }

    auto end = chrono::high_resolution_clock::now();
    double ms = chrono::duration<double, milli>(end - start).count();
    print_table("Step 3: Suffix Array", ms, result);
}

// ==========================================================
// STEP 4 — SEMI-DYNAMIC APSPM (Incremental insertion + rolling double-hash)
// ==========================================================
void apspm_semi_dynamic(const vector<string>& strings) {
    auto start = chrono::high_resolution_clock::now();

    int n = strings.size();
    vector<vector<int>> result(n, vector<int>(n, 0));

    vector<HashedString> hashed; hashed.reserve(n);

    auto match_suffix_prefix = [&](const HashedString &HA, const string &A,
                                   const HashedString &HB, const string &B) -> int {
        int nA = A.size(), nB = B.size();
        int low = 0, high = min(nA, nB), best = 0;
        while (low <= high) {
            int mid = (low + high) >> 1;
            auto h1 = (mid==0 ? make_pair(0LL,0LL) : HA.getHash(nA - mid, nA));
            auto h2 = (mid==0 ? make_pair(0LL,0LL) : HB.getHash(0, mid));
            if (h1 == h2) {
                best = mid;
                low = mid + 1;
            } else high = mid - 1;
        }
        return best;
    };

    for (int k = 0; k < n; ++k) {
        HashedString Hk(strings[k]);
        for (int i = 0; i < (int)hashed.size(); ++i) {
            result[i][k] = match_suffix_prefix(hashed[i], strings[i], Hk, strings[k]);
            result[k][i] = match_suffix_prefix(Hk, strings[k], hashed[i], strings[i]);
        }
        result[k][k] = 0;
        hashed.push_back(std::move(Hk));
    }

    auto end = chrono::high_resolution_clock::now();
    double ms = chrono::duration<double, milli>(end - start).count();
    print_table("Step 4: Semi-Dynamic APSPM (incremental)", ms, result);
}

// ==========================================================
// small verifier (compare matrices) - useful for debugging
// ==========================================================
void verify(const vector<vector<int>>& a, const vector<vector<int>>& b, const string& a_name, const string& b_name) {
    if(a.size() != b.size()) {
        cerr << "VERIFY: size mismatch between " << a_name << " and " << b_name << "\n";
        return;
    }
    int n = a.size();
    int mismatches = 0;
    for(int i=0;i<n;i++) for(int j=0;j<n;j++){
        if(a[i][j] != b[i][j]) {
            // only show a few mismatches
            if(++mismatches <= 10) {
                cerr << "Mismatch [" << i << "][" << j << "]: " << a_name << "=" << a[i][j] << "  " << b_name << "=" << b[i][j] << "\n";
            }
        }
    }
    if(mismatches) cerr << "Verification between " << a_name << " and " << b_name << " failed: " << mismatches << " mismatches (showing up to 10)\n";
    else cerr << "Verification " << a_name << " == " << b_name << " OK\n";
}

// ==========================================================
// RUN ALL ALGORITHMS
// ==========================================================
vector<long long> run_all(int n, int m) {
    vector<string> strings = generate_random_dna(n, m);
    vector<long long> times(4);

    // Run brute and keep result (for verification)
    vector<vector<int>> r0, r1, r2, r3;
    {
        auto t = chrono::high_resolution_clock::now();
        r0 = compute_naive_APSPM(strings);
        times[0] = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
        print_table("Step 1: Naive Brute Force", times[0], r0);
    }
    {
        auto t = chrono::high_resolution_clock::now();
        r1 = compute_hashed_APSPM(strings);
        times[1] = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
        print_table("Step 2: Brute + Double-Hash", times[1], r1);
    }
    {
        auto t = chrono::high_resolution_clock::now();
        apspm_suffix_array(strings); // prints internally (we won't capture its matrix here)
        times[2] = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
    }
    {
        auto t = chrono::high_resolution_clock::now();
        // run semi-dynamic and capture its result by calling the inner function variant:
        // I'll call compute_hashed_APSPM incremental here to capture; for clarity we'll call apspm_semi_dynamic (it prints time).
        // But to verify, we need its matrix. So reproduce similar logic to get the matrix returned:
        int nn = strings.size();
        vector<vector<int>> result(nn, vector<int>(nn, 0));
        vector<HashedString> hashed; hashed.reserve(nn);
        auto match_suffix_prefix = [&](const HashedString &HA, const string &A,
                                       const HashedString &HB, const string &B) -> int {
            int nA = A.size(), nB = B.size();
            int low = 0, high = min(nA, nB), best = 0;
            while (low <= high) {
                int mid = (low + high) >> 1;
                auto h1 = (mid==0 ? make_pair(0LL,0LL) : HA.getHash(nA - mid, nA));
                auto h2 = (mid==0 ? make_pair(0LL,0LL) : HB.getHash(0, mid));
                if (h1 == h2) {
                    best = mid;
                    low = mid + 1;
                } else high = mid - 1;
            }
            return best;
        };
        for (int k = 0; k < nn; ++k) {
            HashedString Hk(strings[k]);
            for (int i = 0; i < (int)hashed.size(); ++i) {
                result[i][k] = match_suffix_prefix(hashed[i], strings[i], Hk, strings[k]);
                result[k][i] = match_suffix_prefix(Hk, strings[k], hashed[i], strings[i]);
            }
            hashed.push_back(std::move(Hk));
        }
        r3 = result;
        apspm_semi_dynamic(strings); // prints time inside
        times[3] = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - t).count();
    }

    // Quick verification for smaller sizes (to avoid huge output)
    if (n <= 200) {
        verify(r0, r1, "Brute", "Hash");
        verify(r0, r3, "Brute", "SemiDynamic");
    } else {
        cerr << "Verification skipped for n > 200\n";
    }

    return times;
}

// ==========================================================
// MAIN: PRINT TABLE
// ==========================================================
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    srand((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());

    vector<pair<int,int>> test_cases = {
        {50,20}, {50,40}, {50,50}, {50,200}, {50,500}, {50,1000},
        {100,20},
    };

    cout << "Strings | Len | Brute(ms) | Hash(ms) | SuffixArray(ms) | SemiDynamic(ms)\n";
    cout << "---------------------------------------------------------------------\n";

    for (auto &tc : test_cases) {
        auto r = run_all(tc.first, tc.second);

        cout << setw(7) << tc.first << " | "
             << setw(4) << tc.second << " | "
             << setw(9) << r[0] << " | "
             << setw(8) << r[1] << " | "
             << setw(15) << r[2] << " | "
             << setw(14) << r[3] << "\n";
    }

    return 0;
}
