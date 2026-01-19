//
// Created by ixiaohu on 2026/1/16.
//

#include <string>
#include <vector>
#include <cassert>
using namespace std;

// Longest prefix-substring match within S
// Let K = LPS[i], S[0:K) = S[i:i+K)
vector<int> calc_lps(const string &s) {
	const int n = s.length();
	vector<int> lps(n); // lps[i] = max{len | S[0:len)=S[i:i+len)}
	lps[0] = 0;
	int len = 0;
	while (1+len < n and s[0+len] == s[1+len]) len++;
	lps[1] = len;

	int k = 1;
	for (int i = 2; i < n; i++) {
		int u = k + lps[k];
		if (u > i) {
			int t = u - i;
			// S[i:u) equals to S[l-t:l), which serves as a bridge
			for (int j = 0; j < t; j++) {
				assert(s[i+j] == s[lps[k]-t+j]);
			}
			// Compare with prefix
			int x = lps[k] - t;
			if (lps[x] < t) {
				lps[i] = lps[x];
			} else {
				// The match can further extend
				len = t;
				while (i + len < n and s[0+len] == s[i+len]) len++;
				lps[i] = len;
			}
		} else {
			len = 0;
			while (i + len < n and s[0+len] == s[i+len]) len++;
			lps[i] = len;
		}
		if (i + lps[i] > k + lps[k]) {
			k = i;
		}
	}

	// Check correctness
	for (int i = 1; i < n; i++) {
//		printf("i=%d, lps=%d\n", i, lps[i]);
		len = lps[i];
		for (int j = 0; j < len; j++) {
			assert(s[0+j] == s[i+j]);
		}
		assert(i+len == n or s[len] != s[i+len]);
	}
	return lps;
}

// Longest substring-prefix match between T and Q
// Let LSP[i] = K, T[i:i+K) = Q[0:K)
void calc_lsp(const string &t, const string &q) {
	const int n = t.length();
	const int m = q.length();
	vector<int> lsp(n);
	for (int i = 0; i < n and i < m; i++) {
		if (q[i] != t[i]) {
			lsp[0] = i;
			break;
		}
	}
	vector<int> tmp = calc_lps(q);
	int k = 0;
	for (int i = 1; i < n; i++) {
		int u = k + lsp[k];
		if (u > i) {
			// Similarly, T[i:i+t) = Q[l-t:l) is the bridge
			int tt = u - i;
			for (int j = 0; j < tt; j++) {
				assert(t[i+j] == q[lsp[k]-tt+j]);
			}
			int x = lsp[k] - tt;
			if (tmp[x] < tt) {
				lsp[i] = tmp[x];
			} else {
				int len = tt;
				while (i + len < n and len < m and q[0+len] == t[i+len]) len++;
				lsp[i] = len;
			}
		} else {
			int len = 0;
			while (i + len < n and len < m and q[0+len] == t[i+len]) len++;
			lsp[i] = len;
		}
		if (i + lsp[i] > k + lsp[k]) {
			k = i;
		}
	}
	// Check correctness
	for (int i = 1; i < n; i++) {
//		printf("i=%d, lsp=%d\n", i, lsp[i]);
		int len = lsp[i];
		for (int j = 0; j < len; j++) {
			assert(q[0+j] == t[i+j]);
		}
		assert(i+len == n or len == m or q[len] != t[i+len]);
	}
}

int main(int argc, char *argv[]) {
	srand(time(nullptr));
	string s, q;
	for (int i = 0; i < 1000; i++) {
		s += "ACGT"[rand() % 4];
	}
	calc_lps(s);
	calc_lps(q);
	for (int i = 0; i < 200; i++) {
		q += "ACGT"[rand() % 4];
	}
	calc_lsp(s, q);
	return 0;
}