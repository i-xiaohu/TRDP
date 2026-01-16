//
// Created by ixiaohu on 2026/1/16.
//

#include <string>
#include <vector>
#include <cassert>
using namespace std;

void kmp_style(const string &s) {
	// KMP
	const int n = s.length();
	vector<int> lps(n); // Longest prefix-suffix match length
	lps[0] = 0;
	for (int i = 1, j = 0; i < n; ) {
		// s[0:j) = s[i-j: i)
		if (s[i] == s[j]) { // Extension of a match
			lps[i++] = ++j;
		} else {
			if (j == 0) { // Mismatch
				lps[i++] = 0;
			} else {
				j = lps[j-1]; // Backward jump
			}
		}
	}
	// Check correctness
	for (int i = 0; i < n; i++) {
		int l = lps[i];
		printf("i=%d, len=%d, s[%d:%d) = s[%d:%d]\n", i, l, 0, l, i-l+1, i);
		// s[0:l) = s[i-l+1:i]
		for (int j = 0; j < l; j++) {
			assert(s[j] == s[i-l+1+j]);
		}
		// I didn't test if a longer match exists
	}

	vector<int> mid(n, 0); // mid(i) = max{l | s[0:l) = s[i,i+l)}
	// Use _lps_ to infer _mid_
	for (int i = 0; i < n; i++) {
		int l = lps[i];
		if (l > 0) {
			mid[i-l+1] = max(mid[i-l+1], l);
		}
	}
	for (int i = 0; i < n; i++) {
		printf("mid[%d]=%d\n", i, mid[i]);
	}

}

int main(int argc, char *argv[]) {
	string s = "ababac";
	kmp_style(s);
	return 0;
}