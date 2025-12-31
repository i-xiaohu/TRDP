//
// Created by ixiaohu on 2025/12/30.
//

#include <vector>
#include <cstdint>
#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <iostream>
using namespace std;

struct TestEntity {
	string motif;
	int period;
	int mutation;
	int flank_l;
	int flank_r;
	string seq;
};

TestEntity input_csv_test_seq(int n, const char *fn)
{
	ifstream in(fn);
	assert(in.is_open());
	string line;
	getline(in, line); // Header
	for (int i = 0; i < n-1; i++) {
		getline(in, line);
	}
	getline(in, line);
	for (char &a : line) {
		if (a == ',') {
			a = ' ';
		}
	}
	stringstream ss(line);
	int id; ss >> id;
	string motif; ss >> motif;
	int period; ss >> period;
	int mutation; ss >> mutation;
	int flank_l; ss >> flank_l;
	int flank_r; ss >> flank_r;
	string seq; ss >> seq;
	in.close();
	return TestEntity{motif, period, mutation, flank_l, flank_r, seq};
}

void solve(int n, const char *seq)
{
	int match_score = 1;
	int mismatch_penalty = -4;
	int open_gap_penalty = -6;
	int extend_gap_penalty = -1;
	vector<vector<int>> E;
	vector<vector<int>> F;
	vector<vector<int>> H;
	E.resize(n + 1);
	F.resize(n + 1);
	H.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		E[i].resize(i + 1, INT32_MIN);
		F[i].resize(i + 1, INT32_MIN);
		H[i].resize(i + 1, INT32_MIN);
	}
	// Global alignment in left-down triangle
	H[0][0] = 0;
	for (int i = 1; i <= n; i++) {
		H[i][0] = E[i][0] = open_gap_penalty + extend_gap_penalty * i;
		for (int j = 1; j < i; j++) {
			F[i][j] = max(H[i][j-1] + open_gap_penalty, F[i][j-1]) + extend_gap_penalty;
			E[i][j] = max(H[i-1][j] + open_gap_penalty, E[i-1][j]) + extend_gap_penalty;
			int M = H[i-1][j-1] + (seq[i-1] == seq[j-1] ?match_score :mismatch_penalty);
			H[i][j] = max(E[i][j], F[i][j]);
			H[i][j] = max(H[i][j], M);
		}
		H[i][i] = max(H[i][i-1] + open_gap_penalty, F[i][i-1]) + extend_gap_penalty;
		H[i][i] = max(H[i-1][i-1] + match_score, H[i][i]);
	}
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j <= i; j++) {
			fprintf(stdout, "%d ", H[i][j]);
		}
		fprintf(stdout, "\n");
	}
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		fprintf(stderr, "TRDP <CSV> <ID>\n");
		return 1;
	}
	const char *fn = argv[1];
	int id = atoi(argv[2]);
	TestEntity te = input_csv_test_seq(id, fn);
	fprintf(stdout, "motif=%s, period=%d, mutation=%d, flank=(%d,%d)\n",
		 te.motif.c_str(), te.period, te.mutation, te.flank_l, te.flank_r);
	fprintf(stdout, "motif_len=%ld, seq_len=%ld\n", te.motif.length(), te.seq.length());

	solve(te.seq.length(), te.seq.c_str());
	return 0;
}