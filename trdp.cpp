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

const int NORMAL = 0;
const int WITHIN_REP = 1;
const int START_REP = 2;
const int END_REP = 3;

struct BackTrace {
	int x, y;
	int t; // 0: normal DP; 1: within rep; 2: start of rep; 3: end of rep
	BackTrace() { x = y = t = -1; }
	BackTrace(int xx, int yy, int tt): x(xx), y(yy), t(tt) {}
};

void solve(int n, const char *seq, const string &out_fn)
{
	const int INF = 100000000;
	const int MIN_UNIT = 5; // Minimum repeat unit size
	const int match_score = 1;
	const int mismatch_penalty = -4;
	const int open_gap_penalty = -6;
	const int extend_gap_penalty = -1;
	const int open_rep = -4;
	const int close_rep = -4;
	vector<vector<int>> E;
	vector<vector<int>> F;
	vector<vector<int>> H;
	vector<vector<int>> D; // Start of duplications
	vector<vector<int>> GE, GF, GH; // Alignment of duplications
	vector<vector<BackTrace>> bt;
	vector<vector<BackTrace>> dbt;
	E.resize(n + 1);
	F.resize(n + 1);
	H.resize(n + 1);
	D.resize(n + 1);
	GE.resize(n + 1);
	GF.resize(n + 1);
	GH.resize(n + 1);
	bt.resize(n + 1);
	dbt.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		E[i].resize(i + 1, -INF);
		F[i].resize(i + 1, -INF);
		H[i].resize(i + 1, -INF);
		D[i].resize(i + 1, -INF);
		GE[i].resize(i + 1, -INF);
		GF[i].resize(i + 1, -INF);
		GH[i].resize(i + 1, -INF);
		bt[i].resize(i + 1);
		dbt[i].resize(i + 1);
	}
	// Global alignment in left-down triangle
	H[0][0] = 0;
	for (int i = 1; i <= n; i++) {
		H[i][0] = E[i][0] = open_gap_penalty + extend_gap_penalty * i;
		for (int j = 1; j < i; j++) {
			F[i][j] = max(H[i][j-1] + open_gap_penalty, F[i][j-1]) + extend_gap_penalty;
			E[i][j] = max(H[i-1][j] + open_gap_penalty, E[i-1][j]) + extend_gap_penalty;
			int M = H[i-1][j-1] + (seq[i-1] == seq[j-1] ?match_score :mismatch_penalty);
			if (E[i][j] > H[i][j]) {
				H[i][j] = E[i][j];
				bt[i][j] = BackTrace(i-1, j, NORMAL);
			}
			if (F[i][j] > H[i][j]) {
				H[i][j] = F[i][j];
				bt[i][j] = BackTrace(i, j-1, NORMAL);
			}
			if (M > H[i][j]) {
				H[i][j] = M;
				bt[i][j] = BackTrace(i-1, j-1, NORMAL);
			}
		}
		int f = max(H[i][i-1] + open_gap_penalty, F[i][i-1]) + extend_gap_penalty;
		int m = H[i-1][i-1] + match_score;
		if (f > H[i][i]) {
			H[i][i] = f;
			bt[i][i] = BackTrace(i, i-1, NORMAL);
		}
		if (m > H[i][i]) {
			H[i][i] = m;
			bt[i][i] = BackTrace(i-1, i-1, NORMAL);
		}

		if (i >= MIN_UNIT + 1) {
			// D comes from the max H value in last row
			// FIXME: is it correct to calculate D this way?
			int max_value = INT32_MIN, max_id = -1;
			for (int j = 0; j <= i-1; j++) {
				if (H[i-1][j] > max_value) {
					max_value = H[i-1][j];
					max_id = j;
				}
				if (GH[i-1][j] > max_value) {
					max_value = GH[i-1][j];
					max_id = j;
				}
			}
			for (int j = 0; j <= max_id - MIN_UNIT + 1; j++) {
				D[i][j] = max_value + open_rep + (seq[i-1] == seq[j-1] ?match_score * 2 :mismatch_penalty);
				dbt[i][j] = BackTrace(i-1, max_id, START_REP);
			}
		}

		// Alignment starting from D
		for (int j = 1; j < i; j++) {
			GE[i][j] = max(D[i-1][j] + open_gap_penalty, GE[i-1][j]) + extend_gap_penalty;
			GF[i][j] = max(D[i][j-1] + open_gap_penalty, GF[i][j-1]) + extend_gap_penalty;
			int GM = max(GH[i-1][j-1], D[i-1][j-1]) + (seq[i-1] == seq[j-1] ?match_score * 2 :mismatch_penalty);
			if (GE[i][j] > GH[i][j]) {
				GH[i][j] = GE[i][j];
				dbt[i][j] = BackTrace(i-1, j, WITHIN_REP);
			}
			if (GF[i][j] > GH[i][j]) {
				GH[i][j] = GF[i][j];
				dbt[i][j] = BackTrace(i, j-1, WITHIN_REP);
			}
			if (GM > GH[i][j]) {
				GH[i][j] = GM;
				dbt[i][j] = BackTrace(i-1, j-1, WITHIN_REP);
			}
		}
		int gf = max(GH[i][i-1] + open_gap_penalty, GF[i][i-1]) + extend_gap_penalty;
		int gm = max(GH[i-1][i-1] + match_score, D[i-1][i-1] + match_score);
		if (gf > GH[i][i]) {
			GH[i][i] = gf;
			dbt[i][i] = BackTrace(i, i-1, WITHIN_REP);
		}
		if (gm > GH[i][i]) {
			GH[i][i] = gm;
			dbt[i][i] = BackTrace(i-1, i-1, WITHIN_REP);
		}

		// B transfer
		int max_value = INT32_MIN, max_id = -1;
		for (int j = 0; j <= i; j++) {
			if (GH[i][j] > max_value) {
				max_value = GH[i][j];
				max_id = j;
			}
		}
		if (max_value + close_rep > H[i][i]) {
			H[i][i] = max_value + close_rep;
			bt[i][i] = BackTrace(i, max_id, END_REP);
		}
	}
	vector<vector<char>> mark;
	mark.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		mark[i].resize(i + 1, ' ');
	}
	int x = n, y = n;
	BackTrace k = bt[x][y];
	while (x > 0 and y > 0) {
		int px = x, py = y;
		x = k.x;
		y = k.y;
		if (k.t == NORMAL) {
			mark[px][py] = '*';
			k = bt[x][y];
		} else if (k.t == END_REP) {
			mark[px][py] = 'B';
			k = dbt[x][y];
		} else if (k.t == WITHIN_REP) {
			mark[px][py] = 'G';
			k = dbt[x][y];
		} else if (k.t == START_REP) {
			mark[px][py] = 'D';
			k = bt[x][y];
		}
	}

	ofstream out(out_fn);
	assert(out.is_open());
	for (int i = 0; i <= n; i++) {
		if (i == 0) out << "0";
		else out << "," << i << "(" << seq[i-1] << ")";
	}
	out << endl;
	for (int i = 0; i <= n; i++) {
		if (i > 0) out << seq[i-1];
		else out << 0;
		for (int j = 1; j <= i; j++) {
			out << ",";
			if (mark[i][j] == ' ')
				out << H[i][j] << mark[i][j];
			else
				out << GH[i][j] << mark[i][j];
		}
		for (int j = i+1; j <= n; j++) {
			out << ",";
		}
		out << endl;
	}
	out.close();
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

	string out_fn = "../self_" + string(argv[2]) + ".csv";
	solve(te.seq.length(), te.seq.c_str(), out_fn);
	return 0;
}