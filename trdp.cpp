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

const int INF = 100000000;

const int MIN_UNIT = 5; // Minimum repeat unit size
const int MAT_SCORE = 1;
const int RMAT_SCORE = 2; // FIXME: would it be too large?
const int MIS_PEN = -4;
const int GAP_O = -6;
const int GAP_E = -1;
const int OPEN_REP = -2; //
const int CLOSE_REP = -4; // TODO: penalize it harder

const int NORMAL = 0;
const int START_REP = 1;
const int NEW_COPY = 2;
const int WITHIN_REP = 3;
const int END_REP = 4;

struct DpCell {
	int E, F, H; // The original SW matrix
	int D_gate; // Gate of duplication
	int D_from; // Where duplication starts
	int de, df, dh; // Sub matrix of duplication
	int pi, pj, event; // Backtrace

	DpCell() {
		E = F = H = -INF;
		D_gate = -INF;
		D_from = -1;
		de = df = dh = -INF;
		pi = pj = event = -1;
	}
};

void solve(int n, const char *seq, const string &out_fn)
{
	vector<vector<DpCell>> dp;
	dp.resize(n + 1);
	for (int i = 0; i <= n; i++) {
		dp[i].resize(i + 1);
	}
	// Global alignment in left-down triangle
	dp[0][0].H = 0;
	for (int i = 1; i <= n; i++) {
		dp[i][0].E = dp[i][0].H = GAP_O + GAP_E * i;
		for (int j = 1; j < i; j++) {
			dp[i][j].E = max(dp[i-1][j].H + GAP_O, dp[i-1][j].E) + GAP_E;
			dp[i][j].F = max(dp[i][j-1].H + GAP_O, dp[i][j-1].F) + GAP_E;
			int M = dp[i-1][j-1].H + (seq[i-1] == seq[j-1] ? MAT_SCORE : MIS_PEN);
			if (dp[i][j].E > dp[i][j].H) {
				dp[i][j].H = dp[i][j].E;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j;
				dp[i][j].event = NORMAL;
			}
			if (dp[i][j].F > dp[i][j].H) {
				dp[i][j].H = dp[i][j].F;
				dp[i][j].pi = i;
				dp[i][j].pj = j-1;
				dp[i][j].event = NORMAL;
			}
			if (M > dp[i][j].H) {
				dp[i][j].H = M;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j-1;
				dp[i][j].event = NORMAL;
			}
		}
		int f = max(dp[i][i-1].H + GAP_O, dp[i][i-1].F) + GAP_E;
		int m = dp[i-1][i-1].H + MAT_SCORE;
		if (f > dp[i][i].H) {
			dp[i][i].H = f;
			dp[i][i].pi = i;
			dp[i][i].pj = i-1;
			dp[i][i].event = NORMAL;
		}
		if (m > dp[i][i].H) {
			dp[i][i].H = m;
			dp[i][i].pi = i-1;
			dp[i][i].pj = i-1;
			dp[i][i].event = NORMAL;
		}

		if (i > MIN_UNIT) {
			// D gate comes from the last row
			int max_value = dp[i-1][i-1].H; // Diagonal must be the maximum in regular matrix
			int from = i-1;
			int event = START_REP;
			for (int j = 1; j < i-1; j++) {
				if (dp[i-1][j].dh > max_value) { // If multiple maximums exist, choose the first one
					max_value = dp[i-1][j].dh;
					from = j;
					event = NEW_COPY;
				}
			}
			// 0 is excluded because a match/mismatch is mandatory
			for (int j = 1; j <= from - MIN_UNIT + 1; j++) {
				// No penalty for new copy
				int tmp = (seq[i-1] == seq[j-1] ?RMAT_SCORE :MIS_PEN) + (event == START_REP ?OPEN_REP :0);
				dp[i][j].D_gate = max_value + tmp;
				dp[i][j].D_from = from;
				dp[i][j].pi = i-1;
				dp[i][j].pj = from;
				dp[i][j].event = event;
			}
		}

		// Repetition alignment starting from D gates
		for (int j = 1; j < i; j++) {
			int v_score = max(max(dp[i-1][j].D_gate, dp[i-1][j].dh) + GAP_O, dp[i-1][j].de) + GAP_E;
			int h_score = -INF;
			int d_score = -INF;
			int tmp = (seq[i-1] == seq[j-1] ?RMAT_SCORE :MIS_PEN);
			if (j <= dp[i][j-1].D_from) {
				h_score = max(max(dp[i][j-1].D_gate, dp[i][j-1].dh) + GAP_O, dp[i][j-1].df) + GAP_E;
				d_score = max(dp[i-1][j-1].D_gate, dp[i-1][j-1].dh) + tmp;
			}
			if (v_score > dp[i][j].dh) {
				dp[i][j].dh = v_score;
				dp[i][j].D_from = dp[i-1][j].D_from;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j;
				dp[i][j].event = WITHIN_REP;
			}
			if (h_score > dp[i][j].dh) {
				dp[i][j].dh = h_score;
				dp[i][j].D_from = dp[i][j-1].D_from;
				dp[i][j].pi = i;
				dp[i][j].pj = j-1;
				dp[i][j].event = WITHIN_REP;
			}
			if (d_score > dp[i][j].dh) {
				dp[i][j].dh = d_score;
				dp[i][j].D_from = dp[i-1][j-1].D_from;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j-1;
				dp[i][j].event = WITHIN_REP;
			}
		}
		// The alignment above can't reach the diagonal

		// B transfer: close a repetition
		int max_value = -INF, max_j = -1;
		for (int j = 1; j < i; j++) {
			if (dp[i][j].dh > max_value) {
				max_value = dp[i][j].dh;
				max_j = j;
			}
		}
		if (max_value + CLOSE_REP > dp[i][i].H) {
			dp[i][i].H = max_value + CLOSE_REP;
			dp[i][i].pi = i;
			dp[i][i].pj = max_j;
			dp[i][i].event = END_REP;
		}
	}

	// Trace back the optimal path
	int ptr_i = n, ptr_j = n;
	while (ptr_i > 0 and ptr_j > 0) {
		const DpCell &t = dp[ptr_i][ptr_j];
		if (t.event == END_REP) {
			assert(ptr_i == ptr_j); // Only main diagonal closes repetitions
			fprintf(stderr, "(%d,%d) -> Diagonal(%d)\n", t.pi, t.pj, ptr_i);
		} else if (t.event == START_REP) {
			fprintf(stderr, "Diagonal(%d) -> (%d, %d)\n", t.D_from, ptr_i, ptr_j);
		} else if (t.event == NEW_COPY) {
			fprintf(stderr, "New copy found at (%d, %d)\n", ptr_i, ptr_j);
		}
		ptr_i = t.pi;
		ptr_j = t.pj;
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
			switch (dp[i][j].event) {
				case NORMAL:
					out << dp[i][j].H;
					break;
				case START_REP:
					out << dp[i][j].D_gate << " D " << dp[i][j].D_from;
					break;
				case NEW_COPY:
					out << dp[i][j].dh << " C " << dp[i][j].D_from;
					break;
				case WITHIN_REP:
					out << dp[i][j].dh << " W ";
					break;
				case END_REP:
					out << dp[i][j].H << " B " << dp[i][j].pj;
					break;
				default:
					abort();
			}
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