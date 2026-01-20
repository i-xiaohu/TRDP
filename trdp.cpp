//
// Created by ixiaohu on 2025/12/30.
//

#include <vector>
#include <cstdint>
#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <iostream>
#include <getopt.h>
#include <algorithm>
#include "utils.h"
using namespace std;


// TRDP parameters
struct TrdpOptions {
	// Minimum repeat unit size
	int min_unit_size;
	// Scoring matrix for self-alignment
	int mat_score;
	int mis_pen;
	int gap_o;
	int gap_e;
	// Scoring matrix for duplication (reward more and/or penalize less than normal matrix)
	int rep_mat_score;
	int rep_mis_pen;
	int rep_gap_o;
	int rep_gap_e;
	// It is necessary to penalize
	int open_rep_pen;
	int close_rep_pen;
	const char *vis_fn;

	TrdpOptions() {
		min_unit_size = 5;
		mat_score = 1;
		mis_pen = -4;
		gap_o = -6;
		gap_e = -1;
		// TODO: would it be too lenient for sub-matrix?
		rep_mat_score = 2;
		rep_mis_pen = -3;
		rep_gap_o = -3;
		rep_gap_e = -1;
		open_rep_pen = -2;
		close_rep_pen = -6; // TODO: how to set this value?
		vis_fn = nullptr;
	}
};

const int INF = 100000000;
const int NORMAL = 0;
const int START_REP = 1;
const int NEW_COPY = 2;
const int WITHIN_REP = 3;
const int END_REP = 4;

int DEBUG = 0;

struct DpCell {
	int E, F, H; // The original SW matrix
	int D_gate; // Gate of duplication
	int D_from; // Where duplication starts
	int C_from; // Where new copy starts
	int de, df, dh; // Sub matrix of duplication
	int pi, pj, event; // Backtrace

	DpCell() {
		E = F = H = -INF;
		D_gate = -INF;
		D_from = -1;
		C_from = -1;
		de = df = dh = -INF;
		pi = pj = event = -1;
	}
};

struct RepUnit {
	int tb, te;
	int qb, qe; // [qb, qe) is a tandem repeat of [tb, te)
	int score; // Alignment score
	int match, mis, gap;

	RepUnit() {
		tb = te = 0;
		qb = qe = 0;
		score = -INF;
		match = mis = gap = 0;
	}
};

// Stage 1: identify breakpoints of tandem repeats using self alignment
vector<RepUnit> self_alignment(const TrdpOptions &o, int n, const char *seq)
{
	double ctime = cputime();
	const int MAT_SCORE = o.mat_score;
	const int MIS_PEN = o.mis_pen;
	const int GAP_O = o.gap_o;
	const int GAP_E = o.gap_e;
	const int MIN_UNIT = o.min_unit_size;
	const int OPEN_REP = o.open_rep_pen;
	const int CLOSE_REP = o.close_rep_pen;
	const int RMAT_SCORE = o.rep_mat_score;
	const int RMIS_PEN = o.rep_mis_pen;
	const int RGAP_O = o.rep_gap_o;
	const int RGAP_E = o.rep_gap_e;

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
			int D_from = i - 1;
			int C_from = -1;
			int event = START_REP;
			for (int j = 1; j < i-1; j++) {
				// FIXME: This prevents a small sub-matrix from forming, which might have better score in global view
				if (dp[i-1][j].dh > max_value and dp[i-1][j].D_from == j) { // If multiple maximums exist, choose the first one
					max_value = dp[i-1][j].dh;
					C_from = dp[i-1][j].C_from;
					D_from = j;
					event = NEW_COPY;
				}
			}

			// 0 is excluded because a match/mismatch is mandatory
			if (event == START_REP) {
				// Open duplication
				for (int j = 1; j <= D_from - MIN_UNIT + 1; j++) {
					int tmp = (seq[i-1] == seq[j-1] ?RMAT_SCORE :RMIS_PEN) + OPEN_REP;
					dp[i][j].D_gate = max_value + tmp;
					dp[i][j].C_from = j;
					dp[i][j].D_from = D_from;
					dp[i][j].pi = i-1;
					dp[i][j].pj = D_from;
					dp[i][j].event = event;
				}
			} else {
				// FIXME: cells in [1, c_from] compete with d_from?
				assert(C_from > 0);
				for (int j = C_from; j <= D_from - MIN_UNIT + 1; j++) {
					// No penalty for new copy
					int tmp = (seq[i-1] == seq[j-1] ?RMAT_SCORE :RMIS_PEN);
					dp[i][j].D_gate = max_value + tmp;
					dp[i][j].C_from = C_from;
					dp[i][j].D_from = D_from;
					dp[i][j].pi = i-1;
					dp[i][j].pj = D_from;
					dp[i][j].event = event;
				}
			}
		}

		// Repetition alignment starting from D gates
		for (int j = 1; j < i; j++) {
			int v_score = max(max(dp[i-1][j].D_gate, dp[i-1][j].dh) + RGAP_O, dp[i-1][j].de) + RGAP_E;
			int h_score = -INF;
			int d_score = -INF;
			int tmp = (seq[i-1] == seq[j-1] ?RMAT_SCORE :RMIS_PEN);
			if (j <= dp[i][j-1].D_from) {
				h_score = max(max(dp[i][j-1].D_gate, dp[i][j-1].dh) + RGAP_O, dp[i][j-1].df) + RGAP_E;
				d_score = max(dp[i-1][j-1].D_gate, dp[i-1][j-1].dh) + tmp;
			}
			DpCell backup = dp[i][j];
			if (v_score > dp[i][j].dh) {
				dp[i][j].dh = v_score;
				dp[i][j].C_from = dp[i-1][j].C_from;
				dp[i][j].D_from = dp[i-1][j].D_from;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j;
				dp[i][j].event = WITHIN_REP;
			}
			if (h_score > dp[i][j].dh) {
				dp[i][j].dh = h_score;
				dp[i][j].C_from = dp[i][j-1].C_from;
				dp[i][j].D_from = dp[i][j-1].D_from;
				dp[i][j].pi = i;
				dp[i][j].pj = j-1;
				dp[i][j].event = WITHIN_REP;
			}
			if (d_score > dp[i][j].dh) {
				dp[i][j].dh = d_score;
				dp[i][j].C_from = dp[i-1][j-1].C_from;
				dp[i][j].D_from = dp[i-1][j-1].D_from;
				dp[i][j].pi = i-1;
				dp[i][j].pj = j-1;
				dp[i][j].event = WITHIN_REP;
			}
			// D gate opens again after sub-matrix calculation
			if (dp[i][j].D_gate > dp[i][j].dh) {
				dp[i][j] = backup;
			}
		}
		// The alignment above can't reach the diagonal

		// B transfer: close a repetition
		int max_value = -INF, max_j = -1;
		for (int j = 1; j < i; j++) {
			// Only return to diagonal if sub-matrix reaches the lower-right corner
			if (dp[i][j].dh > max_value and dp[i][j].D_from == j) {
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
	vector<RepUnit> repetitions;
	while (ptr_i > 0 and ptr_j > 0) {
		DpCell &t = dp[ptr_i][ptr_j];
		if (t.event == END_REP) {
			assert(ptr_i == ptr_j); // Only main diagonal closes repetitions
			if (DEBUG) fprintf(stderr, "Close: %d -> %d (Diagonal)\n", t.pj, ptr_i);

			while (t.event != START_REP) {
				ptr_i = t.pi;
				ptr_j = t.pj;
				t = dp[ptr_i][ptr_j]; // Lower-right corner of the sub-matrix
				RepUnit u;
				u.score = t.dh;
				u.qe = ptr_i + 1;
				u.te = ptr_j + 1;
				while (t.event != NEW_COPY and t.event != START_REP) {
					if (ptr_i == t.pi + 1 and ptr_j == t.pj + 1) {
						if (seq[ptr_i-1] == seq[ptr_j-1]) u.match++;
						else u.mis++;
					} else {
						u.gap++;
					}
					ptr_i = t.pi;
					ptr_j = t.pj;
					t = dp[ptr_i][ptr_j];
				}
				// Pointer is now at the upper-left corner
				u.score -= dp[ptr_i][ptr_j].D_gate;
				u.qb = ptr_i;
				u.tb = ptr_j;
				if (seq[ptr_i-1] == seq[ptr_j-1]) {
					u.match++;
					u.score += o.rep_mat_score;
				} else {
					u.mis++;
					u.score += o.rep_mis_pen;
				}
				repetitions.push_back(u);

				if (t.event == START_REP) {
					if (DEBUG) fprintf(stderr, "Open: %d (Diagonal) -> %d, unit_size=%d\n", ptr_j, t.pj, t.pj - ptr_j + 1);
				} else if (t.event == NEW_COPY) {
					if (DEBUG) fprintf(stderr, "Copy: %d (i=%d) -> %d, unit_size=%d\n", ptr_j, ptr_i, t.pj, t.pj - ptr_j + 1);
				}
			}
		}
		ptr_i = t.pi;
		ptr_j = t.pj;
	}

	reverse(repetitions.begin(), repetitions.end());
	fprintf(stderr, "Self-alignment found %ld duplications in %.3f CPU seconds\n", repetitions.size(), cputime() - ctime);
//	for (int i = 0; i < repetitions.size(); i++) {
//		const RepUnit &u = repetitions[i];
//		fprintf(stdout, "[%d] Tandem repeat between [%d,%d) and [%d,%d), score=%d (matches=%d, mismatches=%d, gaps=%d)\n",
//		  i+1, u.qb, u.qe, u.tb, u.te, u.score, u.match, u.mis, u.gap);
//	}

	if (o.vis_fn) {
		ofstream out(o.vis_fn);
		assert(out.is_open());
		int ti = n, tj = n, te = -1;
		while (ti > 0 and tj > 0) {
			const DpCell &c = dp[ti][tj];
			if (c.event != te) {
				out << ti << "\t" << tj << endl;
			}
			ti = c.pi;
			tj = c.pj;
			te = c.event;
		}
		out << 0 << "\t" << 0 << endl;
		out.close();
	}
	return repetitions;
}

int usage(const TrdpOptions &o) {
	fprintf(stderr, "Usage: TRDP [options] seq1.fa seq2.fa\n");
	fprintf(stderr, "  Scoring options:\n");
	fprintf(stderr, "    -A [INT]  match score [%d]\n", o.mat_score);
	fprintf(stderr, "    -B [INT]  mismatch penalty [%d]\n", o.mis_pen);
	fprintf(stderr, "    -O [INT]  open gap(indel) penalty [%d]\n", o.gap_o);
	fprintf(stderr, "    -E [INT]  extend gap penalty [%d]\n", o.gap_e);
	fprintf(stderr, "  Duplication scoring options in self-alignment:\n");
	fprintf(stderr, "    -u [INT]  minimum repeat unit size [%d]\n", o.min_unit_size);
	fprintf(stderr, "    -d [INT]  open duplication penalty [%d]\n", o.open_rep_pen);
	fprintf(stderr, "    -p [INT]  close duplication penalty [%d]\n", o.close_rep_pen);
	fprintf(stderr, "    -a [INT]  match score [%d]\n", o.rep_mat_score);
	fprintf(stderr, "    -b [INT]  mismatch penalty [%d]\n", o.rep_mis_pen);
	fprintf(stderr, "    -o [INT]  open gap(indel) penalty [%d]\n", o.rep_gap_o);
	fprintf(stderr, "    -e [INT]  extend gap penalty [%d]\n", o.rep_gap_e);
	fprintf(stderr, "Note: duplication scoring matrix must reward more and/or \n"
				    "  penalize less than self-alignment matrix to drive out duplication events.\n");
	return 1;
}

void compare_tr_seqs(const TrdpOptions &opt, const char *fn1, const char *fn2) {
	string seq1 = input_fasta_seq(fn1);
	string seq2 = input_fasta_seq(fn2);
	vector<RepUnit> rep1 = self_alignment(opt, seq1.length(), seq1.data());
	vector<RepUnit> rep2 = self_alignment(opt, seq2.length(), seq2.data());
	for (int i = 1; i < rep1.size(); i++) {
		assert(rep1[i].qe >= rep1[i].te and rep1[i].qe > rep1[i-1].qe);
	}
	for (int i = 1; i < rep2.size(); i++) {
		assert(rep2[i].qe >= rep2[i].te and rep2[i].qe > rep2[i-1].qe);
	}
}

int main(int argc, char *argv[]) {
	double ctime = cputime(), rtime = realtime();
	int c;
	TrdpOptions opt;
	while ((c = getopt(argc, argv, "A:B:O:E:u:p:a:b:o:e:v:")) >= 0) {
		switch (c) {
			case 'A':
				opt.mat_score = str2int(optarg);
				break;
			case 'B':
				opt.mis_pen = str2int(optarg);
				break;
			case 'O':
				opt.gap_o = str2int(optarg);
				break;
			case 'E':
				opt.gap_e = str2int(optarg);
				break;
			case 'u':
				opt.min_unit_size = str2int(optarg);
				break;
			case 'd':
				opt.open_rep_pen = str2int(optarg);
				break;
			case 'p':
				opt.close_rep_pen = str2int(optarg);
				break;
			case 'a':
				opt.rep_mat_score = str2int(optarg);
				break;
			case 'b':
				opt.rep_mis_pen = str2int(optarg);
				break;
			case 'o':
				opt.rep_gap_o = str2int(optarg);
				break;
			case 'e':
				opt.rep_gap_e = str2int(optarg);
				break;
			case 'v':
				opt.vis_fn = optarg;
				break;
			default:
				fprintf(stderr, "Unrecognized option `%c`\n", c);
				return 1;
		}
	}

	if (argc - optind == 2) {
		compare_tr_seqs(opt, argv[optind], argv[optind+1]);
	} else {
		fprintf(stderr, "Two FASTA files are required\n");
		return 1;
	}
	fprintf(stderr, "Program finishes in %.3f CPU seconds, %.3f real seconds\n", cputime()-ctime, realtime()-rtime);
	return 0;
}