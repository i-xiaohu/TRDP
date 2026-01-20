//
// Created by ixiaohu on 2026/1/20.
//

#include <fstream>
#include <cassert>
#include <sstream>
#include "utils.h"

double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

string input_fasta_seq(const char *fn)
{
	gzFile f = gzopen(fn, "r");
	assert(f != nullptr);
	kseq_t *ks = kseq_init(f);
	// Only process the first sequence
	string ret;
	if (kseq_read(ks) >= 0) {
		fprintf(stderr, "Input sequence %s (length=%ld) from `%s`\n", ks->name.s, ks->seq.l, fn);
		int len = ks->seq.l;
		char *seq = ks->seq.s;
		ret.resize(len);
		for (int i = 0; i < len; i++) {
			seq[i] = (char)toupper(seq[i]);
			if (seq[i] !=  'A' and seq[i] != 'C' and seq[i] != 'G' and seq[i] != 'T') {
				fprintf(stderr, "Unsupported character `%c`\n", seq[i]);
				abort();
			}
			ret[i] = seq[i];
		}
	} else {
		fprintf(stderr, "No sequence found in `%s`\n", fn);
		abort();
	}
	kseq_destroy(ks);
	gzclose(f);
	return ret;
}

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
