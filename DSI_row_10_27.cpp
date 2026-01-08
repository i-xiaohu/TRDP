#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <cstdint>
#include <cstdio>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
#include <iostream>
#include <getopt.h>
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

typedef struct {
    int w; // 0: 向下, 1: 斜下, 2: 向右, 3: D走位
    int last_i;
    int last_j;
    int D_flag;
    int dx;
    int gx, gy;
    int last_D_i;
    int last_D_j;
    int last_DD_i;
    int last_DD_j;
    int last_non_DD_i;
    int last_non_DD_j;
    int v; // 累计积分
} Node;

int tran(char a) {
    if (a == 'A') return 0;
    if (a == 'G') return 1;
    if (a == 'C') return 2;
    if (a == 'T') return 3;
    return -1;
}

int main(int argc, char *argv[]) {
	const char *fn = argv[optind];
	int id = atoi(argv[optind + 1]);
	TestEntity te = input_csv_test_seq(id, fn);
	fprintf(stdout, "motif=%s, period=%d, mutation=%d, flank=(%d,%d)\n",
	        te.motif.c_str(), te.period, te.mutation, te.flank_l, te.flank_r);
	fprintf(stdout, "motif_len=%ld, seq_len=%ld\n", te.motif.length(), te.seq.length());
	const char *s1_buf = te.seq.c_str();
	const char *s2_buf = te.seq.c_str();

    /* --- 输入序列（这里用固定字符串，若要从 stdin 读可改为 fgets/scanf） --- */
//    char s1_buf[] = "ATCGGAATCGGAATCGGAGGTTTACCCTAAGGTTTACCCTAA";
//    char s2_buf[] = "ATCGGAATCGGAATCGGAGGTTTACCCTAAGGTTTACCCTAA";
    /* 也可以用 scanf 从键盘读入：
       char s1_buf[1000], s2_buf[1000];
       scanf("%999s %999s", s1_buf, s2_buf);
    */

    int m = (int)strlen(s1_buf);
    int n = (int)strlen(s2_buf);

    /* 评分矩阵与间隙/dup 分数 */
    int s[4][4];
    int gap = -30;
    int dup = 30; /* 注意：原代码未使用 dup 变量，但保留 */

    /* 初始化对角、非对角 */
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            s[i][j] = -30; /* 默认 mismatch -30 */

    s[0][0] = 50; s[1][1] = 50; s[2][2] = 50; s[3][3] = 50;

    /* 保证对称（已经对称，但保留此步以防修改） */
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            if (s[i][j] != 0)
                s[j][i] = s[i][j];

    /* dp 数组大小（原程序使用 100x100），这里同样固定上限 100。
       若序列更长，请调大 MAXN */
    #define MAXN 1000
    static Node dp[MAXN][MAXN];
    static int pathm[MAXN][MAXN];

    if (m+1 > MAXN || n+1 > MAXN) {
        fprintf(stderr, "Error: sequence too long for MAXN=%d\n", MAXN);
        return 1;
    }

    /* 初始化 dp 全域 */
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            dp[i][j].v = 0;
            dp[i][j].w = 0;
            dp[i][j].D_flag = 0;
            dp[i][j].last_i = -1;
            dp[i][j].last_j = -1;
            dp[i][j].last_D_i = -1;
            dp[i][j].last_D_j = -1;
            dp[i][j].last_DD_i = -1;
            dp[i][j].last_DD_j = -1;
            dp[i][j].last_non_DD_i = -1;
            dp[i][j].last_non_DD_j = -1;
        }
    }

    dp[0][0].v = 0;
    dp[0][0].w = 0;
    dp[0][0].D_flag = 0;
    dp[0][0].last_D_i = -1;
    dp[0][0].last_D_j = -1;
    dp[0][0].last_DD_i = -1;
    dp[0][0].last_DD_j = -1;
    dp[0][0].last_non_DD_i = -1;
    dp[0][0].last_non_DD_j = -1;

    /* 初始化第一行第一列 */
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (i == 0 && j != 0) {
                dp[i][j].v = dp[i][j-1].v + gap;
                dp[i][j].w = 2;
                dp[i][j].D_flag = 0;
                dp[i][j].last_D_i = -1;
                dp[i][j].last_D_j = -1;
                dp[i][j].last_DD_i = -1;
                dp[i][j].last_DD_j = -1;
                dp[i][j].last_non_DD_i = 0;
                dp[i][j].last_non_DD_j = 0;
            }
            if (j == 0 && i != 0) {
                dp[i][j].v = dp[i-1][j].v + gap;
                dp[i][j].w = 0;
                dp[i][j].D_flag = 0;
                dp[i][j].last_D_i = -1;
                dp[i][j].last_D_j = -1;
                dp[i][j].last_DD_i = -1;
                dp[i][j].last_DD_j = -1;
                dp[i][j].last_non_DD_i = 0;
                dp[i][j].last_non_DD_j = 0;
            }
        }
    }

    /* 核心动态规划 */
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int x = tran(s1_buf[i-1]);
            int y = tran(s2_buf[j-1]);

            int score_0 = dp[i-1][j].v + gap;
            int score_1 = dp[i-1][j-1].v + s[x][y];
            if (dp[i-1][j-1].D_flag == 1 && x == y && j <= dp[i-1][j-1].last_D_j) score_1 += 1;
            else if (dp[i-1][j-1].D_flag == 1 && x == y && j > dp[i-1][j-1].last_D_j) score_1 -= 1;

            int score_2 = dp[i][j-1].v + gap;
            if (dp[i][j-1].D_flag == 1 && j > dp[i][j-1].last_D_j) score_2 = dp[i][j-1].v + 0;

            /* 先比较 0 和 1 */
            if (score_0 > score_1) {
                dp[i][j].v = score_0;
                dp[i][j].w = 0;
                dp[i][j].last_i = i - 1;
                dp[i][j].last_j = j;
                dp[i][j].D_flag = dp[i-1][j].D_flag;
                dp[i][j].last_D_i = dp[i-1][j].last_D_i;
                dp[i][j].last_D_j = dp[i-1][j].last_D_j;
                dp[i][j].last_DD_i = dp[i-1][j].last_DD_i;
                dp[i][j].last_DD_j = dp[i-1][j].last_DD_j;
                dp[i][j].last_non_DD_i = dp[i-1][j].last_non_DD_i;
                dp[i][j].last_non_DD_j = dp[i-1][j].last_non_DD_j;
            } else {
                dp[i][j].v = score_1;
                dp[i][j].w = 1;
                dp[i][j].last_i = i - 1;
                dp[i][j].last_j = j - 1;
                dp[i][j].D_flag = dp[i-1][j-1].D_flag;
                dp[i][j].last_D_i = dp[i-1][j-1].last_D_i;
                dp[i][j].last_D_j = dp[i-1][j-1].last_D_j;
                dp[i][j].last_DD_i = dp[i-1][j-1].last_DD_i;
                dp[i][j].last_DD_j = dp[i-1][j-1].last_DD_j;
                dp[i][j].last_non_DD_i = dp[i-1][j-1].last_non_DD_i;
                dp[i][j].last_non_DD_j = dp[i-1][j-1].last_non_DD_j;
            }

            /* 再比较 2 */
            if (dp[i][j].v < score_2) {
                dp[i][j].w = 2;
                dp[i][j].v = score_2;
                dp[i][j].last_i = i;
                dp[i][j].last_j = j - 1;
                dp[i][j].D_flag = dp[i][j-1].D_flag;
                dp[i][j].last_D_i = dp[i][j-1].last_D_i;
                dp[i][j].last_D_j = dp[i][j-1].last_D_j;
                dp[i][j].last_DD_i = dp[i][j-1].last_DD_i;
                dp[i][j].last_DD_j = dp[i][j-1].last_DD_j;
                dp[i][j].last_non_DD_i = dp[i][j-1].last_non_DD_i;
                dp[i][j].last_non_DD_j = dp[i][j-1].last_non_DD_j;
            }

            /* 检查可能的 D 走位（遍历 k） */
            for (int k = j + 1; k <= n; ++k) {
                if (dp[i-1][k].w == 3) continue;
                int ver_dist, hor_dist;
                if (dp[i-1][k].D_flag == 1) {
                    ver_dist = (i-1) - dp[i-1][k].last_DD_i;
                    hor_dist = k - dp[i-1][k].last_DD_j;
                } else {
                    ver_dist = (i-1) - dp[i-1][k].last_non_DD_i;
                    hor_dist = k - dp[i-1][k].last_non_DD_j;
                }
                int dist = ver_dist <= hor_dist ? ver_dist : hor_dist;

                /* 计算 score_dup */
                int score_dup = dp[i-1][k].v + s[x][y];
                if (dp[i][j].v < score_dup) {
                    dp[i][j].w = 3;
                    dp[i][j].v = score_dup;
                    dp[i][j].last_i = i - 1;
                    dp[i][j].last_j = k;
                    dp[i][j].D_flag = 1;
                    dp[i][j].last_D_i = i - 1;
                    dp[i][j].last_D_j = k;
                    dp[i][j].last_DD_i = i;
                    dp[i][j].last_DD_j = j;
                    dp[i][j].last_non_DD_i = dp[i-1][k].last_non_DD_i;
                    dp[i][j].last_non_DD_j = dp[i-1][k].last_non_DD_j;
                }
            }

            /* 回到对角线时复位 D_flag 等信息 */
            if (i == j && dp[i][j].last_i != dp[i][j].last_j) {
                dp[i][j].D_flag = 0;
                dp[i][j].last_D_i = -1;
                dp[i][j].last_D_j = -1;
                dp[i][j].last_DD_i = -1;
                dp[i][j].last_DD_j = -1;
                dp[i][j].last_non_DD_i = i;
                dp[i][j].last_non_DD_j = j;
            }
        }
    }

//    /* 输出最终得分矩阵（数值） */
//    printf("最终得分矩阵为\n");
//    for (int i = 0; i <= m; ++i) {
//        for (int j = 0; j <= n; ++j) {
//            printf("%d ", dp[i][j].v);
//        }
//        printf("\n");
//    }
//
//    /* 输出路径标记 w */
//    printf("输出标记的路径（w 值）:\n");
//    for (int i = 0; i <= m; ++i) {
//        for (int j = 0; j <= n; ++j) {
//            printf("%d ", dp[i][j].w);
//        }
//        printf("\n");
//    }

    /* 找最大累计分值的位置（原代码遍历所有格） */
    int temp = dp[0][0].v;
    int i_max = 0, j_max = 0;
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= n; ++j) {
            if (dp[i][j].v > temp) {
                temp = dp[i][j].v;
                i_max = i;
                j_max = j;
            }
        }
    }
    printf("最大累计分值位置: i_max=%d, j_max=%d, value=%d\n", i_max, j_max, temp);

    /* 回溯标记路径（从矩阵右下角回溯，原代码把 i_max 改为 m, j_max 改为 n，保留此逻辑） */
    for (int i = 0; i <= m; ++i)
        for (int j = 0; j <= n; ++j)
            pathm[i][j] = 0;

    i_max = m;
    j_max = n;
    int ii = i_max;
    int jj = j_max;
    while (ii >= 0 && jj >= 0 && dp[ii][jj].v != 0) {
        pathm[ii][jj] = 1;
        if (dp[ii][jj].w == 1) {
            ii -= 1; jj -= 1;
        } else if (dp[ii][jj].w == 0) {
            ii -= 1;
        } else if (dp[ii][jj].w == 2) {
            jj -= 1;
        } else if (dp[ii][jj].w == 3) {
            int new_i = dp[ii][jj].last_i;
            int new_j = dp[ii][jj].last_j;
            ii = new_i;
            jj = new_j;
        } else {
            fprintf(stderr, "ERROR! Only 0,1,2,3 allowed\n");
            break;
        }
    }

    /* 打印 pathm */
//    printf("pathm matrix:\n");
//    for (int i = 0; i <= m; ++i) {
//        for (int j = 0; j <= n; ++j) {
//            printf("%d ", pathm[i][j]);
//        }
//        printf("\n");
//    }

    /* 打印带星号的矩阵到屏幕 */
//    printf("最终得分矩阵（含*表示路径）:\n\t");
//    for (int j = 1; j <= n; ++j) printf("%c\t", s2_buf[j-1]);
//    printf("\n");
//    for (int i = 0; i <= m; ++i) {
//        if (i == 0) printf("\t");
//        else printf("%c\t", s1_buf[i-1]);
//        for (int j = 0; j <= n; ++j) {
//            if (pathm[i][j]) printf("%d*", dp[i][j].v);
//            else printf("%d", dp[i][j].v);
//            printf("\t");
//        }
//        printf("\n");
//    }

    /* 输出 CSV 文件（final.csv）*/
    FILE *fout = fopen("../final.csv", "w");
    if (fout) {
        fprintf(fout, ",");
        for (int j = 1; j <= n; ++j) {
            if (j > 1) fprintf(fout, ",");
            fprintf(fout, "%c", s2_buf[j-1]);
        }
        fprintf(fout, "\n");
        for (int i = 0; i <= m; ++i) {
            if (i == 0) fprintf(fout, ",");
            else fprintf(fout, "%c,", s1_buf[i-1]);
            for (int j = 1; j <= n; ++j) {
                if (pathm[i][j]) fprintf(fout, "%d*", dp[i][j].v);
                else fprintf(fout, "%d", dp[i][j].v);
                if (j != n) fprintf(fout, ",");
            }
            fprintf(fout, "\n");
        }
        fclose(fout);
//        printf("Wrote final.csv\n");
    } else {
        fprintf(stderr, "Failed to open final.csv for writing\n");
    }

    return 0;
}
