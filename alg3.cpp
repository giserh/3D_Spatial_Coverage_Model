#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;
#define N 100000
#define M 10000
double x[2], y[2], z[2];
double pts[N][5][3];
int sections[N];
int n;
int faces[6][3] = {0, 1, 2, 3, 0, 2, 4, 3, 2, 1, 4, 2, 1, 0, 4, 4, 0, 3};
int cells;
double weights[32][32][32];

double pts_[8][N][5][3];
int pts_cnt[8];



double det(double a[3][3]) {
    return a[0][0] * a[1][1] * a[2][2] - a[0][0] * a[1][2] * a[2][1] +
           a[0][1] * a[1][2] * a[2][0] - a[0][1] * a[1][0] * a[2][2] +
           a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0];
}

inline double sqr(double x) {
    return x * x;
}

int main() {
    srand(time(0));
    freopen("../Python code/coverage.in", "r", stdin);
    scanf("%lf%lf%lf%lf%lf%lf", &x[0], &x[1], &y[0], &y[1], &z[0], &z[1]);

    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 3; k++) {
                scanf("%lf", &pts[i][j][k]);
            }
    }
    for (int i = 0; i < n; i++) {
        scanf("%d", &sections[i]);
    }
    scanf("%d", &cells);
    double Rv = 0.05;

    for (int i = 0; i < cells; i++)
        for (int j  = 0; j < cells; j++)
            for (int k = 0; k < cells; k++) {
                double center[3];
                center[0] = x[0] + (x[1] - x[0]) * (i + 0.5) / cells;
                center[1] = y[0] + (y[1] - y[0]) * (j + 0.5) / cells;
                center[2] = z[0] + (z[1] - z[0]) * (k + 0.5) / cells;
                double sum = 0;
                int cnt = 0;
                for (int in = 0; in < n; in++)
                {
                    double dis = sqrt(sqr(center[0] - pts[in][2][0]) + sqr(center[1] - pts[in][2][1])
                            + sqr(center[2] - pts[in][2][2]));
                    if (dis < Rv) {
                        sum += dis;
                        cnt++;
                    }
                }


                weights[i][j][k] = cnt? (1 - sum / cnt / Rv): 1;


            }


    for (int i = 0; i < n; i++)
    {

        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 3; k++)
                pts_[sections[i]][pts_cnt[sections[i]]][j][k] = pts[i][j][k];
        pts_cnt[sections[i]] ++;

    }
    double cnt = 0;
    int iter = 1000.0 / cells / cells / cells / 8 + 1;

    for (int i = 0; i < cells; i++)
        for (int j = 0; j < cells; j++)
            for (int k = 0; k < cells; k++) {
                double cnt_ = 0;
                double x_[2], y_[2], z_[2], size = (x[1] - x[0]) / cells;
                x_[0] = x[0] + (x[1] - x[0]) * i / cells;
                x_[1] = x_[0] + size;
                y_[0] = y[0] + (y[1] - y[0]) * j / cells;
                y_[1] = y_[0] + size;
                z_[0] = z[0] + (z[1] - z[0]) * k / cells;
                z_[1] = z_[0] + size;




                for (int angles = 0; angles < 8; angles++) {

                    for (int it = 0; it < iter; it++) {
                        double pt[3];
                        pt[0] = x_[0] + (x_[1] - x_[0]) * (rand() % M) / M;
                        pt[1] = y_[0] + (y_[1] - y_[0]) * (rand() % M) / M;
                        pt[2] = z_[0] + (z_[1] - z_[0]) * (rand() % M) / M;

                        for (int in = 0; in < pts_cnt[angles]; in++) {
                            bool flag = true;
                            for (int jn = 0; jn < 6; jn++) {
                                double a[3][3];
                                for (int kn = 0; kn < 3; kn++)
                                    for (int p = 0; p < 3; p++)
                                        a[kn][p] = pts_[angles][in][faces[jn][kn]][p] - pt[p];
                                double v = det(a);
                                if (v < 0) {
                                    flag = false;
                                    break;
                                }
                            }
                            if (flag) {
                                cnt_++;
                                break;
                            }
                        }
                    }

                }
                cnt += cnt_ * weights[i][j][k];

            }
    double ans1 = cnt * 1.0 / iter / 8 / cells / cells / cells;

    printf("%lf\n", ans1);
    freopen("../Python code/coverage.out", "w", stdout);
    printf("%lf\n", ans1);



}