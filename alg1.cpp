#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;
#define N 1000000
#define M 10000

struct Point {
    double x, y, z;

    Point() {}

    Point(double x, double y, double z) : x(x), y(y), z(z) {}

    double &operator[](int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    Point operator-(Point p) {
        return Point(x - p.x, y - p.y, z - p.z);
    }

    Point cross(Point p) {
        return Point(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
    }
};

double x[2], y[2], z[2];
double pts[N][5][3];
int n;
int faces[6][3] = {0, 1, 2, 3, 0, 2, 4, 3, 2, 1, 4, 2, 1, 0, 4, 4, 0, 3};

double det(double a[3][3]) {
    return a[0][0] * a[1][1] * a[2][2] - a[0][0] * a[1][2] * a[2][1] +
           a[0][1] * a[1][2] * a[2][0] - a[0][1] * a[1][0] * a[2][2] +
           a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0];
}

bool planeLine(Point plane[3], Point line[2], Point &ans) {
    double a[3][3], v[2];
    for (int i = 0; i < 2; i++) {
        for (int k = 0; k < 3; k++)
            for (int p = 0; p < 3; p++)
                a[k][p] = plane[k][p] - line[i][p];
        v[i] = det(a);
    }

    for (int i = 0; i < 3; i++) ans[i] = (v[1] * line[0][i] - v[0] * line[1][i]) / (v[1] - v[0]);

    if ((v[0] > 0) == (v[1] > 0)) return false;


    Point plane_[3];
    for (int k = 0; k < 3; k++) plane_[k] = plane[k] - ans;
    Point normal[3];
    for (int k = 0; k < 3; k++) {
        normal[k] = plane_[k].cross(plane_[(k + 1) % 3]);
    }
    if ((normal[0][0] > 0) == (normal[1][0] > 0) && (normal[1][0] > 0) == (normal[2][0] > 0)) return true;
    return false;

}

int cnttt = 0;

void plane_poly(int p, int face, int q) {
    Point plane[3];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            plane[i][j] = pts[p][faces[face][i]][j];
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++) {
            Point line[2];
            for (int k = 0; k < 3; k++) {
                line[0][k] = pts[q][faces[i][j]][k];
                line[1][k] = pts[q][faces[i][(j + 1) % 3]][k];
            }
            Point ans;
            if (planeLine(plane, line, ans)) cnttt++;
        }
}

int main() {
    srand(time(0));
    freopen("/Users/mazeyu/Downloads/coverage_3D/coverage.in", "r", stdin);
    scanf("%lf%lf%lf%lf%lf%lf", &x[0], &x[1], &y[0], &y[1], &z[0], &z[1]);
    scanf("%d", &n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 5; j++)
            for (int k = 0; k < 3; k++) {
                scanf("%lf", &pts[i][j][k]);
            }
    }
//    for (int i = 0; i < n; i++)
//        for (int j = 0; j < 6; j++) {
//            cnttt = 0;
//            for (int k = 0; k < n; k++)
//                if (k != i) plane_poly(i, j, k);
//            cout << cnttt << endl;
//
//        }
    int iter = 1000 / 8;
    int cnt = 0;
    for (int it = 0; it < iter; it++) {
        double pt[3];
        pt[0] = x[0] + (x[1] - x[0]) * (rand() % M) / M;
        pt[1] = y[0] + (y[1] - y[0]) * (rand() % M) / M;
        pt[2] = z[0] + (z[1] - z[0]) * (rand() % M) / M;

        for (int i = 0; i < n; i++) {
            bool flag = true;
            for (int j = 0; j < 6; j++) {
                double a[3][3];
                for (int k = 0; k < 3; k++)
                    for (int p = 0; p < 3; p++)
                        a[k][p] = pts[i][faces[j][k]][p] - pt[p];
                double v = det(a);
                if (v < 0) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                cnt += 1;
                break;
            }
        }
    }
    printf("%lf\n", cnt * 1.0 / iter);
    freopen("/Users/mazeyu/Downloads/coverage_3D/coverage.out", "w", stdout);
    printf("%lf\n", cnt * 1.0 / iter);


}