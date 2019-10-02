//
// Created by 马泽余 on 2019-07-26.
//

#include "database.h"

#define SUBSETTING // This MACRO is for test of the effect of subsetting

//#define KDTREE // This MACRO is to enable KD-tree, but we found that in current dataset, it doesn't have much effect because this dataset is still not big enough


#define M 10000 // divide the region into M * M * M when generating random points

#include <ctime>

// KD-tree building, independent of external library, but slower
node *build(vector<point3> &pts, vector<int> &index, int l, int r, box B) {
    node *cur = new node;
    cur->b = B;
    cur->l = l;
    cur->r = r;


    if (r >= l + 1) {
//        cout << l << " " << r << endl;
        box bound;
        for (int i = 0; i < 3; i++) {
            bound[i][0] = 1e9;
            bound[i][1] = -1e9;
        }
        for (int i = l; i < r; i++) {
            auto point = pts[index[i]];
            for (int k = 0; k < 3; k++) {
                if (point[k] < bound[k][0]) bound[k][0] = point[k];
                if (point[k] > bound[k][1]) bound[k][1] = point[k];
            }
        }
        cur->b = B = bound;

        if (fabs(B[0][1] - B[0][0]) < 1e-8 || fabs(B[1][1] - B[1][0]) < 1e-8 || fabs(B[2][1] - B[2][0]) < 1e-8) {
            cur->lc = cur->rc = NULL;
            return cur;
        }
    }

    if (r <= l + 1) {
        cur->lc = cur->rc = NULL;
        return cur;
    }


    double maxSpan = 0;
    int maxDim = -1;
    for (int i = 0; i < 3; i++) {
        if (B[i][1] - B[i][0] > maxSpan) {
            maxSpan = B[i][1] - B[i][0];
            maxDim = i;
        }
    }
    box b1 = B, b2 = B;
    b2[maxDim][0] = b1[maxDim][1] = (B[maxDim][0] + B[maxDim][1]) / 2;
    int p = l, q = r - 1;
    while (p <= q) {
        while (p < r && pts[index[p]][maxDim] < b1[maxDim][1]) p++;
        while (q >= l && pts[index[q]][maxDim] > b1[maxDim][1]) q--;
        if (p < q) {
            swap(index[p], index[q]);
            p++;
            q--;
        }
    }
    cur->dim = maxDim;
    cur->lc = build(pts, index, l, p, b1);
    cur->rc = build(pts, index, p, r, b2);

    return cur;
}


// class database, which host all data and perform calculations
database::database(string dataFile, double R, double alpha, int st, int ed, bool extrinsic) : Rv(R), alpha(alpha) {
    FILE *fp;
    fp = fopen(dataFile.c_str(), "r");
    char x[100];
    fscanf(fp, "%s", x);

    cnt = 0;

    int start = clock();
    double lat_, lng_, hgt_, yaw_, pitch_, roll_, timestamp_;

    while (fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &lat_, &lng_, &hgt_, &yaw_, &pitch_,
                  &roll_, &timestamp_) != -1) {
        if (cnt >= st) {
            lat.push_back(lat_);
            lng.push_back(lng_);
            hgt.push_back(hgt_);
            yaw.push_back(yaw_);
            pitch.push_back(pitch_);
            roll.push_back(roll_);
            timestamp.push_back(timestamp_);
        }

//        cout << lat << " " << lng << " " << hgt << " " << yaw << " " << pitch << " " << roll << " " << t << endl;
        cnt++;
        if (cnt == ed) break;
    }

    cnt = lat.size();
    cout << "Load data: " << cnt << endl;
    lat_min = 1e9;
    lat_max = -1e9;
    lng_min = 1e9;
    lng_max = -1e9;
    hgt_min = 1e9;
    hgt_max = -1e9;

    for (int i = 0; i < cnt; i++) {
        if (lat[i] < lat_min) lat_min = lat[i];
        if (lat[i] > lat_max) lat_max = lat[i];
        if (lng[i] < lng_min) lng_min = lng[i];
        if (lng[i] > lng_max) lng_max = lng[i];
        if (hgt[i] < hgt_min) hgt_min = hgt[i];
        if (hgt[i] > hgt_max) hgt_max = hgt[i];
    }

    cout << "mean latitude: " << (lat_min + lat_max) / 2 << endl;
    cout << "mean longitude: " << (lng_min + lng_max) / 2 << endl;

    lat_scale = 40075.017 / 360 * 1000;
    lng_scale = 40075.017 / 360 * 1000 * cos((lat_min + lat_max) / 2 * acos(-1) / 180);
    hgt_scale = 1;
    x_center = (lng_max + lng_min) / 2 * lng_scale;
    y_center = (lat_max + lat_min) / 2 * lat_scale;
    z_center = (hgt_max + hgt_min) / 2 * hgt_scale;

//    cout << lng_scale << " " << lat_scale << endl;
//    cout << x_center << " " << y_center << " " << z_center;
    for (int i = 0; i < 3; i++) {
        bound[i][0] = 1e9;
        bound[i][1] = -1e9;
    }
    for (int i = 0; i < cnt; i++) {
        points.push_back(points53());
        bds.push_back(box());
        for (int k = 0; k < 3; k++) points[i][0][k] = 0;
        double phi = yaw[i], theta = pitch[i], psi = roll[i];
        if (!extrinsic) {
            phi = -phi;
            theta = -theta;
            psi = -psi;
        }
        double Rot[3][3];
        Rot[0][0] = cos(theta) * cos(phi);
        Rot[0][1] = cos(theta) * sin(phi);
        Rot[0][2] = -sin(theta);
        Rot[1][0] = sin(psi) * sin(theta) * cos(phi) - cos(psi) * sin(phi);
        Rot[1][1] = sin(psi) * sin(theta) * sin(phi) + cos(psi) * cos(phi);
        Rot[1][2] = sin(psi) * cos(theta);
        Rot[2][0] = cos(psi) * sin(theta) * cos(phi) + sin(psi) * sin(phi);
        Rot[2][1] = cos(psi) * sin(theta) * sin(phi) - sin(psi) * cos(phi);
        Rot[2][2] = cos(psi) * cos(theta);
        if (!extrinsic) {
            for (int j = 0; j < 3; j++)
                for (int k = j; k < 3; k++)
                    swap(Rot[j][k], Rot[k][j]);

        }
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++) {
                points[i][j * 2 + k + 1][0] = (j * 2 - 1) * R * sin(alpha / 2 * acos(-1) / 180);
                points[i][j * 2 + k + 1][1] = (k * 2 - 1) * R * sin(alpha / 2 * acos(-1) / 180);
                points[i][j * 2 + k + 1][2] = -R * cos(alpha / 2 * acos(-1) / 180);
            }
        for (int j = 0; j < 5; j++) {
            double p[3] = {0};
            for (int k = 0; k < 3; k++)
                for (int l = 0; l < 3; l++)
                    p[k] += Rot[l][k] * points[i][j][l];
            for (int k = 0; k < 3; k++) points[i][j][k] = p[k];
        }
        for (int j = 0; j < 3; j++) {
            bds[i][j][0] = 1e9;
            bds[i][j][1] = -1e9;
        }
        for (int j = 0; j < 5; j++) {
            points[i][j][0] += lng[i] * lng_scale - x_center;
            points[i][j][1] += lat[i] * lat_scale - y_center;
            points[i][j][2] += hgt[i] * hgt_scale - z_center;
            for (int k = 0; k < 3; k++) {
                if (points[i][j][k] < bound[k][0]) bound[k][0] = points[i][j][k];
                if (points[i][j][k] > bound[k][1]) bound[k][1] = points[i][j][k];
                if (points[i][j][k] < bds[i][k][0]) bds[i][k][0] = points[i][j][k];
                if (points[i][j][k] > bds[i][k][1]) bds[i][k][1] = points[i][j][k];
            }
//            cout << points[i][j][2] + z_center << " " << endl;
        }

    }

    int k = 0;
    for (auto p5: points) {
        point3 mid = p5[0];
        for (int i = 0; i < 3; i++)
            for (int j = 1; j < 5; j++)
                mid[i] += p5[j][i] / 4;
        for (int i = 0; i < 3; i++) mid[i] /= 2;
        pointOne.push_back(mid);
        ind.push_back(k++);
    }
    root = build(pointOne, ind, 0, k, bound);

    cout << "using time: " << clock() - start << endl;


//    start = 0;
//    cout << "coverage of the bounding box: " << query(bound) << endl;
//    cout << "using time: " << clock() - start << endl;

}


double det(double a[3][3]) {
    return a[0][0] * a[1][1] * a[2][2] - a[0][0] * a[1][2] * a[2][1] +
           a[0][1] * a[1][2] * a[2][0] - a[0][1] * a[1][0] * a[2][2] +
           a[0][2] * a[1][0] * a[2][1] - a[0][2] * a[1][1] * a[2][0];
}


vector<int> findInd(node *cur, box B, vector<point3> &points, vector<int> &index) {
    if (cur == NULL || cur->l == cur->r) return vector<int>();
    if (cur->lc == NULL) {
        vector<int> ans;
        for (int i = cur->l; i < cur->r; i++) ans.push_back(i);
        return ans;
    }

    for (int j = 0; j < 3; j++)
        if (cur->b[j][0] > B[j][1] || cur->b[j][1] < B[j][0])
            return vector<int>();

    bool full = 1;
    for (int i = 0; i < 3; i++) {
        full &= B[i][1] >= cur->b[i][1];
        B[i][1] = min(B[i][1], cur->b[i][1]);
    }
    for (int i = 0; i < 3; i++) {
        full &= B[i][0] <= cur->b[i][0];
        B[i][0] = max(B[i][0], cur->b[i][0]);
    }
    if (full) {
        vector<int> ans;
        for (int i = cur->l; i < cur->r; i++) ans.push_back(i);
        return ans;
    }
    auto v1 = findInd(cur->lc, B, points, index);
    auto v2 = findInd(cur->rc, B, points, index);
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
}


double database::query_alwayssubsetting(box q, int threadNum) {
    int faces[6][3] = {2, 1, 0, 4, 2, 0, 3, 4, 0, 1, 3, 0, 1, 2, 3, 3, 2, 4};


    thread threads[threadNum];
    int iter = int(1000.0 / threadNum + 0.99);
    double ans[threadNum];
    vector<int> index;

    auto start = high_resolution_clock::now();

    for (int i = 0; i < cnt; i++) {

        bool flag[3];
        for (int j = 0; j < 3; j++) flag[j] = bds[i][j][0] > q[j][1] || bds[i][j][1] < q[j][0];
        if (flag[0] || flag[1] || flag[2]) continue;
        index.push_back(i);
    }


    cout << "size: " << index.size() << endl;
//    cout << "time in subsetting: " << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << endl;

//    cout << index.size() << " of " << cnt << " intersected" << endl;
    start = high_resolution_clock::now();

    for (int t = 0; t < threadNum; t++)
        threads[t] = thread([=, &ans, &index]() {
            int hit = 0;
            for (int it = 0; it < iter; it++) {
                double pt[3];
                pt[0] = q.x[0] + (q.x[1] - q.x[0]) * (rand() % M) / M;
                pt[1] = q.y[0] + (q.y[1] - q.y[0]) * (rand() % M) / M;
                pt[2] = q.z[0] + (q.z[1] - q.z[0]) * (rand() % M) / M;

                for (int i: index) {
                    bool flag = true;
                    for (int j = 0; j < 6; j++) {
                        double a[3][3];
                        for (int k = 0; k < 3; k++)
                            for (int p = 0; p < 3; p++)
                                a[k][p] = points[i][faces[j][k]][p] - pt[p];
                        double v = det(a);
                        if (v < 0) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        hit += 1;
                        break;
                    }
                }
            }
            ans[t] = hit * 1.0 / iter;
        });
    for (int t = 0; t < threadNum; t++) threads[t].join();
    double ans_ = 0;
    for (int i = 0; i < threadNum; i++) ans_ += ans[i];
//    cout << "time in MC: " << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << endl;

    return ans_ / threadNum;

}


double database::query(box q, int threadNum) {
    int faces[6][3] = {2, 1, 0, 4, 2, 0, 3, 4, 0, 1, 3, 0, 1, 2, 3, 3, 2, 4};


    thread threads[threadNum];
    int iter = int(1000.0 / threadNum + 0.99);
    double ans[threadNum];
    vector<int> index;

    auto start = high_resolution_clock::now();

#ifdef KDTREE
    auto qq = q;
    for (int i = 0; i < 3; i++) {
        qq[i][0] -= Rv;
        qq[i][1] += Rv;
    }
    auto index_ = findInd(root, qq, pointOne, ind);
    cout << "size_: " << index_.size() << endl;
    for (int i:index_) {
        bool flag[3];
        for (int j = 0; j < 3; j++) flag[j] = bds[i][j][0] > q[j][1] || bds[i][j][1] < q[j][0];
        if (flag[0] || flag[1] || flag[2]) continue;
        index.push_back(i);
    }
//    cout << "found " << index.size() << endl;
#else

    for (int i = 0; i < cnt; i++) {
#ifdef SUBSETTING
        bool flag[3];
        for (int j = 0; j < 3; j++) flag[j] = bds[i][j][0] > q[j][1] || bds[i][j][1] < q[j][0];
        if (flag[0] || flag[1] || flag[2]) continue;
#endif
        index.push_back(i);
    }
#endif

    cout << "size: " << index.size() << endl;
//    cout << "time in subsetting: " << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << endl;

//    cout << index.size() << " of " << cnt << " intersected" << endl;
    start = high_resolution_clock::now();

    for (int t = 0; t < threadNum; t++)
        threads[t] = thread([=, &ans, &index]() {
            int hit = 0;
            for (int it = 0; it < iter; it++) {
                double pt[3];
                pt[0] = q.x[0] + (q.x[1] - q.x[0]) * (rand() % M) / M;
                pt[1] = q.y[0] + (q.y[1] - q.y[0]) * (rand() % M) / M;
                pt[2] = q.z[0] + (q.z[1] - q.z[0]) * (rand() % M) / M;

                for (int i: index) {
                    bool flag = true;
                    for (int j = 0; j < 6; j++) {
                        double a[3][3];
                        for (int k = 0; k < 3; k++)
                            for (int p = 0; p < 3; p++)
                                a[k][p] = points[i][faces[j][k]][p] - pt[p];
                        double v = det(a);
                        if (v < 0) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        hit += 1;
                        break;
                    }
                }
            }
            ans[t] = hit * 1.0 / iter;
        });
    for (int t = 0; t < threadNum; t++) threads[t].join();
    double ans_ = 0;
    for (int i = 0; i < threadNum; i++) ans_ += ans[i];
//    cout << "time in MC: " << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << endl;

    return ans_ / threadNum;

}


double database::query2(box q, int angles, int threadNum) {
    int faces[6][3] = {2, 1, 0, 4, 2, 0, 3, 4, 0, 1, 3, 0, 1, 2, 3, 3, 2, 4};


    thread threads[threadNum];
    int iter = int(1000.0 / threadNum / angles + 0.99);
    double ans[threadNum];
    vector<int> index[angles];

#ifdef KDTREE
    auto qq = q;
    for (int i = 0; i < 3; i++) {
        qq[i][0] -= Rv;
        qq[i][1] += Rv;
    }
    auto index_ = findInd(root, qq, pointOne, ind);
    for (int i: index_) {
        bool flag[3];
        for (int j = 0; j < 3; j++) flag[j] = bds[i][j][0] > q[j][1] || bds[i][j][1] < q[j][0];
        if (flag[0] || flag[1] || flag[2]) continue;
        int discAngle = int((yaw[i] + acos(-1)) / (2 * acos(-1)) * angles);
        index[discAngle].push_back(i);
    }
#else

    for (int i = 0; i < cnt; i++) {
#ifdef SUBSETTING
        bool flag[3];
        for (int j = 0; j < 3; j++) flag[j] = bds[i][j][0] > q[j][1] || bds[i][j][1] < q[j][0];
        if (flag[0] || flag[1] || flag[2]) continue;
#endif
        int discAngle = int((yaw[i] + acos(-1)) / (2 * acos(-1)) * angles);
        index[discAngle].push_back(i);
    }

#endif

    for (int t = 0; t < threadNum; t++)
        threads[t] = thread([=, &ans, &index]() {
            int hit = 0;
            for (int it = 0; it < iter; it++) {
                double pt[3];
                pt[0] = q.x[0] + (q.x[1] - q.x[0]) * (rand() % M) / M;
                pt[1] = q.y[0] + (q.y[1] - q.y[0]) * (rand() % M) / M;
                pt[2] = q.z[0] + (q.z[1] - q.z[0]) * (rand() % M) / M;

                for (int angle = 0; angle < angles; angle++)
                    for (int i: index[angle]) {
                        bool flag = true;
                        for (int j = 0; j < 6; j++) {
                            double a[3][3];
                            for (int k = 0; k < 3; k++)
                                for (int p = 0; p < 3; p++)
                                    a[k][p] = points[i][faces[j][k]][p] - pt[p];
                            double v = det(a);
                            if (v < 0) {
                                flag = false;
                                break;
                            }
                        }
                        if (flag) {
                            hit += 1;
                            break;
                        }
                    }
            }
            ans[t] = hit * 1.0 / iter / angles;
        });
    for (int t = 0; t < threadNum; t++) threads[t].join();
    double ans_ = 0;
    for (int i = 0; i < threadNum; i++) ans_ += ans[i];
    return ans_ / threadNum;

}

inline double sqr(double x) {
    return x * x;
}

double database::query3(box q, int angles, int cells, int threadNum) {
    int faces[6][3] = {2, 1, 0, 4, 2, 0, 3, 4, 0, 1, 3, 0, 1, 2, 3, 3, 2, 4};

    vector<int> index[angles];

#ifdef KDTREE
    auto qq = q;
    for (int i = 0; i < 3; i++) {
        qq[i][0] -= Rv;
        qq[i][1] += Rv;
    }
    auto index_ = findInd(root, qq, pointOne, ind);
    for (int i: index_) {
        bool flag[3];
        for (int j = 0; j < 3; j++) flag[j] = bds[i][j][0] > q[j][1] || bds[i][j][1] < q[j][0];
        if (flag[0] || flag[1] || flag[2]) continue;
        int discAngle = int((yaw[i] + acos(-1)) / (2 * acos(-1)) * angles);
        index[discAngle].push_back(i);
    }
#else

    for (int i = 0; i < cnt; i++) {
#ifdef SUBSETTING
        bool flag[3];
        for (int j = 0; j < 3; j++) flag[j] = bds[i][j][0] > q[j][1] || bds[i][j][1] < q[j][0];
        if (flag[0] || flag[1] || flag[2]) continue;
#endif
        int discAngle = int((yaw[i] + acos(-1)) / (2 * acos(-1)) * angles);
        index[discAngle].push_back(i);
    }

#endif

    double weights[16][16][16];
    for (int i = 0; i < cells; i++)
        for (int j = 0; j < cells; j++)
            for (int k = 0; k < cells; k++) {
                double center[3];
                center[0] = q.x[0] + (q.x[1] - q.x[0]) * (i + 0.5) / cells;
                center[1] = q.y[0] + (q.y[1] - q.y[0]) * (j + 0.5) / cells;
                center[2] = q.z[0] + (q.z[1] - q.z[0]) * (k + 0.5) / cells;
                double sum = 0;
                int cnt_ = 0;
                for (int angle = 0; angle < angles; angle++)
                    for (int in: index[angle]) {
                        double dis = sqrt(sqr(center[0] - points[in][0][0]) + sqr(center[1] - points[in][0][1])
                                          + sqr(center[2] - points[in][0][2]));
                        if (dis < Rv) {
                            sum += dis;
                            cnt_++;
                        }
                    }
                weights[i][j][k] = cnt_ ? (1 - sum / cnt_ / Rv) : 0;
            }
//    for (int i = 0; i < cells; i++)
//        for (int j = 0; j < cells; j++)
//            for (int k = 0; k < cells; k++)
//                cout << weights[i][j][k] << " ";
//    cout << endl;
    int iter = int(1000.0 / cells / cells / cells / angles + 0.99);
    thread threads[threadNum];
    double ans[threadNum];

    for (int t = 0; t < threadNum; t++)
        threads[t] = thread([=, &ans, &index]() {

            for (int i = 0; i < cells; i++)
                for (int j = 0; j < cells; j++)
                    for (int k = 0; k < cells; k++) {

                        if (k + j * cells + i * cells * cells < cells * cells * cells * t / threadNum) continue;
                        if ((i != cells - 1 || j != cells - 1 || k != cells - 1) &&
                            k + j * cells + i * cells * cells >= cells * cells * cells * (t + 1) / threadNum)
                            continue;

                        double x_[2], y_[2], z_[2], size = (q.x[1] - q.x[0]) / cells;
                        x_[0] = q.x[0] + (q.x[1] - q.x[0]) * i / cells;
                        x_[1] = x_[0] + size;
                        y_[0] = q.y[0] + (q.y[1] - q.y[0]) * j / cells;
                        y_[1] = y_[0] + size;
                        z_[0] = q.z[0] + (q.z[1] - q.z[0]) * k / cells;
                        z_[1] = z_[0] + size;

                        int hit = 0;
                        for (int it = 0; it < iter; it++) {
                            double pt[3];
                            pt[0] = x_[0] + (x_[1] - x_[0]) * (rand() % M) / M;
                            pt[1] = y_[0] + (y_[1] - y_[0]) * (rand() % M) / M;
                            pt[2] = z_[0] + (z_[1] - z_[0]) * (rand() % M) / M;

                            for (int angle = 0; angle < angles; angle++)
                                for (int ii: index[angle]) {
                                    bool flag = true;
                                    for (int jj = 0; jj < 6; jj++) {
                                        double a[3][3];
                                        for (int kk = 0; kk < 3; kk++)
                                            for (int p = 0; p < 3; p++)
                                                a[kk][p] = points[ii][faces[jj][kk]][p] - pt[p];
                                        double v = det(a);
                                        if (v < 0) {
                                            flag = false;
                                            break;
                                        }
                                    }
                                    if (flag) {
                                        hit++;
                                        break;
                                    }
                                }
                        }


                        ans[t] += hit * 1.0 / iter / angles * weights[i][j][k];

                    }
        });
    for (int i = 0; i < threadNum; i++) threads[i].join();
    double ans_ = 0;
    for (int i = 0; i < threadNum; i++) ans_ += ans[i];
    ans_ /= cells * cells * cells;
    return ans_;


}

void database::generateQueries(string file, int n) {
    FILE *fp;
    fp = fopen(file.c_str(), "w");
    for (int i = 0; i < n; i++) {
        double pt[3];
        auto q = bound;
        pt[0] = q.x[0] + (q.x[1] - q.x[0]) * (rand() % M) / M;
        pt[1] = q.y[0] + (q.y[1] - q.y[0]) * (rand() % M) / M;
        pt[2] = q.z[0] + (q.z[1] - q.z[0]) * (rand() % M) / M;
        fprintf(fp, "%lf %lf %lf\n", pt[0], pt[1], pt[2]);
    }

}

void database::generateContinuousQueries(string file, int n, double size) {
    FILE *fp;
    fp = fopen(file.c_str(), "w");
    vector<box> boxes;
    for (int k = 0; k < n; k++) {
        box cur;
        for (int i = 0; i < 3; i++) {
            cur[i][0] = bound[i][0] + (bound[i][1] - bound[i][0]) * k / n - size / 2;
            cur[i][1] = bound[i][0] + (bound[i][1] - bound[i][0]) * k / n + size / 2;
        }
        fprintf(fp, "%lf %lf %lf %lf %lf %lf ", ((cur[1][0] + cur[1][1]) / 2 + y_center) / lat_scale,
                ((cur[0][0] + cur[0][1]) / 2 + x_center) / lng_scale,
                ((cur[2][0] + cur[2][1]) / 2 + z_center) / hgt_scale,
                size, size, size
        );
        boxes.push_back(cur);

        double ans1 = query(cur), ans2 = query2(cur), ans3 = query3(cur);
//        cout << ans1 << " " << ans2 << " " << ans3 << endl;

        fprintf(fp, "%.4lf %.4lf %.4lf\n", ans1, ans2, ans3);
    }
    fclose(fp);
}

void database::generateExpandingQueries(string file, int n, double size0) {
    FILE *fp;
    fp = fopen(file.c_str(), "w");
    vector<box> boxes;
    for (int k = 0; k < n; k++) {
        box cur;
        double size = size0 + 5.0 * k / n * size0;
        for (int i = 0; i < 3; i++) {
            cur[i][0] = (bound[i][1] + bound[i][0]) / 2 - size / 2;
            cur[i][1] = (bound[i][1] + bound[i][0]) / 2 + size / 2;
        }
        cur[2][0] = (bound[2][1] + bound[2][0]) / 2 - size0 / 2;
        cur[2][1] = (bound[2][1] + bound[2][0]) / 2 + size0 / 2;

        fprintf(fp, "%lf %lf %lf %lf %lf %lf ", ((cur[1][0] + cur[1][1]) / 2 + y_center) / lat_scale,
                ((cur[0][0] + cur[0][1]) / 2 + x_center) / lng_scale,
                ((cur[2][0] + cur[2][1]) / 2 + z_center) / hgt_scale,
                size, size, size0
        );
        boxes.push_back(cur);

        double ans1 = query(cur), ans2 = query2(cur), ans3 = query3(cur);
//        cout << ans1 << " " << ans2 << " " << ans3 << endl;

        fprintf(fp, "%.4lf %.4lf %.4lf\n", ans1, ans2, ans3);
    }
    fclose(fp);
}

void database::generateRisingQueries(string file, int n, double size) {
    FILE *fp;
    fp = fopen(file.c_str(), "w");
    vector<box> boxes;

    auto x = points[0][0][0];
    auto y = points[0][0][1];
    auto z = points[0][0][2];

    for (int k = 0; k < n; k++) {
        box cur;


        cur[0][0] = x - size / 2;
        cur[0][1] = x + size / 2;

        cur[1][0] = y - size / 2;
        cur[1][1] = y + size / 2;

        cur[2][0] = z - size / 2 - size * k;
        cur[2][1] = z + size / 2 - size * k;




        fprintf(fp, "%lf %lf %lf %lf %lf %lf ", ((cur[1][0] + cur[1][1]) / 2 + y_center) / lat_scale,
                ((cur[0][0] + cur[0][1]) / 2 + x_center) / lng_scale,
                ((cur[2][0] + cur[2][1]) / 2 + z_center) / hgt_scale,
                size, size, size
        );

        boxes.push_back(cur);

        double ans1 = query(cur), ans2 = query2(cur), ans3 = query3(cur);
//        cout << ans1 << " " << ans2 << " " << ans3 << endl;

        fprintf(fp, "%.4lf %.4lf %.4lf\n", ans1, ans2, ans3);
    }
    fclose(fp);
}

void database::generateBoundQueries(string file) {
    FILE *fp;
    fp = fopen(file.c_str(), "w");

    box cur = bound;
    cur[2][0] = 414 - z_center;
    auto sk = 250;
    cur[0][0] += sk;
    cur[1][0] += sk;
    cur[0][1] -= sk;
    cur[1][1] -= sk;


    fprintf(fp, "%lf %lf %lf %lf %lf %lf ", ((cur[1][0] + cur[1][1]) / 2 + y_center) / lat_scale,
            ((cur[0][0] + cur[0][1]) / 2 + x_center) / lng_scale,
            ((cur[2][0] + cur[2][1]) / 2 + z_center) / hgt_scale,
            cur[1][1] - cur[1][0],
            cur[0][1] - cur[0][0],
            cur[2][1] - cur[2][0]
    );



    double ans1 = query(cur), ans2 = query2(cur), ans3 = query3(cur);
//        cout << ans1 << " " << ans2 << " " << ans3 << endl;

    fprintf(fp, "%.4lf %.4lf %.4lf\n", ans1, ans2, ans3);

    fclose(fp);
}



void database::generateQueriesGuaranteeIntersection(string file, int n, double size) {
    FILE *fp;
    fp = fopen(file.c_str(), "w");
    int tot = 0;
    for (;;) {
        double pt[3];
        auto q = bound;
        pt[0] = q.x[0] + (q.x[1] - q.x[0]) * (rand() % M) / M;
        pt[1] = q.y[0] + (q.y[1] - q.y[0]) * (rand() % M) / M;
        pt[2] = q.z[0] + (q.z[1] - q.z[0]) * (rand() % M) / M;
        for (int k = 0; k < 3; k++) {
            q[k][0] = pt[k] - size / 2;
            q[k][1] = pt[k] + size / 2;
        }
        double c = query_alwayssubsetting(q);
        if (c > 1e-5) {
//            cout << c << " " << query2(q) << " " << query3(q) << endl;
            fprintf(fp, "%lf %lf %lf\n", pt[0], pt[1], pt[2]);
            if (++tot == n) break;
        }
    }
    fclose(fp);

}
