//
// Created by 马泽余 on 2019-07-26.
//

#ifndef COVERAGE_DATABASE_H
#define COVERAGE_DATABASE_H
//#define N 1000000
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <thread>
#include <vector>
#include <chrono>




using namespace std;
using namespace std::chrono;

// bounding box

struct box {
    double x[2], y[2], z[2];
    double *operator [] (int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
};

// points struct
struct point3 {
    double a[3];
    double &operator [] (int i) {
        return a[i];
    }
};

struct points53 {
    point3 a[5];
    point3 &operator [] (int i) {
        return a[i];
    }
};


// KD-tree node
struct node {
    node *lc, *rc;
    box b;
    int l, r;
    int dim;
};


// class database, which host all data and perform calculations
class database {


public:

    vector<double> lat, lng, hgt, yaw, pitch, roll, timestamp;
    vector<points53> points;

    vector<point3> pointOne;
    vector<int> ind;
    node *root;

    double lat_min, lat_max, lat_scale;
    double lng_min, lng_max, lng_scale;
    double hgt_min, hgt_max, hgt_scale;
    double x_center, y_center, z_center;
    double Rv, alpha;
    int cnt;
    box bound;
    vector<box> bds;

    database(string dataFile, double R=50, double alpha=45, int st=0, int ed=-1, bool extrinsic=true);
    double query(box q, int threadNum=10);
    double query_alwayssubsetting(box q, int threadNum=10);

    double query2(box q, int angles=8, int threadNum=10);
    double query3(box q, int angles=8, int cells=8, int threadNum=10);


    void generateQueries(string file, int n);
    void generateQueriesGuaranteeIntersection(string file, int n, double size);

    void generateContinuousQueries(string file, int n, double size);
    void generateExpandingQueries(string file, int n, double size);
    void generateRisingQueries(string file, int n = 3, double size = 100);
    void generateBoundQueries(string file);


};


#endif //COVERAGE_DATABASE_H
