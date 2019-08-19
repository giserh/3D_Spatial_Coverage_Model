#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <emscripten/emscripten.h>

// 一旦WASM模块被加载，main()中的代码就会执行
int main(int argc, char ** argv) {
    printf("WebAssembly module loaded\n");
}

struct box {
    double x[2], y[2], z[2];
    double *operator [] (int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
};

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



double VCM(double *lat, double *lng, double hgt*, double *yaw, double *pitch, double *roll, int cnt, double R, double alpha, int camera) {

//    bool extrinsic = true;


    points53 *points = (struct points53*) malloc(cnt * sizeof(struct points53));
    box *bds = (struct box*) malloc(cnt * sizeof(struct box));

    double lat_min, lat_max, lat_scale;
    double lng_min, lng_max, lng_scale;
    double hgt_min, hgt_max, hgt_scale;
    double x_center, y_center, z_center;
    double Rv, alpha;
    int cnt;
    box bound;


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

//    cout << "mean latitude: " << (lat_min + lat_max) / 2 << endl;
//    cout << "mean longitude: " << (lng_min + lng_max) / 2 << endl;

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
        for (int k = 0; k < 3; k++) points[i][0][k] = 0;
        double phi = yaw[i], theta = pitch[i], psi = roll[i];
//        if (!extrinsic) {
//            phi = -phi;
//            theta = -theta;
//            psi = -psi;
//        }
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
//        if (!extrinsic) {
//            for (int j = 0; j < 3; j++)
//                for (int k = j; k < 3; k++)
//                    swap(Rot[j][k], Rot[k][j]);
//
//        }
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

    return 0;


}

