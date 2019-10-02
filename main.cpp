#include "database.h"
#include <iomanip>


// Given configurations, this function can automatically generate result experiment files as the paper needs

void generateFile(string infoFile) {
    ifstream info(infoFile);

    string _, title, dataset, alg;
    int angle, cell;
    vector<string> algs;
    vector<int> angles, cells;

    info >> _ >> title;
    assert(_ == "title:");
    info >> _ >> dataset;
    assert(_ == "dataset:");
    info >> _ >> alg;
    assert(_ == "alg:");
    while (alg.length() == 2) {
        algs.push_back(alg.substr(0, 1));
        info >> alg;
    }
    algs.push_back(alg);

    info >> _;
    assert(_ == "angles:");

    for (int i = 0; i < algs.size(); i++) {
        info >> angle;
        if (i != algs.size() - 1) info >> _;
        angles.push_back(angle);
    }

    info >> _;
    assert(_ == "cells:");

    for (int i = 0; i < algs.size(); i++) {
        info >> cell;
        if (i != algs.size() - 1) info >> _;
        cells.push_back(cell);
    }

    string s;
    vector<double> Rs, alphas, cubeSizes;

    double maxR = 0;

    info >> _ >> s;
    assert(_ == "R:");
    while (s[s.length() - 1] == ',') {
        double R = atof(s.substr(0, s.length() - 1).c_str());
        maxR = max(maxR, R);
        Rs.push_back(R);
        info >> s;
    }
    maxR = max(maxR, atof(s.c_str()));
    Rs.push_back(atof(s.c_str()));

    info >> _ >> s;
    assert(_ == "alpha:");
    while (s[s.length() - 1] == ',') {
        alphas.push_back(atof(s.substr(0, s.length() - 1).c_str()));
        info >> s;
    }
    alphas.push_back(atof(s.c_str()));


    double maxCubeSize = 0;
    info >> _ >> s;
    assert(_ == "cubeSize:");
    while (s[s.length() - 1] == ',') {
        double cubeSize = atof(s.substr(0, s.length() - 1).c_str());
        maxCubeSize = max(maxCubeSize, cubeSize);
        cubeSizes.push_back(cubeSize);
        info >> s;
    }
    maxCubeSize = max(maxCubeSize, atof(s.c_str()));
    cubeSizes.push_back(atof(s.c_str()));


    vector<string> names;

    info >> _ >> s;
    assert(_ == "names:");
    while (s[s.length() - 1] == ',') {
        names.push_back(s.substr(0, s.length() - 1).c_str());
        info >> s;
    }
    names.push_back(s.substr(0, s.length()).c_str());

    int num;
    info >> _ >> num;
    assert(_ == "num:");


    info.close();

    database db(dataset, maxR, alphas[0]);
    db.generateQueriesGuaranteeIntersection("queries_" + title + ".txt", num, maxCubeSize);

    vector<double> ansList[100];
    vector<double> timeList[100];

    int cur = 0;
    for (int iAlg = 0; iAlg < algs.size(); iAlg++)
        for (auto R: Rs)
            for (auto alpha: alphas)
                for (auto cubeSize: cubeSizes) {
                    db = database(dataset, R, alpha);

                    ifstream queryFile("queries_" + title + ".txt");

                    for (int it = 0; it < num; it++) {
                        cout << "it: " << it << endl;
                        double x, y, z;
                        queryFile >> x >> y >> z;
//                        cout << x << " " << y << " " << z << endl;
                        box q;
                        q.x[0] = x - cubeSize / 2;
                        q.x[1] = x + cubeSize / 2;
                        q.y[0] = y - cubeSize / 2;
                        q.y[1] = y + cubeSize / 2;
                        q.z[0] = z - cubeSize / 2;
                        q.z[1] = z + cubeSize / 2;
                        auto start = high_resolution_clock::now();
                        double ans;
                        if (algs[iAlg] == "1") {
                            ans = db.query(q);
                        }
                        else if (algs[iAlg] == "2") {
                            ans = db.query2(q, angles[iAlg]);
                        }
                        else if (algs[iAlg] == "3") {
                            ans = db.query3(q, angles[iAlg], cells[iAlg]);
                        }
                        int duration = duration_cast<microseconds>(high_resolution_clock::now() - start).count();
                        ansList[it].push_back(ans);
                        timeList[it].push_back(duration);
                    }
                    queryFile.close();

                }
    ofstream out(title + ".csv");

    for (string name: names)
        out << name + "_result" << ", ";
    for (int i = 0; i < names.size(); i++)
    {
        string name = names[i];
        out << name + "_time";
        if (i != names.size() - 1) out << ", ";
    }
    out << endl;
    for (int i = 0; i < num; i++) {
        for (double ans: ansList[i])
            out << fixed << setprecision(4) << ans << ", ";
        for (int j = 0; j < names.size(); j++)
        {
            out << timeList[i][j];
            if (j != names.size() - 1) out << ", ";
        }
        out << endl;

    }
    for (int j = 0; j < names.size(); j++) {
        double sum = 0;
        for (int i = 0; i < num; i++) sum += ansList[i][j];
        sum /= num;
        cout << sum * 100 << endl;
    }
    out.close();
}



int main(int argc, char **argv) {

    // initialize a database with a file path and the parameters
//    database db("dataset_large.csv", 500, 45, 40000, 42000, true);
//
//
//    auto start = high_resolution_clock::now();
//    cout << db.query(db.bound) << endl;
//    cout << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << endl;
//    start = high_resolution_clock::now();
//    cout << db.query2(db.bound) << endl;
//    cout << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << endl;
//    start = high_resolution_clock::now();
//    cout << db.query3(db.bound) << endl;
//    cout << duration_cast<microseconds>(high_resolution_clock::now() - start).count() << endl;
//


//    db.generateContinuousQueries("demoQuery.txt", 30, 100);
//    db.generateExpandingQueries("demoQuery2.txt", 30, 100);
//    db.generateBoundQueries("demoQuery3.txt");
    string s = argv[1];
    generateFile(s);
}