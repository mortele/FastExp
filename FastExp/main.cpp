#include <iostream>
#include <iomanip>
#include <fstream>
#include "fastexp.h"

using namespace std;


int main(int, char**) {
    double      dx   =  0.01;
    //double      xMin =  -709;
    double      xMin =  -30;
    double      xMax =   0;
    long int    N    = (xMax - xMin) / dx;

    double values  [11];
    double maxError[11];
    double xvalues [11];
    for (int i=0; i<11; i++) values  [i] = 0;
    for (int i=0; i<11; i++) maxError[i] = 0;

    ofstream outFile;
    outFile.open("../FastExp/data.dat", ios::out);

    for (long int i = 0; i < N; i++) {
        const double x   = xMin + i * dx;
        const double EXP = exp(x);

        values[2]  = expFast2 (x);
        values[3]  = expFast3 (x);
        values[4]  = expFast4 (x);
        values[5]  = expFast5 (x);
        values[6]  = expFast6 (x);
        values[7]  = expFast7 (x);
        values[8]  = expFast8 (x);
        values[9]  = expFast9 (x);
        values[10] = expFast10(x);

        outFile << setprecision(20) << x << " ";
        for (int i = 2; i < 11; i++) {
            const double diff = fabs(EXP-values[i]);
            outFile << setprecision(20) << diff << " ";
            if (diff > maxError[i]) {
                maxError[i] = diff;
                xvalues [i] = x;
            }
        }
        outFile << EXP << endl;
    }
    outFile.close();

    cout << "MAX ERRORS:" << endl;
    for (int i = 2; i < 11; i++) {
        cout << i << " : " << setprecision(20) << maxError[i] << "    x: " << xvalues[i] << endl;
    }


    return 0;
}
