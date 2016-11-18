#pragma once
#include <cmath>

#define COEFF_LOG2E   1.442695040888963
#define COEFF_A        4503599627370496
#define COEFF_B     4607182418800017408

#define COEFF_P2_A   -0.002473940320180
#define COEFF_P2_B    0.348948919055000
#define COEFF_P2_C   -0.344001038415000

#define COEFF_P3_A    0.000106959669219
#define COEFF_P3_B    0.303543209079000
#define COEFF_P3_C   -0.224339000371000
#define COEFF_P3_D   -0.079204208707200

#define COEFF_P4_A   -3.70239070165e-06
#define COEFF_P4_B    0.307033839653000
#define COEFF_P4_C   -0.241638340487000
#define COEFF_P4_D   -0.051690438287400
#define COEFF_P4_E   -0.013697656096900

#define COEFF_P5_A    1.06838616541e-07
#define COEFF_P5_B    0.306845249040000
#define COEFF_P5_C   -0.240139718656000
#define COEFF_P5_D   -0.055866234366400
#define COEFF_P5_E   -0.008942834565720
#define COEFF_P5_F   -0.001896461451990

//#define COEFF_P6_A   -2.64320006772e-09
//#define COEFF_P6_B    0.306853075389000
//#define COEFF_P6_C   -0.240230549817000
//#define COEFF_P6_D   -0.055480222052700
//#define COEFF_P6_E   -0.009684975075760
//#define COEFF_P6_F   -0.001238430894380
//#define COEFF_P6_G   -0.000218892263012

#define COEFF_P6_A   -2.64303273610414963822e-09
#define COEFF_P6_B    3.06853075372807815313e-01
#define COEFF_P6_C   -2.40230549677691723742e-01
#define COEFF_P6_D   -5.54802224547989303316e-02
#define COEFF_P6_E   -9.68497459444197204836e-03
#define COEFF_P6_F   -1.23843111224273085859e-03
#define COEFF_P6_G   -2.18892247566917477666e-04

#define COEFF_P7_A    5.72283550146e-11
#define COEFF_P7_B    0.306852812183000
#define COEFF_P7_C   -0.240226356055000
#define COEFF_P7_D   -0.055505302286300
#define COEFF_P7_E   -0.009613506234810
#define COEFF_P7_F   -0.001343024385490
#define COEFF_P7_G   -0.000142962480244
#define COEFF_P7_H   -2.16607410310e-05

#define COEFF_P8_A   -1.10148979590e-12
#define COEFF_P8_B    0.306852819617000
#define COEFF_P8_C   -0.240226511645000
#define COEFF_P8_D   -0.055504060970800
#define COEFF_P8_E   -0.009618371833940
#define COEFF_P8_F   -0.001332664048900
#define COEFF_P8_G   -0.000155186861569
#define COEFF_P8_H   -1.41484303383e-05
#define COEFF_P8_I   -1.87582398427e-06

#define COEFF_P9_A    1.90763813177e-14
#define COEFF_P9_B    0.306852819436000
#define COEFF_P9_C   -0.240226506835000
#define COEFF_P9_D   -0.055504110246200
#define COEFF_P9_E   -0.009618118912090
#define COEFF_P9_F   -0.001333393498320
#define COEFF_P9_G   -0.000153950778902
#define COEFF_P9_H   -1.53694079792e-05
#define COEFF_P9_I   -1.22534748028e-06
#define COEFF_P9_J   -1.44410574187e-07

#define COEFF_P10_A   7.32388148129e-13
#define COEFF_P10_B   3.06852819216e-01
#define COEFF_P10_C  -2.40226499275e-01
#define COEFF_P10_D  -5.55042073859e-02
#define COEFF_P10_E  -9.61749102796e-03
#define COEFF_P10_F  -1.33571753728e-03
#define COEFF_P10_G  -1.48718480159e-04
#define COEFF_P10_H  -2.26598047213e-05
#define COEFF_P10_I   4.91492761180e-06
#define COEFF_P10_J  -3.00875847392e-06
#define COEFF_P10_K   5.68126156224e-07

double expFast2(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P2_A + xf *
         (COEFF_P2_B + xf *
          COEFF_P2_C );
    long int castInt = (long int) COEFF_A * x + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast3(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P3_A + xf *
         (COEFF_P3_B + xf *
         (COEFF_P3_C + xf *
          COEFF_P3_D ));
    long int castInt = (long int) COEFF_A * x + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast4(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P4_A + xf *
         (COEFF_P4_B + xf *
         (COEFF_P4_C + xf *
         (COEFF_P4_D + xf *
          COEFF_P4_E )));
    long int castInt = (long int) COEFF_A * x + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast5(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P5_A + xf *
         (COEFF_P5_B + xf *
         (COEFF_P5_C + xf *
         (COEFF_P5_D + xf *
         (COEFF_P5_E + xf *
          COEFF_P5_F ))));
    long int castInt = (long int) COEFF_A * x + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast6(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P6_A + xf *
         (COEFF_P6_B + xf *
         (COEFF_P6_C + xf *
         (COEFF_P6_D + xf *
         (COEFF_P6_E + xf *
         (COEFF_P6_F + xf *
          COEFF_P6_G )))));
    long int castInt = (long int) COEFF_A * x + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast7(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P7_A + xf *
         (COEFF_P7_B + xf *
         (COEFF_P7_C + xf *
         (COEFF_P7_D + xf *
         (COEFF_P7_E + xf *
         (COEFF_P7_F + xf *
         (COEFF_P7_G + xf *
          COEFF_P7_H ))))));
    long int castInt = (long int) COEFF_A * x + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast8(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P8_A + xf *
         (COEFF_P8_B + xf *
         (COEFF_P8_C + xf *
         (COEFF_P8_D + xf *
         (COEFF_P8_E + xf *
         (COEFF_P8_F + xf *
         (COEFF_P8_G + xf *
         (COEFF_P8_H + xf *
          COEFF_P8_I )))))));
    long int castInt = (long int) COEFF_A * x + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast9(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P9_A + xf *
         (COEFF_P9_B + xf *
         (COEFF_P9_C + xf *
         (COEFF_P9_D + xf *
         (COEFF_P9_E + xf *
         (COEFF_P9_F + xf *
         (COEFF_P9_G + xf *
         (COEFF_P9_H + xf *
         (COEFF_P9_I + xf *
          COEFF_P9_J ))))))));
    long int castInt = (long int) COEFF_A * static_cast<long double>(x) + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}

double expFast10(double x) {
    x *= COEFF_LOG2E;
    const double xf = x - floor(x);
    x -=  COEFF_P10_A + xf *
         (COEFF_P10_B + xf *
         (COEFF_P10_C + xf *
         (COEFF_P10_D + xf *
         (COEFF_P10_E + xf *
         (COEFF_P10_F + xf *
         (COEFF_P10_G + xf *
         (COEFF_P10_H + xf *
         (COEFF_P10_I + xf *
         (COEFF_P10_J + xf *
          COEFF_P10_K )))))))));
    long int castInt = (long int) COEFF_A * static_cast<long double>(x) + COEFF_B;
    return reinterpret_cast<double&>(castInt);
}
