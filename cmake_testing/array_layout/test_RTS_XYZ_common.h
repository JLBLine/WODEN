#define LST_BEFORE 0.0077340
#define LST_AFTER 0.0042459
#define JD_DATE 2457278.2010995


// Some expected X,Y,Z values when calculating from e,n,h in
// example_array_layout.txt, with no precession
float expec_X_noprec[] = {53.00956, 93.45293 , 173.87949, 106.47388,
                    258.35040, 270.51578, 31.45595, 242.61703};

float expec_Y_noprec[] = {84.00000, 12.00000, 86.00000, 780.00000,
                    813.00000, 899.00000, 460.00000, 810.00000};

float expec_Z_noprec[] = {98.70657, 179.10765, 334.54434, 200.54254,
                    498.02115, 535.55786, 62.53418, 464.51801};

// Some expected X,Y,Z values when calculating from e,n,h in
// example_array_layout.txt, with precession
float expec_X_prec[] = {53.16233, 93.72652, 174.39215, 106.80007,
                  259.13132, 271.35620, 31.56328, 243.34676};

float expec_Y_prec[] = {83.99356, 11.98839, 85.97833, 779.98694,
                  812.96777, 898.96552, 459.99597, 809.96991};

float expec_Z_prec[] = {98.62984, 178.96541, 334.28296, 200.41980,
                  497.66794, 535.19043, 62.50969, 464.18869};


// Expected baseline lengths when using expec_X_noprec, expec_Y_noprec,
// expec_Z_noprec
float expec_X_diffs_noprec[] = { -40.44337, -120.86994, -53.46433, -205.34085,
                                -217.50626, 21.55361, -189.60748, -80.42657,
                                -13.02096, -164.89746, -177.06287, 61.99698,
                                -149.16412, 67.40561, -84.47090, -96.63631,
                                142.42355, -68.73755, -151.87651, -164.04192,
                                75.01794, -136.14316, -12.16541, 226.89445,
                                15.73335, 239.05986, 27.89876, -211.16110 };

float expec_Y_diffs_noprec[] = { 72.00000, -2.00000, -696.00000, -729.00000,
                                 -815.00000, -376.00000, -726.00000, -74.00000,
                                 -768.00000, -801.00000, -887.00000, -448.00000,
                                 -798.00000, -694.00000, -727.00000, -813.00000,
                                 -374.00000, -724.00000, -33.00000, -119.00000,
                                 320.00000, -30.00000, -86.00000, 353.00000,
                                 3.00000, 439.00000, 89.00000, -350.00000 };

float expec_Z_diffs_noprec[] = { -80.40108, -235.83777, -101.83598, -399.31458,
                                 -436.85129, 36.17239, -365.81143, -155.43669,
                                 -21.43489, -318.91351, -356.45020, 116.57347,
                                 -285.41034, 134.00180, -163.47681, -201.01352,
                                 272.01016, -129.97366, -297.47861, -335.01532,
                                 138.00836, -263.97546, -37.53671, 435.48697,
                                 33.50314, 473.02368, 71.03986, -401.98383 };

// Expected baseline lengths when using expec_X_prec, expec_Y_prec,
// expec_Z_prec
float expec_X_diffs_prec[] = { -40.56419, -121.22983, -53.63774, -205.96898,
                               -218.19389, 21.59905, -190.18443, -80.66564,
                               -13.07355, -165.40479, -177.62970, 62.16324,
                               -149.62024, 67.59209, -84.73915, -96.96407,
                               142.82889, -68.95461, -152.33124, -164.55615,
                               75.23679, -136.54669, -12.22491, 227.56804,
                               15.78455, 239.79295, 28.00946, -211.78349 };

float expec_Y_diffs_prec[] = { 72.00517, -1.98477, -695.99341, -728.97424,
                               -814.97192, -376.00241, -725.97632, -73.98994,
                               -767.99854, -800.97937, -886.97711, -448.00757,
                               -797.98151, -694.00861, -726.98944, -812.98718,
                               -374.01764, -723.99158, -32.98083, -118.97858,
                               319.99097, -29.98297, -85.99774, 352.97180,
                               2.99786, 438.96954, 88.99561, -349.97394 };

float expec_Z_diffs_prec[] = { -80.33556, -235.65311, -101.78996, -399.03809,
                               -436.56058, 36.12016, -365.55884, -155.31755,
                               -21.45439, -318.70251, -356.22504, 116.45572,
                               -285.22327, 133.86316, -163.38498, -200.90747,
                               271.77325, -129.90573, -297.24814, -334.77063,
                               137.91011, -263.76889, -37.52249, 435.15826,
                               33.47925, 472.68073, 71.00174, -401.67902 };