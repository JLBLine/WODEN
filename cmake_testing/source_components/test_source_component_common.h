#include <math.h>
#include <unity.h>
#include <stdlib.h>
#include <complex.h>
#include <erfa.h>

#include "constants.h"
#include "woden_struct_defs.h"
#include "woden_precision_defs.h"
#include "FEE_primary_beam.h"

/*******************************************************************************
Here are a bunch of predefined arrays to be used
Some are inputs for tests, others are expected outputs
*******************************************************************************/

user_precision_t azs[] = {4.71238899, 4.71238899, 4.71238899, 4.71238899, 4.71238899,
               4.71238899, 4.71238899, 4.71238899, 4.71238899, 4.71238899,
               4.71238899, 4.71238899, 0.00000000, 0.00000000, 0.00000000,
               1.57079637, 1.57079637, 1.57079637, 1.57079637, 1.57079637,
               1.57079637, 1.57079637, 1.57079637, 1.57079637, 1.57079637,
               1.57079637, 1.57079637};
user_precision_t zas[] = {1.57079631, 1.57079631, 1.57079631, 1.04719766, 1.04719766,
               1.04719766, 0.78539832, 0.78539832, 0.78539832, 0.52359898,
               0.52359898, 0.52359898, 0.00000000, 0.00000000, 0.00000000,
               0.52359875, 0.52359875, 0.52359875, 0.78539820, 0.78539820,
               0.78539820, 1.04719760, 1.04719760, 1.04719760, 1.57079637,
               1.57079637, 1.57079637};

user_precision_t sin_para_angs[] = {1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         1.00000000, 1.00000000, 1.00000000, 1.00000000,
                         0.00000000, 0.00000000, 0.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000, -1.00000000,
                         -1.00000000, -1.00000000, -1.00000000};

user_precision_t cos_para_angs[] = {0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         1.00000000, 1.00000000, 1.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000, 0.00000000,
                         0.00000000, 0.00000000, 0.00000000};

double analy_expec_J00[] = { 0.00000002, 0.52576818, 0.73128019,
  0.88075472, 1.00000000, 0.88075483, 0.73128027, 0.52576824, -0.00000005,
  0.00000002, 0.61822932, 0.81629589, 0.93152102, 1.00000000, 0.93152109,
  0.81629596, 0.61822938, -0.00000006, 0.00000002, 0.52576818, 0.73128019,
  0.88075472, 1.00000000, 0.88075483, 0.73128027, 0.52576824, -0.00000005,
  0.00000002, 0.61822932, 0.81629589, 0.93152102, 1.00000000, 0.93152109,
  0.81629596, 0.61822938, -0.00000006, 0.00000002, 0.52576818, 0.73128019,
  0.88075472, 1.00000000, 0.88075483, 0.73128027, 0.52576824, -0.00000005,
  0.00000002, 0.61822932, 0.81629589, 0.93152102, 1.00000000, 0.93152109,
  0.81629596, 0.61822938, -0.00000006 };

double analy_expec_J11[] = { 0.00000000, 0.26288404, 0.51709310,
  0.76275588, 1.00000000, 0.76275607, 0.51709322, 0.26288410, -0.00000000,
  0.00000000, 0.30911460, 0.57720827, 0.80672078, 1.00000000, 0.80672094,
  0.57720839, 0.30911466, -0.00000000, 0.00000000, 0.26288404, 0.51709310,
  0.76275588, 1.00000000, 0.76275607, 0.51709322, 0.26288410, -0.00000000,
  0.00000000, 0.30911460, 0.57720827, 0.80672078, 1.00000000, 0.80672094,
  0.57720839, 0.30911466, -0.00000000, 0.00000000, 0.26288404, 0.51709310,
  0.76275588, 1.00000000, 0.76275607, 0.51709322, 0.26288410, -0.00000000,
  0.00000000, 0.30911460, 0.57720827, 0.80672078, 1.00000000, 0.80672094,
  0.57720839, 0.30911466, -0.00000000 };

double fee_expec_J00_re[] = { -0.00000375, -0.00007259, 0.00001992, 0.00008827,
  0.00025053, 0.00008904, 0.00002027, -0.00007286, -0.00000362, -0.00000375, -0.00007259,
  0.00001992, 0.00008827, 0.00025053, 0.00008904, 0.00002027, -0.00007286, -0.00000362,
  -0.00000375, -0.00007259, 0.00001992, 0.00008827, 0.00025053, 0.00008904, 0.00002027,
  -0.00007286, -0.00000362, -0.00000375, -0.00007259, 0.00001992, 0.00008827,
  0.00025053, 0.00008904, 0.00002027, -0.00007286, -0.00000362, -0.00000375,
  -0.00007259, 0.00001992, 0.00008827, 0.00025053, 0.00008904, 0.00002027, -0.00007286,
  -0.00000362, -0.00000375, -0.00007259, 0.00001992, 0.00008827, 0.00025053,
  0.00008904, 0.00002027, -0.00007286, -0.00000362 };

double fee_expec_J00_im[] = { 0.00000086, -0.00003368, -0.00001715,
  0.00002948, 0.00012354, 0.00002986, -0.00001680, -0.00003321, 0.00000057,
  0.00000086, -0.00003368, -0.00001715, 0.00002948, 0.00012354, 0.00002986,
  -0.00001680, -0.00003321, 0.00000057, 0.00000086, -0.00003368, -0.00001715,
  0.00002948, 0.00012354, 0.00002986, -0.00001680, -0.00003321, 0.00000057,
  0.00000086, -0.00003368, -0.00001715, 0.00002948, 0.00012354, 0.00002986,
  -0.00001680, -0.00003321, 0.00000057, 0.00000086, -0.00003368, -0.00001715,
  0.00002948, 0.00012354, 0.00002986, -0.00001680, -0.00003321, 0.00000057,
  0.00000086, -0.00003368, -0.00001715, 0.00002948, 0.00012354, 0.00002986,
  -0.00001680, -0.00003321, 0.00000057 };

double fee_expec_J01_re[] = { -0.00299196, 0.05242207, 0.21119083,
  0.09754611, 0.97041277, 0.09754551, 0.21118965, 0.05242191, -0.00299213,
  -0.00299196, 0.05242207, 0.21119083, 0.09754611, 0.97041277, 0.09754551,
  0.21118965, 0.05242191, -0.00299213, -0.00299196, 0.05242207, 0.21119083,
  0.09754611, 0.97041277, 0.09754551, 0.21118965, 0.05242191, -0.00299213,
  -0.00299196, 0.05242207, 0.21119083, 0.09754611, 0.97041277, 0.09754551,
  0.21118965, 0.05242191, -0.00299213, -0.00299196, 0.05242207, 0.21119083,
  0.09754611, 0.97041277, 0.09754551, 0.21118965, 0.05242191, -0.00299213,
  -0.00299196, 0.05242207, 0.21119083, 0.09754611, 0.97041277, 0.09754551,
  0.21118965, 0.05242191, -0.00299213 };

double fee_expec_J01_im[] = { -0.00004542, 0.00969854, 0.03443863,
  0.01937939, 0.24145197, 0.01937985, 0.03443851, 0.00969827, -0.00004550,
  -0.00004542, 0.00969854, 0.03443863, 0.01937939, 0.24145197, 0.01937985,
   0.03443851, 0.00969827, -0.00004550, -0.00004542, 0.00969854, 0.03443863,
   0.01937939, 0.24145197, 0.01937985, 0.03443851, 0.00969827, -0.00004550,
   -0.00004542, 0.00969854, 0.03443863, 0.01937939, 0.24145197, 0.01937985,
   0.03443851, 0.00969827, -0.00004550, -0.00004542, 0.00969854, 0.03443863,
    0.01937939, 0.24145197, 0.01937985, 0.03443851, 0.00969827, -0.00004550,
    -0.00004542, 0.00969854, 0.03443863, 0.01937939, 0.24145197, 0.01937985,
     0.03443851, 0.00969827, -0.00004550 };

double fee_expec_J10_re[] = {0.00169129, -0.02178988, -0.13555926,
  -0.07058298, -0.97021026, -0.07058250, -0.13556084, -0.02179252, 0.00168970,
  0.00169129, -0.02178988, -0.13555926, -0.07058298, -0.97021026, -0.07058250,
  -0.13556084, -0.02179252, 0.00168970, 0.00169129, -0.02178988, -0.13555926,
  -0.07058298, -0.97021026, -0.07058250, -0.13556084, -0.02179252, 0.00168970,
  0.00169129, -0.02178988, -0.13555926, -0.07058298, -0.97021026, -0.07058250,
  -0.13556084, -0.02179252, 0.00168970, 0.00169129, -0.02178988, -0.13555926,
  -0.07058298, -0.97021026, -0.07058250, -0.13556084, -0.02179252, 0.00168970,
  0.00169129, -0.02178988, -0.13555926, -0.07058298, -0.97021026, -0.07058250,
  -0.13556084, -0.02179252, 0.00168970  };

double fee_expec_J10_im[] = { 0.00086604, 0.01911517, -0.01697326,
  -0.03916865, -0.24226442, -0.03917003, -0.01697404, 0.01911435, 0.00086739,
  0.00086604, 0.01911517, -0.01697326, -0.03916865, -0.24226442, -0.03917003,
  -0.01697404, 0.01911435, 0.00086739, 0.00086604, 0.01911517, -0.01697326,
  -0.03916865, -0.24226442, -0.03917003, -0.01697404, 0.01911435, 0.00086739,
  0.00086604, 0.01911517, -0.01697326, -0.03916865, -0.24226442, -0.03917003,
  -0.01697404, 0.01911435, 0.00086739, 0.00086604, 0.01911517, -0.01697326,
  -0.03916865, -0.24226442, -0.03917003, -0.01697404, 0.01911435, 0.00086739,
   0.00086604, 0.01911517, -0.01697326, -0.03916865, -0.24226442, -0.03917003,
   -0.01697404, 0.01911435, 0.00086739 };

double fee_expec_J11_re[] = { -0.00000147, -0.00000903, -0.00001642,
  -0.00006570, -0.00024691, -0.00006504, -0.00001520, -0.00000766, -0.00000182,
  -0.00000147, -0.00000903, -0.00001642, -0.00006570, -0.00024691, -0.00006504,
  -0.00001520, -0.00000766, -0.00000182, -0.00000147, -0.00000903, -0.00001642,
  -0.00006570, -0.00024691, -0.00006504, -0.00001520, -0.00000766, -0.00000182,
  -0.00000147, -0.00000903, -0.00001642, -0.00006570, -0.00024691, -0.00006504,
  -0.00001520, -0.00000766, -0.00000182, -0.00000147, -0.00000903, -0.00001642,
  -0.00006570, -0.00024691, -0.00006504, -0.00001520, -0.00000766, -0.00000182,
  -0.00000147, -0.00000903, -0.00001642, -0.00006570, -0.00024691, -0.00006504,
  -0.00001520, -0.00000766, -0.00000182 };

double fee_expec_J11_im[] = { 0.00000424, 0.00007635, 0.00002101,
  -0.00010817, -0.00011722, -0.00010677, 0.00002321, 0.00007781, 0.00000384,
  0.00000424, 0.00007635, 0.00002101, -0.00010817, -0.00011722, -0.00010677,
   0.00002321, 0.00007781, 0.00000384, 0.00000424, 0.00007635, 0.00002101,
   -0.00010817, -0.00011722, -0.00010677, 0.00002321, 0.00007781, 0.00000384,
    0.00000424, 0.00007635, 0.00002101, -0.00010817, -0.00011722, -0.00010677,
    0.00002321, 0.00007781, 0.00000384, 0.00000424, 0.00007635, 0.00002101,
    -0.00010817, -0.00011722, -0.00010677, 0.00002321, 0.00007781, 0.00000384,
    0.00000424, 0.00007635, 0.00002101, -0.00010817, -0.00011722, -0.00010677,
    0.00002321, 0.00007781, 0.00000384 };



double fee_expec_interp_J00_re[] = { -0.00000190, -0.00007459, -0.00002645,
  0.00007861, 0.00024664, 0.00007857, -0.00002519, -0.00007411, -0.00000322,
  0.00000392, 0.00002112, -0.00003946, -0.00000576, 0.00016112, -0.00000579,
  -0.00003943, 0.00002129, 0.00000283, -0.00000190, -0.00007459, -0.00002645,
  0.00007861, 0.00024664, 0.00007857, -0.00002519, -0.00007411, -0.00000322,
  0.00000392, 0.00002112, -0.00003946, -0.00000576, 0.00016112, -0.00000579,
  -0.00003943, 0.00002129, 0.00000283, -0.00000190, -0.00007459, -0.00002645,
  0.00007861, 0.00024664, 0.00007857, -0.00002519, -0.00007411, -0.00000322,
  0.00000392, 0.00002112, -0.00003946, -0.00000576, 0.00016112, -0.00000579,
  -0.00003943, 0.00002129, 0.00000283 };

double fee_expec_interp_J00_im[] = { 0.00000691, 0.00001597, -0.00001522,
  -0.00001145, 0.00001658, -0.00001049, -0.00001483, 0.00001387, 0.00000682,
  -0.00000116, 0.00000937, 0.00002198, -0.00001854, -0.00017165, -0.00001998,
  0.00002280, 0.00001084, -0.00000343, 0.00000691, 0.00001597, -0.00001522,
  -0.00001145, 0.00001658, -0.00001049, -0.00001483, 0.00001387, 0.00000682,
  -0.00000116, 0.00000937, 0.00002198, -0.00001854, -0.00017165, -0.00001998,
  0.00002280, 0.00001084, -0.00000343, 0.00000691, 0.00001597, -0.00001522,
  -0.00001145, 0.00001658, -0.00001049, -0.00001483, 0.00001387, 0.00000682,
  -0.00000116, 0.00000937, 0.00002198, -0.00001854, -0.00017165, -0.00001998,
  0.00002280, 0.00001084, -0.00000343  };

double fee_expec_interp_J01_re[] = { -0.00551808, -0.05004035, 0.15862949,
  0.19153106, 0.94742141, 0.19153010, 0.15863058, -0.05003839, -0.00551712,
  -0.00091142, -0.13863655, -0.02048749, 0.20174565, 0.99971698, 0.20174408,
  -0.02048994, -0.13863725, -0.00091196, -0.00551808, -0.05004035, 0.15862949,
  0.19153106, 0.94742141, 0.19153010, 0.15863058, -0.05003839, -0.00551712,
  -0.00091142, -0.13863655, -0.02048749, 0.20174565, 0.99971698, 0.20174408,
  -0.02048994, -0.13863725, -0.00091196, -0.00551808, -0.05004035, 0.15862949,
  0.19153106, 0.94742141, 0.19153010, 0.15863058, -0.05003839, -0.00551712,
  -0.00091142, -0.13863655, -0.02048749, 0.20174565, 0.99971698, 0.20174408,
  -0.02048994, -0.13863725, -0.00091196 };

double fee_expec_interp_J01_im[] = { -0.00062418, -0.00817662, 0.02841667,
  0.03704812, 0.31998856, 0.03704975, 0.02841827, -0.00817651, -0.00062417,
  0.00066865, 0.03425132, -0.00134228, -0.04199880, -0.02378985, -0.04200162,
  -0.00134435, 0.03425103, 0.00066920, -0.00062418, -0.00817662, 0.02841667,
  0.03704812, 0.31998856, 0.03704975, 0.02841827, -0.00817651, -0.00062417,
  0.00066865, 0.03425132, -0.00134228, -0.04199880, -0.02378985, -0.04200162,
  -0.00134435, 0.03425103, 0.00066920, -0.00062418, -0.00817662, 0.02841667,
  0.03704812, 0.31998856, 0.03704975, 0.02841827, -0.00817651, -0.00062417,
  0.00066865, 0.03425132, -0.00134228, -0.04199880, -0.02378985, -0.04200162,
  -0.00134435, 0.03425103, 0.00066920 };

double fee_expec_interp_J10_re[] = { 0.00312704, 0.03543797, -0.09541299,
  -0.15001436, -0.94700971, -0.15001204, -0.09541293, 0.03543539, 0.00312805,
  -0.00008473, 0.07794426, 0.02281578, -0.21326641, -0.99976741, -0.21326888,
  0.02281533, 0.07794430, -0.00008479, 0.00312704, 0.03543797, -0.09541299,
  -0.15001436, -0.94700971, -0.15001204, -0.09541293, 0.03543539, 0.00312805,
  -0.00008473, 0.07794426, 0.02281578, -0.21326641, -0.99976741, -0.21326888,
  0.02281533, 0.07794430, -0.00008479, 0.00312704, 0.03543797, -0.09541299,
  -0.15001436, -0.94700971, -0.15001204, -0.09541293, 0.03543539, 0.00312805,
  -0.00008473, 0.07794426, 0.02281578, -0.21326641, -0.99976741, -0.21326888,
  0.02281533, 0.07794430, -0.00008479 };

double fee_expec_interp_J10_im[] = { -0.00000324, 0.01859633, -0.01843399,
  -0.06486252, -0.32120493, -0.06486467, -0.01843614, 0.01859592, -0.00000058,
  -0.00260255, -0.02655038, 0.01283152, 0.01974520, 0.02156680, 0.01974410,
  0.01283253, -0.02655104, -0.00260283, -0.00000324, 0.01859633, -0.01843399,
  -0.06486252, -0.32120493, -0.06486467, -0.01843614, 0.01859592, -0.00000058,
  -0.00260255, -0.02655038, 0.01283152, 0.01974520, 0.02156680, 0.01974410,
  0.01283253, -0.02655104, -0.00260283, -0.00000324, 0.01859633, -0.01843399,
  -0.06486252, -0.32120493, -0.06486467, -0.01843614, 0.01859592, -0.00000058,
  -0.00260255, -0.02655038, 0.01283152, 0.01974520, 0.02156680, 0.01974410,
  0.01283253, -0.02655104, -0.00260283 };

double fee_expec_interp_J11_re[] = { -0.00000362, 0.00000946, 0.00002380,
  -0.00005941, -0.00024983, -0.00005855, 0.00002271, 0.00000699, -0.00000172,
  0.00000063, -0.00000679, 0.00003557, 0.00000351, -0.00016745, 0.00000231,
  0.00003550, -0.00000641, 0.00000197, -0.00000362, 0.00000946, 0.00002380,
  -0.00005941, -0.00024983, -0.00005855, 0.00002271, 0.00000699, -0.00000172,
  0.00000063, -0.00000679, 0.00003557, 0.00000351, -0.00016745, 0.00000231,
  0.00003550, -0.00000641, 0.00000197, -0.00000362, 0.00000946, 0.00002380,
  -0.00005941, -0.00024983, -0.00005855, 0.00002271, 0.00000699, -0.00000172,
  0.00000063, -0.00000679, 0.00003557, 0.00000351, -0.00016745, 0.00000231,
  0.00003550, -0.00000641, 0.00000197 };

double fee_expec_interp_J11_im[] = { 0.00000210, 0.00005938, 0.00003546,
  -0.00004659, -0.00001182, -0.00004497, 0.00003518, 0.00005756, 0.00000348,
  -0.00000258, 0.00002352, 0.00004401, -0.00000665, 0.00016438, -0.00000648,
  0.00004322, 0.00002389, -0.00000270, 0.00000210, 0.00005938, 0.00003546,
  -0.00004659, -0.00001182, -0.00004497, 0.00003518, 0.00005756, 0.00000348,
  -0.00000258, 0.00002352, 0.00004401, -0.00000665, 0.00016438, -0.00000648,
  0.00004322, 0.00002389, -0.00000270, 0.00000210, 0.00005938, 0.00003546,
  -0.00004659, -0.00001182, -0.00004497, 0.00003518, 0.00005756, 0.00000348,
  -0.00000258, 0.00002352, 0.00004401, -0.00000665, 0.00016438, -0.00000648,
  0.00004322, 0.00002389, -0.00000270 };

double MWA_analy_expec_J00_re[] = { -0.000000004365, -0.105094352084,
  -0.026591963099, 0.267847379368, 0.893345358913, 0.267847728174,
  -0.026591881348, -0.105094352255, 0.000000011230, 0.000000001347,
  0.150284040279, 0.042977910381, -0.226477767063, 0.893345358913,
  -0.226477788010, 0.042977774708, 0.150284054600, -0.000000003466,
  -0.000000004365, -0.105094352084, -0.026591963099, 0.267847379368,
  0.893345358913, 0.267847728174, -0.026591881348, -0.105094352255,
  0.000000011230, 0.000000001347, 0.150284040279, 0.042977910381,
  -0.226477767063, 0.893345358913, -0.226477788010, 0.042977774708,
  0.150284054600, -0.000000003466, -0.000000004365, -0.105094352084,
  -0.026591963099, 0.267847379368, 0.893345358913, 0.267847728174,
  -0.026591881348, -0.105094352255, 0.000000011230, 0.000000001347,
  0.150284040279, 0.042977910381, -0.226477767063, 0.893345358913,
  -0.226477788010, 0.042977774708, 0.150284054600, -0.000000003466 };

double MWA_analy_expec_J11_re[] = { -0.000000004886, -0.058820673906,
  -0.021048251099, 0.259656170546, 1.000000000000, 0.259656508684,
  -0.021048186391, -0.058820674002, 0.000000000000, 0.000000001508,
  0.084113069363, 0.034018167295, -0.219551708320, 1.000000000000,
  -0.219551728627, 0.034018059906, 0.084113077379, -0.000000000000,
  -0.000000004886, -0.058820673906, -0.021048251099, 0.259656170546,
  1.000000000000, 0.259656508684, -0.021048186391, -0.058820674002,
  0.000000000000, 0.000000001508, 0.084113069363, 0.034018167295,
  -0.219551708320, 1.000000000000, -0.219551728627, 0.034018059906,
  0.084113077379, -0.000000000000, 0.000000000000, -0.058820673906,
  -0.021048251099, 0.259656170546, 1.000000000000, 0.259656508684,
  -0.021048186391, -0.058820674002, 0.000000000000, -0.000000000000,
  0.084113069363, 0.034018167295, -0.219551708320, 1.000000000000,
  -0.219551728627, 0.034018059906, 0.084113077379, -0.000000000000 };

double MWA_analy_expec_J01_re[] = { -0.000000000000, 0.045782069703,
  0.009458468353, -0.067366319296, 0.000000000000, 0.067366407025,
  -0.009458439275, -0.045782069777, 0.000000005649, 0.000000000000,
  -0.065467974928, -0.015286769304, 0.056961444258, 0.000000000000,
  -0.056961449526, 0.015286721047, 0.065467981167, -0.000000001744,
  -0.000000000000, 0.045782069703, 0.009458468353, -0.067366319296,
  0.000000000000, 0.067366407025, -0.009458439275, -0.045782069777,
  0.000000005649, 0.000000000000, -0.065467974928, -0.015286769304,
  0.056961444258, 0.000000000000, -0.056961449526, 0.015286721047,
  0.065467981167, -0.000000001744, 0.000000002196, 0.045782069703,
  0.009458468353, -0.067366319296, 0.000000000000, 0.067366407025,
  -0.009458439275, -0.045782069777, 0.000000005649, -0.000000000678,
  -0.065467974928, -0.015286769304, 0.056961444258, 0.000000000000,
  -0.056961449526, 0.015286721047, 0.065467981167, -0.000000001744};