#pragma once
#include "constants.h"
#include "woden_precision_defs.h"
#include "primary_beam_cpu.h"

void test_GaussBeam_GivesCorrectlValues(int do_gpu);

void test_GaussBeam_GivesCorrectmValues(int do_gpu);

void test_GaussBeam_GivesCorrectlValuesByFreq(int do_gpu);

void test_GaussBeam_GivesCorrectmValuesByFreq(int do_gpu);