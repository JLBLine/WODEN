#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "shapelet_basis.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#ifdef DOUBLE_PRECISION
  double TOL = 1e-15;
#else
  double TOL = 1e-7;
#endif

/*
create_sbf just creates a massive array of shapelet basis functions
Here just test that a few array locations return the correct values
*/

void test_create_sbf_GivesCorrectValues(void)
{
  //Create the shapelet basis function array
  user_precision_t *sbf;
  sbf = NULL;
  sbf = malloc( sbf_N * sbf_L * sizeof(user_precision_t) );
  sbf = create_sbf(sbf);

  int num_values = 20;

  //Indexes to test in sbf
  int test_indexes[] = {5000, 15011, 25022, 35033, 45044, 55055,
                        65066, 75077, 85088, 95099, 105110, 115121,
                        125132, 135143, 145154, 155165, 165176,
                        175187, 185198, 195209};

  //Matching expected values
  double expected[] = {1.00000000, 0.14071601, -0.63765672, -0.46694581,
                       0.22280003, 0.58389053, 0.31120116, -0.2327714,
                       -0.52281663, -0.35279668, 0.081828458, 0.42677188,
                       0.45464097, 0.18577076, -0.18055274, -0.42282655,
                       -0.42453389, -0.21182068, 0.09026158, 0.33683096};

  //Check that they are within tolerance
  for (int ind = 0; ind < num_values; ind++) {
    TEST_ASSERT_DOUBLE_WITHIN(TOL, expected[ind], sbf[test_indexes[ind]]);
  }
  free(sbf);
}

//Run test using unity
int main(void)
{
    UNITY_BEGIN();

    RUN_TEST(test_create_sbf_GivesCorrectValues);

    return UNITY_END();
}
