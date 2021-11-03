#include <unity.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "shapelet_basis.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */

#define UNITY_INCLUDE_FLOAT

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
  user_precision_t expected[] = { 1.000000, 0.140716, -0.637657, -0.466946, 0.222800,
                       0.583891, 0.311201, -0.232771, -0.522817, -0.352797,
                       0.081828, 0.426772, 0.454641, 0.185771, -0.180553,
                       -0.422827, -0.424534, -0.211821, 0.090262, 0.336831};

  //Check that they are equal
  for (int ind = 0; ind < num_values; ind++) {
    TEST_ASSERT_EQUAL_FLOAT((float)expected[ind], (float)sbf[test_indexes[ind]]);
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
