/*

*/


#include "calculate_visibilities_everybeam_common.h"

void setUp (void) {} /* Is run before every test, put unit init calls here. */
void tearDown (void) {} /* Is run after every test, put unit clean-up calls here. */


void do_test(){
  profile_lofar_everybeam(0, EB_LOFAR,
                   "../../../../test_installation/everybeam/LOFAR_HBA_MOCK.ms");
}


//Run the test with unity
int main(void)
{
  

  UNITY_BEGIN();

  RUN_TEST(do_test);

  return UNITY_END();

}
