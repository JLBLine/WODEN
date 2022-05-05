/*******************************************************************************
*  Methods to read in the hyperdrive style sky model. Yaml code based on this
* useful example here https://www.wpsoftware.net/andrew/pages/libyaml.html
* big shoutout to Andrew Poelstra!
*
*  @author J.L.B. Line
*
*  Please see documentation in https://woden.readthedocs.io/en/latest/index.html
*******************************************************************************/
#pragma once

#include <stdio.h>
#include <yaml.h>
#include "constants.h"
#include "woden_struct_defs.h"


/*!
A struct to contain how many elements we have malloc'd for each set of COMPONENTS
*/
typedef struct _track_comp_malloc_t {
  int n_comps; /*!< How many components have been malloced for */
  int n_powers; /*!< How many POWER_LAW have been malloced for */
  int n_curves; /*!< How many CURVED_POWER_LAW have been malloced for */
  int n_lists; /*!< How many LIST components have been malloced for */
  int n_list_values; /*!< How many LIST flux values have been malloced for */
  int n_shape_coeffs; /*!< How many SHAPELET coeffs have been malloced for */

} track_comp_malloc_t;



/**
 @brief Takes a path to `hyperdrive`-style sky model and populates a
 `source_catalogue_t` struct with the contents of `filename`.

 @details The `hyperdrive` sourcelist at `filename` should contain a number of
 sources, with the basic structure:

            singlepoint_power:
            - ra: 30.0
              dec: -30.0
              comp_type: point
              flux_type:
                power_law:
                  si: -0.8
                  fd:
                    freq: 150000000.0
                    i: 2.0
                    q: 0.0
                    u: 0.0
                    v: 0.0

For a more detailed explanation of the sky model, please see the
`hyperdrive` documentation https://github.com/MWATelescope/mwa_hyperdrive/wiki/Source-lists

Note that `srccat` should be memory intialised, so declare with something
like `source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );`
before feeding into this function.

@param[in] *yaml_path Path to a hyperdrive-style sky model
@param[in] *srccat Struct to contain sky model information.

@return Integer where 0 if read was successful, 1 if failed
 */
int read_yaml_skymodel(const char *yaml_path, source_catalogue_t *srccat);
