/*! \file
  Methods to read in and crop a WODEN style sky model to everything above
  the horizon.

  @author J.L.B. Line
*/
#pragma once
#include "woden_struct_defs.h"

//Somthing WEIRD is up with the way I'm using the documentation package
//'breathe', so having to include the definition value in the documentation
//string to get values to appear (sigh)

/**"SOURCE" - Line beginning / containing SOURCE information (includes number of COMPONENTS) */
#define SRC_KEY         "SOURCE"
/**"ENDSOURCE" - Line ending SOURCE information */
#define SRC_END         "ENDSOURCE"
/**"COMPONENT" - Line beginning / containing COMPONENT information */
#define COMP_KEY        "COMPONENT"
/**"ENDCOMPONENT" - Line ending COMPONENT information */
#define COMP_END        "ENDCOMPONENT"
/**"FREQ" (Deprecated) - Lines contains FREQ information */
#define FREQ_KEY        "FREQ"
/**"LINEAR" - Line contains simple exponential SED information */
#define LINEAR_KEY      "LINEAR"
/**"POINT" - Line contains POINT RA, Dec information */
#define POINT_KEY       "POINT"
/**"GAUSSIAN" - Line contains GAUSSIAN RA, Dec information */
#define GAUSSIAN_KEY    "GAUSSIAN"
/**"GPARAMS" - Line contains GAUSSIAN major, minor, PA information */
#define GPARAMS_KEY     "GPARAMS"
/**"SHAPELET" - Line containing SHAPELET RA, Dec information */
#define SHAPELET_KEY    "SHAPELET"
/**"SPARAMS" - Line contains SHAPELET beta1 (major), beta2 (minor), PA information  */
#define SPARAMS_KEY     "SPARAMS"
/**"SCOEFF" - Line contains SHAPELET basis numbers n1, n2m coefficient information */
#define SCOEFF_KEY      "SCOEFF"

/**
 @brief Takes a path to WODEN-style sky model and populates a
 `source_catalogue_t` struct with the contents of `filename`.

 @details The WODEN sourcelist at `filename` should contain a number of
 sources, with the basic structure:

            SOURCE source_name P 1 G 0 S 0 0
            COMPONENT POINT 4.0 -27.0
            LINEAR 1.8e+08 10.0 0 0 0 -0.8
            ENDCOMPONENT
            ENDSOURCE

For a more detailed explanation of the sky model, please see the
documentation at @todo Link the online documentation when there is a link.

Note that `srccat` should be memory intialised, so declare with something
like `source_catalogue_t *raw_srccat = malloc( sizeof(source_catalogue_t) );`
before feeding into this function.

@param[in] *filename Path to a WODEN-style sky model
@param[in] *srccat Struct to contain sky model information.

@return Integer where 0 if read was successful, 1 if failed
 */
int read_text_skymodel(const char *filename, source_catalogue_t *srccat);
