#define N_LNA_FREQS 451
#define MAX_ZMATRIX_FREQS 256
#define NUM_DIPOLES 16
#define SIZ_Z_MATRIX (NUM_DIPOLES*2*NUM_DIPOLES*2)
#define VEL_LIGHT 299792458.0
#define DQ (435e-12*VEL_LIGHT)  // delay quantum of the MWA beamformer in meters.
#define MAX_POLS 4
#define N_COPOL  2
#define NUM_PARAMS_PER_DIPOLE 2

#include "hdf5.h"
#include "hdf5_hl.h"

typedef struct _copy_primary_beam {
  int type;                     //!< id of the primary beam function type
  int dipole_type;              //!< id of the primary beam function type
  int num_parameters;
  // coord_PNT_t coord_pnt;        //!< un-phased boresight coordinates for the tile [radians]
  float ref_az,ref_za;          //!< pointing direction of the primary beam [radians]
  float dpl_sep;                //!< ideal spacing between dipoles (centre to centre) in tiles [metres]
  float dpl_hgt;                //!< ideal height of dipole above the ground plane [metres]
  float delay[NUM_DIPOLES];     //!< delay for each dipole (sec). Same for both pol
  float *parameters;
  float ref_freq;               //!< if used to generate reference J matrices, some applications need a frequency.
  double _Complex **Q1, **Q2;      //!< Beam modes used for Spherical Harmonic model
  double _Complex **p_T, **p_P;    //!< Some pre-computed values used in tile response
  double **M,**N;
  int nmax;
  int nMN;
  float _Complex norm_fac[MAX_POLS];

  // BP 2019: All the Spherical Harmonic Beam data are double
  // so we will use them on the GPUs as well or there will be all kinds
  // of issues with copying

  float _Complex *d_Q1, *d_Q2;
  float _Complex *d_p_T, *d_p_P; // Precomputed values used in CPU optimization
  float *d_M, *d_N;
  int d_nmax;
  int d_nMN;
  float _Complex d_norm_fac[MAX_POLS];

  float _Complex *emn_P;
  float _Complex *emn_T;

  float _Complex *d_emn_T_sum;
  float _Complex *d_emn_P_sum;

  float _Complex *rts_P1;
  float _Complex *rts_P_sin;

} copy_primary_beam_t;

typedef struct _copy_cal_context_t {

  // vis_set_t *bandvis;
  // source_catalogue_t *srccat;
  // source_list_t *source_list;
  //
  // int nSources;
  // int nStations;
  //
  // // These contain PP and QQ gains for each source and tile.
  // // *** Note: these have been superceeded by the GainMatrices, and may become defunct in future versions.
  // double _Complex ***current_gain_models;	// the current gain estimates (after averaging, etc.).
  // double _Complex ***updated_gain_models;	// the working estimates.
  // double _Complex ***true_gain_models;		// the true estimates, if available (for simulated data tests).

  // // These contain a Jones matrix for each source and tile.
  // float _Complex ***TileGainMatrices;		// the current Jones matrices (after averaging, etc.).
  // float _Complex ***NewGainMatrices;		// the working Jones matrices.
  // float _Complex ***VelocityMatrices;		// extra derivties for the Jones matrices (a la kalman filters)
  // float _Complex ***TrueGainMatrices;		// the true Jones matrices, if available (for simulated data tests).
  //
  // // This contains a coherency matrix for each source.
  // float _Complex **PeeledCoherencyMatrices;      // the acucmulated coherency matrices.
  // float _Complex **PeeledCoherencyWeights;       // the acucmulated coherency matrices.
  //
  // // Direction-independent gain correction matrices.
  // float _Complex *PostAlignmentMatrix;		// Model Jones matrix to align to in the reference direction.
  // float _Complex **PreAlignmentMatrices;		// An estimate of each tiles Jones matrix in some reference direction.
  // float _Complex **AlignmentChanges;		// Ratios of new to old PreAlignmentMatrices for quick updates.
  // float _Complex **DIJonesMatrices;		// An estimate of each tiles direction-independent Jones matrix.
  // float _Complex **NewDIMatrices;		// Updated DIJonesMatrices. Leave update until after unpeeling.
  // float AlignmentFluxDensity;			// Flux Density in the reference direction. Used when printing stats.
  //
  // int **fit_flag;

  // These contain primary beam model parameters, etc., for each tile.
  /* Each is an array with elements for each frequency which point to an array of primary_beam_t for each antenna */
  copy_primary_beam_t *primary_beams, *true_beams;

  // Note that these are only single precisoin
  float _Complex *d_CurrGainModels, *d_UpdateGainModels, *d_TrueGainModels;
  float _Complex *d_TileGainMatrices, *d_NewGainMatrices, *d_VelocityMatrices, *d_TrueGainMatrices, *d_TmpMatrices;
  float _Complex *d_PostAlignmentMatrix, *d_PreAlignmentMatrices, *d_AlignmentChanges;
  float _Complex *d_DIJonesMatrices, *d_NewDIMatrices;
  //! Mirror of primary_beams[*]->dpl_sep
  float *d_primaryBeamDplSep;
  //! Mirror of primary_beams[*]->dpl_hgt
  float *d_primaryBeamDplHgt;
  //! Mirror of primary_beams[*]->coord_pnt.ha
  float *d_primaryBeamCoordPntHA;
  //! Mirror of primary_beams[*]->coord_pnt.dec
  float *d_primaryBeamCoordPntDec;
  //! Mirror of primary_beams[*]->parameters[*]
  float *d_primaryBeamParams;
  // FEE Beam model parameters;
  //! Mirror of primary_beams[*]->M[*]
  float *d_M;
  //! Mirror of primary_beams[*]->N[*]
  float *d_N;
  //! Mirror of primary_beams[*]->Q1[*]
  float *d_Q1;
  //! Mirror of primary_beams[*]->Q2[*]
  float *d_Q2;
  //! Mirror of primary_beams[*].norm_fac. This is the same for all beams
  float _Complex *d_norm_fac;


} copy_cal_context_t;

//! RTS high level options that can be changed by command line input
typedef struct _copy_rts_options {

    /* global settings */
    int n_chan,n_iter,do_rts;
    char *logfile_path;
    float max_freq;             //!< max freq (Hz) that the RTS is currently running with
    float bandwidth;            //!< bandwidth here, rather than pass a separate variable around

    float cal_baseline_min;
    float cal_baseline_max;
    float cal_short_baseline_sigma;
    float cal_long_baseline_sigma;

    float iono_baseline_min;
    float iono_baseline_max;
    float iono_short_baseline_sigma;
    float iono_long_baseline_sigma;

    int update_cal_amplitudes;

    float baseline_min;
    float baseline_max;
    float short_baseline_sigma;
    float uv_taper_sigma;
    float image_oversampling;

    int undo_AA_function;
    int image_PSF;
    int do_uv_sampling_accum;
    int use_accum_uv_sampling;
    int do_uniform_weighting;
    int do_robust_weighting;
    float robustness;

    int assume_conjugate_symmetry;

    /* general options */
    float import_base_freq;
    int import_base_timestep;
    char *import_base_filename;
    char *TileFlags_filename;       // file contain tiles to skip. Tiles start at 0.
    char *ChannelFlags_filename;    // file contain frequency channels to skip. Channels start at 0.
    int obs_type;
    float vis_scale;

    int do_raw_data_corrections;
    int do_MWA_rx_corrections;
    int do_RFI_flagging;
    int use_fast_PB_models;
    int use_extended_dipoles_models;
    int use_stored_cal_files;
    int apply_DI_calibration;
    int generate_DI_jones;
    int atomic_MAPS_vis;

    int array_type;

    int beam_type;
    char *hdf5_filename;

    /* Simulator options */
    int sim_mode;
    int n_sim;
    char *simcat_filename;

    /* CML options */
    int do_cml,n_cal,n_iono,n_peel,n_prepeel,fit_iono,fit_beams;
    char *sourcecat_filename, *beam_param_files, *primary_calibrator;
    int store_cal_data;
    int disable_srclist_vetos;

    /* gridding & imaging options */
    int make_image,do_grid;
    float array_size_meters;
    float fov_degs;         //!< field of view in degrees. 2 radians (114 degs) does all-sky
    int gridding_method;

    /* regridding options */
    int make_stokes_snapshots;      //!< convert the intrumental pol output into Stokes.
    int make_weighted_snapshots;    //!< convert the intrumental pol output into weighted snapshots.
    int do_regrid;                  //!< regrid the data into a storage frame.
    int projection;                 //!< The regridding projection enum in regrid_engine.h.
    int wrap;                       //!< if image goes across 12hrs 0 images the <12 1 images > 12
    int regrid_method;              //!< enum in regrid_engine.h
    int gen_pix_matrices;           //!< generate the matrices needed for either make_stokes_snapshots or integrate
    int store_pix_matrices;         //!< store the weight matrices in the output fits files
    int iterate_over_all_gains;     //!< iterate around the entire calibrator list when finding Jones matrices (iter=0)
    int start_integration_at;       //!< cadence at which to start integrating images
    int start_processing_at;        //!< cadence at which to start processing calibration and imaging

    /* analysis options */
    int analyse_image;      //!< calc statistics in UV plane.

    /* wide field correction options */
    double point_centre_ha;
    double point_centre_dec;
    double image_centre_ra;
    double image_centre_dec;
    //
    // /* shared memory */
    // int use_shared;
    // /* Integrate */
    // int integrate;
    // int fscrunch;
    // /* NSIDE */
    // int nside;
    // /* debugging */
    // int get_tiles;
    //
    // int do_mfs_gridding;
    // int do_initial_bandpass_fit;
    //
    // int swap_pols;
    //
    // /* read time/freq data from a single file */
    // int read_all;
    // int add_node_number_to_filename;
    // int realtime_read_mwac;

    /* CUDA */
    // int cudaDevice;
    // int number_of_devices;
    //
    // /* precomputation */
    // char *master_ipv4;
    //
    // /* Adding options that were separate variables */
    // int n_CorrDumpsPerCadence;
    // float CorrDumpTime;
    // int n_IntegrationLengths;
    // float *integ_times;

    /** Flags for whether to use packets, and whether they should be generated or received from correlator */
    // int use_packets;
    // int use_correlator;
    // int correlator_port;
    // int use_threaded_vi;
    // int skip_primary_header;

    //Sets the acceptable time decorrelation limit - defaults to 1%
    // float set_timedecor;

    /** Extra parameters that were taken from uvfits headers, to be transferred to nodes */
    // array_spec_t arr_spec;
    // obs_context_t context;
    // char *array_file;

    int delays_given;
    float *dipole_delays;
    int *dead_dip;
  int disable_dipole_flags;

  //   int mpi_node;
  //
  //   // these are used to read non-standard mpi_node numbers, but not distributed. Just used to set mpi_node.
  //   int num_IDs;
  //   int *nodeIDs;
  //
  // int importCotterFlags; // flag visibilities using Cotter Flags
  // char *importCotterfilename; // Path for Cotter produced mwaf files
  //
  //
  // int interleaveTimeSamples; // Image even/odd time samples for Bryna Hazelton's image based Power Spectrum code (0:no interleaving, 1: odd samples, 2: even samples)
  //
  // // For direct reading of visibility from correlator gpubox files
  // int read_gpubox_direct;


    // // shapelet basis functions.
    // int sbf_N;
    // int sbf_L;
    // int sbf_c;
    // float sbf_dx;
    // float *sbf;
    //
    // // shapelet2 basis functions.
    // int sbf_N2;
    // int sbf_L2;
    // int sbf_c2;
    // float *sbf2;

    // uvfits output options
    // int writeVisToUVFITS;

    // write peeled models as images or uvfits
  // int write_vis_models;

    //metafits file options

    int readMetafitsFile;
    char *metafitsFilename;
    int *cott_ant_map;

} copy_rts_options_t;

/***************************************************************
Data Structure used by H5Literate.
***************************************************************/
struct opdata {
  float freq_in;
  float freq_out;
  float least_diff;
};

// enumerated type for antenna tile models
// typedef enum {
//     REF_CROSSED_DIPOLES_ON_GROUNDPLANE ,
//     ANT_CROSSED_DIPOLES_ON_GROUNDPLANE ,
//     ANT_CROSSED_DIPOLES_ON_GROUNDPLANE_32T ,
//     ANT_CROSSED_DIPOLES_LWA_MEMO_175 ,
//     ANT_SIMULATED
// } ant_types;

/*
 * Operator function to be called by H5Literate.
 */
herr_t RTS_op_func (hid_t loc_id, const char *name, const H5L_info_t *info,
            void *operator_data);

// void RTS_reorderDelays(float in[NUM_DIPOLES], float out[NUM_DIPOLES]);
//
// void RTS_reorderDelays2RTS(float in[NUM_DIPOLES], float out[NUM_DIPOLES]);

int RTS_HDFBeamInit(char *h5filename, float freq_Hz, copy_primary_beam_t *pb,
                float *FEE_delays, int stn);

void RTS_P1SIN(int nmax, float theta, double _Complex **P_sin, double _Complex **P1);
//
int RTS_get_FF2(float phi, float theta, copy_primary_beam_t *pb, float _Complex result[4], int scaling, int clean_up);
//
int RTS_getJonesSphHarm(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);
//
int RTS_getJonesDipole(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4],int scaling);
//
int RTS_getTileResponse(float freq_Hz, float az, float za, copy_primary_beam_t *pb, float _Complex result[4], int scaling, float rotation);
//
void RTS_freeHDFBeam(copy_primary_beam_t *pb);
//
int RTS_getHDFBeamNormalization(char *h5filename, float freq, copy_primary_beam_t *pb);
//
void RTS_HDFBeamCleanUp();
