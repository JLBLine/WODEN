#include "test_kern_calc_visi_common.h"

user_precision_t k_ref_stokesI[] = {28.0850901732358480,47.6022001610990344,13.1446050508080035,92.6177510299835234,62.9457456456704634,32.6574321930856470,9.4798732914494561,57.4176088590129510,6.6942051856086575,7.1798666254149817};

double k_ref_freqs[] = {5.0000000000000000e+07,7.7777777777777776e+07,1.0555555555555555e+08,1.3333333333333333e+08,1.6111111111111110e+08,1.8888888888888890e+08,2.1666666666666666e+08,2.4444444444444442e+08,2.7222222222222221e+08,3.0000000000000000e+08};

user_precision_t k_ref_curve_SIs[] = {34.0661738177012410,27.9416345804748012,4.0737385180043200,10.2020475527417478,-10.1852352751980622,-11.4724548440797207,15.9910003869060642,9.2318619436768898,-12.8687564415199187,21.6193252919712791};

user_precision_t k_ref_qs[] = {-0.9221355653482664,-0.7491819789226084,-0.1080190860621295,-0.2736834792726799,0.2724408259487379,0.3011619228028337,-0.4208711192002694,-0.2421777019147586,0.3390914894603383,-0.5779299448039859};

int k_num_list_values[] = {10,13,26,25,22};

int k_list_start_indexes[] = {0,10,23,49,74};

double k_list_freqs[] = {5.0000000000000000e+07,7.7777777777777776e+07,1.0555555555555555e+08,1.3333333333333333e+08,1.6111111111111110e+08,1.8888888888888890e+08,2.1666666666666666e+08,2.4444444444444442e+08,2.7222222222222221e+08,3.0000000000000000e+08,5.0000000000000000e+07,7.0833333333333328e+07,9.1666666666666657e+07,1.1250000000000000e+08,1.3333333333333333e+08,1.5416666666666666e+08,1.7500000000000000e+08,1.9583333333333331e+08,2.1666666666666666e+08,2.3750000000000000e+08,2.5833333333333331e+08,2.7916666666666663e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.0000000000000000e+07,7.0000000000000000e+07,8.0000000000000000e+07,9.0000000000000000e+07,1.0000000000000000e+08,1.1000000000000000e+08,1.2000000000000000e+08,1.3000000000000000e+08,1.4000000000000000e+08,1.5000000000000000e+08,1.6000000000000000e+08,1.7000000000000000e+08,1.8000000000000000e+08,1.9000000000000000e+08,2.0000000000000000e+08,2.1000000000000000e+08,2.2000000000000000e+08,2.3000000000000000e+08,2.4000000000000000e+08,2.5000000000000000e+08,2.6000000000000000e+08,2.7000000000000000e+08,2.8000000000000000e+08,2.9000000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.0416666666666664e+07,7.0833333333333328e+07,8.1250000000000000e+07,9.1666666666666657e+07,1.0208333333333333e+08,1.1250000000000000e+08,1.2291666666666666e+08,1.3333333333333333e+08,1.4375000000000000e+08,1.5416666666666666e+08,1.6458333333333331e+08,1.7500000000000000e+08,1.8541666666666666e+08,1.9583333333333331e+08,2.0625000000000000e+08,2.1666666666666666e+08,2.2708333333333331e+08,2.3750000000000000e+08,2.4791666666666666e+08,2.5833333333333331e+08,2.6875000000000000e+08,2.7916666666666663e+08,2.8958333333333331e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.1904761904761903e+07,7.3809523809523806e+07,8.5714285714285716e+07,9.7619047619047612e+07,1.0952380952380952e+08,1.2142857142857143e+08,1.3333333333333333e+08,1.4523809523809522e+08,1.5714285714285713e+08,1.6904761904761904e+08,1.8095238095238096e+08,1.9285714285714287e+08,2.0476190476190478e+08,2.1666666666666666e+08,2.2857142857142857e+08,2.4047619047619048e+08,2.5238095238095239e+08,2.6428571428571430e+08,2.7619047619047618e+08,2.8809523809523809e+08,3.0000000000000000e+08};

user_precision_t k_list_stokesI[] = {7.3334322214859675,7.7439115769393183,-2.4534036420851169,-4.5142903521680608,3.6718811417663204,7.7284896788444382,8.2806063500796672,7.0318252981874707,2.7650572489062561,1.2916153218546018,2.7527293921225109,-2.3478809511927770,-1.5692534332700312,9.3910716007455939,-4.4744602905959079,-2.7289403469062474,2.9095835881454919,8.4138548633403119,-1.8980918177073325,-1.9074070919935093,3.3648367959793077,3.1428506865350734,1.5457308789467525,-0.1166075235375956,6.5824117059065426,1.2645992730550812,7.2592316167161357,-1.8369880122656785,8.3628632204260605,9.4137692873846959,5.8178096799645758,6.7045601487781443,-4.3139903968134430,0.0539220351165817,8.0433550159070677,2.8891541033169679,-3.9222850468221830,0.0492599014758106,6.4322186092069451,-1.6039532195669288,9.6884131091036423,-3.4245532827523970,-1.1715363283273250,5.9723854889927015,-3.4667426279798748,8.6457937477458877,4.5189165618079041,-1.5894796125725295,-1.6472195603995479,9.2754971510147968,-1.3325889805345130,5.8324230292527020,1.8146888310350135,1.0759229229404301,3.6388767655208376,0.3014408460890126,0.8540327877351865,9.3862258030276671,7.4403840101851415,3.6713825186574116,0.5381824193284723,1.1359786562565644,4.5672056802461309,8.4334877890223474,-3.7600136320675932,-4.9457878671395754,5.9418372820862508,6.1478601670520145,9.9677798924598200,-2.0672473980617219,6.7251339366465661,-4.4042928427047610,2.9198170392191329,-3.8315296796642482,1.5184137771727482,0.2698018715063046,-3.0907195963789977,7.2747166436584827,2.7707235032590072,-0.2949571826543966,1.2878418299646697,5.4865893056358832,3.8238143150651496,-0.2391695870813031,3.5333605440492484,7.3178501228721373,2.1797345036560394,-2.6996842556340623,6.0927295044453817,1.5884155269966067,-4.3839047123429262,-2.1276676410443862,-0.4625634877686222,8.7682803241222729,-0.6898544822683439,6.7028203874275736};

user_precision_t k_list_stokesQ[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

user_precision_t k_list_stokesU[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

user_precision_t k_list_stokesV[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

int k_curve_comp_inds[] = {10,11,12,13,14,15,16,17,18,19};

int k_list_comp_inds[] = {20,21,22,23,24};

user_precision_t k_stokesV_pol_fracs[] = {-0.9817449998168577,-0.0800747001243891,-0.4667981531901007};

user_precision_t k_ref_stokesV[] = {0.2807806318184438,2.8094661033600721,2.1885406751495444,8.7206196693947611};

user_precision_t k_stokesV_power_SIs[] = {-0.4486011022372798,-0.3276087061731809,-0.0836999986307041,-0.7603315926664758};

user_precision_t k_stokesV_qs[] = {-1.3162875161337708,-1.4562530888604293,-1.4172677481297624,0.2495887949058320};

user_precision_t k_stokesV_curve_SIs[] = {-0.8192710027456049,-0.3742928352681743,0.8752868228104469,0.6798016360482528};

user_precision_t k_linpol_pol_fracs[] = {0.3215661231695712,0.1190820457509076,0.2662650520290260,0.0186837669462436,-0.2156996335573080};

user_precision_t k_ref_linpol[] = {2.5433630189430398,4.1962859091444091};

user_precision_t k_linpol_power_SIs[] = {-0.8464591879794217,-0.9220433157465378};

user_precision_t k_linpol_qs[] = {-1.6982837869692762,-0.4451547411803758};

user_precision_t k_linpol_curve_SIs[] = {0.6267541058497474,0.4389211873774903};

user_precision_t k_intr_pol_angle[] = {2.3566520273940732,4.7403995532856200,3.3713210828577527,2.0057071683313787,3.4729202270283879,0.6643848655602123,6.2604029929435363,5.5253508707736199,1.9911782891353793};

user_precision_t k_rms[] = {23.1701039856028430,53.1937314074958962,47.4008569609604677,12.7846061419278989,10.6669718594821994,72.2765065894846401,67.3684515811262514,68.5056511656220977,4.8686344972062390};

int k_stokesV_pol_frac_comp_inds[] = {11,16,17};

int k_stokesV_power_comp_inds[] = {22,3,10,2};

int k_stokesV_curve_comp_inds[] = {14,16,0,15};

int k_linpol_pol_frac_comp_inds[] = {18,8,16,2,1};

int k_linpol_power_comp_inds[] = {6,3};

int k_linpol_curve_comp_inds[] = {24,22};

int k_linpol_angle_inds[] = {18,8,16,2,1,6,3,24,22};

int total_num_flux_entires = 96;




void malloc_args_for_testing(args_for_testing_t *args_ft,
                             components_t *components,
                             int num_baselines,  int num_times,
                             int num_freqs, int num_components,
                             int n_powers, int n_curves, int n_lists,
                             int num_coeffs,
                             e_component_type component_type) {

  int num_visis = num_baselines*num_times*num_freqs;
  int num_beam_values = num_freqs*num_times*num_components;

  args_ft->num_baselines = num_baselines;
  args_ft->num_times = num_times;
  args_ft->num_freqs = num_freqs;
  args_ft->num_visis = num_visis;
  args_ft->num_beam_values = num_beam_values;

  args_ft->us = malloc(num_visis*sizeof(user_precision_t));
  args_ft->vs = malloc(num_visis*sizeof(user_precision_t));
  args_ft->ws = malloc(num_visis*sizeof(user_precision_t));
  args_ft->allsteps_wavelengths = malloc(num_visis*sizeof(user_precision_t));

  args_ft->extrap_freqs = malloc(num_freqs*sizeof(double));

  args_ft->primay_beam_J00 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  args_ft->primay_beam_J01 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  args_ft->primay_beam_J10 = malloc(num_beam_values*sizeof(user_precision_complex_t));
  args_ft->primay_beam_J11 = malloc(num_beam_values*sizeof(user_precision_complex_t));

  //Make sure arrays to hold summed visis are initialised to zero
  args_ft->sum_visi_XX_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_XX_imag = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_XY_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_XY_imag = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YX_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YX_imag = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YY_real = calloc(num_visis, sizeof(user_precision_t));
  args_ft->sum_visi_YY_imag = calloc(num_visis, sizeof(user_precision_t));

  //Setup the component information

  components->ras = malloc(num_components*sizeof(double));
  components->decs = malloc(num_components*sizeof(double));

  components->power_ref_stokesI = malloc(n_powers*sizeof(user_precision_t));
  // components->power_ref_stokesQ = malloc(n_powers*sizeof(user_precision_t));
  // components->power_ref_stokesU = malloc(n_powers*sizeof(user_precision_t));
  // components->power_ref_stokesV = malloc(n_powers*sizeof(user_precision_t));
  components->power_SIs = malloc(n_powers*sizeof(user_precision_t));
  components->power_ref_freqs = malloc(n_powers*sizeof(double));
  components->power_comp_inds = malloc(n_powers*sizeof(int));

  //GAUSS/SHAPELET STUFF
  components->pas = malloc(num_components*sizeof(user_precision_t));
  components->majors = malloc(num_components*sizeof(user_precision_t));
  components->minors = malloc(num_components*sizeof(user_precision_t));
  //
  components->ls = malloc(num_components*sizeof(double));
  components->ms = malloc(num_components*sizeof(double));
  components->ns = malloc(num_components*sizeof(double));

  //
  if (component_type == SHAPELET) {
    args_ft->u_shapes = malloc(num_components*num_baselines*num_times*sizeof(user_precision_t));
    args_ft->v_shapes = malloc(num_components*num_baselines*num_times*sizeof(user_precision_t));
    args_ft->w_shapes = malloc(num_components*num_baselines*num_times*sizeof(user_precision_t));
    args_ft->num_coeffs = num_coeffs;
    args_ft->sbf = malloc( sbf_N * sbf_L * sizeof(user_precision_t) );
    args_ft->sbf = create_sbf(args_ft->sbf);

    components->n1s = malloc(num_coeffs*sizeof(user_precision_t));
    components->n2s = malloc(num_coeffs*sizeof(user_precision_t));
    components->shape_coeffs = malloc(num_coeffs*sizeof(user_precision_t));
    components->param_indexes = malloc(num_coeffs*sizeof(user_precision_t));
  }

  components->n_stokesV_pol_frac = 0;
  components->n_stokesV_power = 0;
  components->n_stokesV_curve = 0;
  components->n_linpol_pol_frac = 0;
  components->n_linpol_power = 0;
  components->n_linpol_curve = 0;
  components->n_linpol_angles = 0;

}

void free_args_for_testing(args_for_testing_t *args_ft,
                           components_t components,
                           e_component_type component_type) {

  free( args_ft->us );
  free( args_ft->vs );
  free( args_ft->ws );
  free( args_ft->allsteps_wavelengths );
  free( args_ft->extrap_freqs );

  free( args_ft->sum_visi_XX_real );
  free( args_ft->sum_visi_XX_imag );
  free( args_ft->sum_visi_XY_real );
  free( args_ft->sum_visi_XY_imag );
  free( args_ft->sum_visi_YX_real );
  free( args_ft->sum_visi_YX_imag );
  free( args_ft->sum_visi_YY_real );
  free( args_ft->sum_visi_YY_imag );

  free( args_ft->primay_beam_J00 );
  free( args_ft->primay_beam_J01 );
  free( args_ft->primay_beam_J10 );
  free( args_ft->primay_beam_J11 );

  free( components.ras );
  free( components.decs );

  free( components.power_ref_stokesI );
  // free( components.power_ref_stokesQ );
  // free( components.power_ref_stokesU );
  // free( components.power_ref_stokesV );
  free( components.power_SIs );
  free( components.power_ref_freqs );
  free( components.power_comp_inds );

  free( components.pas );
  free( components.majors );
  free( components.minors );

  free( components.ls );
  free( components.ms );
  free( components.ns );

  if (component_type == SHAPELET) {
    free( args_ft->sbf );
    free( args_ft->u_shapes );
    free( args_ft->v_shapes );
    free( args_ft->w_shapes );
    free( components.n1s );
    free( components.n2s );
    free( components.shape_coeffs );
    free( components.param_indexes );
  }

  free( args_ft );

}

void create_lmn(components_t components) {
  int count = 0;
  for (int l_ind = 0; l_ind < 5; l_ind++) {
    for (int m_ind = 0; m_ind < 5; m_ind++) {
      components.ls[count] = -0.5 + 0.25*l_ind;
      components.ms[count] = -0.5 + 0.25*m_ind;
      components.ns[count] = sqrt(1 - components.ls[count]*components.ls[count] - components.ms[count]*components.ms[count]);

      count ++;
    }
  }
}


void setup_uvw_and_freqs(args_for_testing_t *args_ft, int num_times,
                         int num_freqs, int num_baselines){
  //Make up some u,v,w values and scale by wavelength in correct order
  int count = 0;
  float freq_base = 150e+6;
  float freq_inc = 25e+6;
  float wavelength, frequency;
  for ( int time_step = 0; time_step < num_times; time_step++ ) {
    for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
      frequency = freq_base + freq_step*freq_inc;
      wavelength = VELC / frequency;
      for (int baseline = 0; baseline < num_baselines; baseline++) {
        args_ft->us[count] = ((baseline + 1)*10) / wavelength;
        args_ft->vs[count] = ((baseline + 1)*10) / wavelength;
        args_ft->ws[count] = ((baseline + 1)*10) / wavelength;

        args_ft->allsteps_wavelengths[count] = wavelength;

        count ++;
      }
    }
  }

  for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
    args_ft->extrap_freqs[freq_step] = frequency = freq_base + freq_step*freq_inc;
  }
}





double _Complex visi_env_shape(int comp, int visi, int coeff,
                              args_for_testing_t *args_ft,
                              components_t components  ) {


  int num_baselines = args_ft->num_baselines;
  int num_times = args_ft->num_times;
  int num_freqs = args_ft->num_freqs;
  // int num_visis = args_ft->num_visis;

  int mod_baseline = visi - num_baselines*floorf((float)visi / (float)num_baselines);
  int time_ind = (int)floorf( (float)visi / ((float)num_baselines * (float)num_freqs));
  // int freq_ind = (int)floorf( ((float)visi - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);

  user_precision_t wavelength = args_ft->allsteps_wavelengths[visi];

  int uv_stripe = num_baselines*num_times*comp + time_ind*num_baselines + mod_baseline;

  double pa = (double)components.pas[comp];
  double sinpa = sin(pa);
  double cospa = cos(pa);

  double u_shape = args_ft->u_shapes[uv_stripe] / wavelength;
  double v_shape = args_ft->v_shapes[uv_stripe] / wavelength;

  double x = (cospa*v_shape + sinpa*u_shape); // major axis
  double y = (-sinpa*v_shape + cospa*u_shape); // minor axis

  //Scales the FWHM to std to match basis functions, and account for the
  //basis functions being stored with beta = 1.0
  //Basis functions have been stored in such a way that x is in the same
  //direction as on sky, but y is opposite, so include negative here
  double const_x = (components.majors[comp]*SQRT_M_PI_2_2_LN_2)/sbf_dx;
  double const_y = -(components.minors[comp]*SQRT_M_PI_2_2_LN_2)/sbf_dx;

  double _Complex Ipow_lookup[] = { 1.0 + I*0.0,
                           0.0 + I*1.0,
                          -1.0 + I*0.0,
                           0.0 + I*-1.0 };

  double xlow, xhigh, ylow, yhigh, u_value, v_value, f_hat;
  user_precision_t *sbf_n;

  // find the indices in the basis functions for u*beta_u and v*beta_v

  double xpos = x*const_x + sbf_c;
  double ypos = y*const_y + sbf_c;

  int xindex = (int)floor(xpos);
  int yindex = (int)floor(ypos);

  int n1 = (int)components.n1s[coeff];
  int n2 = (int)components.n2s[coeff];

  // if ( n1 < 0 || n2 < 0 || n1 >= sbf_N || n2 >= sbf_N ) continue;

  f_hat = components.shape_coeffs[coeff];
  //
  sbf_n = &args_ft->sbf[n1*sbf_L];
  xlow  = (double)sbf_n[xindex];
  xhigh = (double)sbf_n[xindex+1];
  u_value = xlow + (xhigh-xlow)*(xpos-xindex);

  sbf_n = &args_ft->sbf[n2*sbf_L];
  ylow  = (double)sbf_n[yindex];
  yhigh = (double)sbf_n[yindex+1];
  v_value = ylow + (yhigh-ylow)*(ypos-yindex);

  // accumulate the intensity model for baseline pair (u,v)
  double _Complex visi_env = 0.0 + I*0.0;
  visi_env = visi_env + Ipow_lookup[(n1+n2) % 4] * f_hat * u_value*v_value;

  return visi_env;
}

double _Complex visi_env_gauss(int comp, int visi,
                              args_for_testing_t *args_ft,
                              components_t components  ) {
  double pa = components.pas[comp];
  double u = args_ft->us[visi];
  double v = args_ft->vs[visi];
  double sinpa = sin(pa);
  double cospa = cos(pa);

  double x =  cospa*v + sinpa*u; // major axis
  double y = -sinpa*v + cospa*u; // minor axis
  double invsig_x = components.majors[comp];
  double invsig_y = components.minors[comp];

  double _Complex visi_env = exp( -0.5 * ( x*x*invsig_x*invsig_x*M_PI_2_2_LN_2 + y*y*invsig_y*invsig_y*M_PI_2_2_LN_2 ) ) + I*0.0;

  return visi_env;
}

/*
Basic implementation of the measurement equation to get expected visibilities
Loops over components, gets expected flux and beam gain and sum. DO everything
in double precision to test against
*/
void get_expected(int visi, int num_powers, int num_curves, int num_lists,
                  int num_baselines, int num_freqs, double *extrap_freqs,
                  int beamtype,
                  args_for_testing_t *args_ft,
                  components_t components,
                  e_component_type component_type,
                  double * expec_re_xx, double * expec_im_xx,
                  double * expec_re_xy, double * expec_im_xy,
                  double * expec_re_yx, double * expec_im_yx,
                  double * expec_re_yy, double * expec_im_yy) {
  double expec_re_inc, expec_im_inc;
  double temp; //, visi_freq;
  double gain, leak; //flux_ratio, ,
  int time_ind, freq_ind, beam_ind;

  double flux_I,flux_Q,flux_U,flux_V;

  * expec_re_xx = 0;
  * expec_im_xx = 0;
  * expec_re_xy = 0;
  * expec_im_xy = 0;
  * expec_re_yx = 0;
  * expec_im_yx = 0;
  * expec_re_yy = 0;
  * expec_im_yy = 0;

  int num_components = num_powers + num_curves + num_lists;

  double *expec_flux_I = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_Q = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_U = malloc(num_freqs*num_components*sizeof(double));
  double *expec_flux_V = malloc(num_freqs*num_components*sizeof(double));

  //Extrapolate all the fluxes for all components, it's a cheap enough calc
  CPU_extrapolate_fluxes_in_components(&components,
    num_powers, num_curves, num_lists,
    extrap_freqs, num_freqs,
    expec_flux_I, expec_flux_Q,
    expec_flux_U, expec_flux_V);

  for (int comp = 0; comp < num_components; comp++) {
    temp = 2*M_PI*( args_ft->us[visi]*components.ls[comp] + args_ft->vs[visi]*components.ms[comp] + args_ft->ws[visi]*(components.ns[comp]-1) );
    // sincos(temp, &(expec_im_inc), &(expec_re_inc));

    expec_im_inc = sin(temp);
    expec_re_inc = cos(temp);

    //Do some indexing so we can call up correct instrumental gain
    time_ind = (int)floorf( (float)visi / ((float)num_baselines * (float)num_freqs));
    freq_ind = (int)floorf( ((float)visi - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + comp;

    //Only testing Stokes I here, we test other pols in
    // flux_extrap = expec_flux_I[num_freqs*comp + freq_ind];
    flux_I = expec_flux_I[num_freqs*comp + freq_ind];
    flux_Q = expec_flux_Q[num_freqs*comp + freq_ind];
    flux_U = expec_flux_U[num_freqs*comp + freq_ind];
    flux_V = expec_flux_V[num_freqs*comp + freq_ind];

    // printf("%.5f %.5f %.5f\n",components.power_ref_freqs[comp],
    //                           components.power_SIs[comp],
    //                           components.power_ref_stokesI[comp] );

    double _Complex visi_env = 0.0 + I*0.0;

    if (component_type == POINT) {
      visi_env = 1.0 + I*0.0;
    }
    else if (component_type == GAUSSIAN) {
      visi_env = visi_env_gauss(comp, visi, args_ft, components);
      // visi_env = 1.0 + I*0.0;
    }
    else if (component_type == SHAPELET) {
      for (int coeff = 0; coeff < args_ft->num_coeffs; coeff++) {
        if (components.param_indexes[coeff] == comp) {
          // visi_env = visi_env_gauss(comp, visi, args_ft, components);
          visi_env += visi_env_shape(comp, visi, coeff, args_ft, components);

        }
      }
    }

    //Combine visi envelope with the measurement equation from above
    visi_env *= expec_re_inc + I*expec_im_inc;

    // printf("comp I Q U V %d %f %f %f %f\n", comp, flux_I, flux_Q, flux_U, flux_V);

    if (beamtype == NO_BEAM) {
      gain = 1.0;
      leak = 0.0;
    } else {
      //Righto in these tests, we've set all gains to the
      //same value. If the beam has leakages, they are same as well. That just
      //doubles the values; this is taken into account in `test_visi_outputs`
      //also have made the gains purely real to make life easier
      gain = creal(args_ft->primay_beam_J00[beam_ind]);
      leak = creal(args_ft->primay_beam_J01[beam_ind]);
    }

    // double _Complex expec_re_inc_xx;
    // double _Complex expec_im_inc_xx;
    // double _Complex expec_re_inc_xy;
    // double _Complex expec_im_inc_xy;
    // double _Complex expec_re_inc_yx;
    // double _Complex expec_im_inc_yx;
    // double _Complex expec_re_inc_yy;
    // double _Complex expec_im_inc_yy;
    double _Complex expec_inc_xx;
    double _Complex expec_inc_xy;
    double _Complex expec_inc_yx;
    double _Complex expec_inc_yy;
    double _Complex mult_xx, mult_xy, mult_yx, mult_yy;

    mult_xx = (gain*gain + leak*leak)*flux_I;
    mult_xx += (gain*gain - leak*leak)*flux_Q;
    mult_xx += (gain*leak + leak*gain)*flux_U;
    mult_xx += I*flux_V*(gain*leak - leak*gain);

    mult_xy = (gain*leak + leak*gain)*flux_I;
    mult_xy += (gain*leak - leak*gain)*flux_Q;
    mult_xy += (gain*gain + leak*leak)*flux_U;
    mult_xy += (I*flux_V)* (gain*gain - leak*leak);

    mult_yx = (leak*gain + gain*leak)*flux_I;
    mult_yx += (leak*gain - gain*leak)*flux_Q;
    mult_yx += (leak*leak + gain*gain)*flux_U;
    mult_yx += (I*flux_V)* (leak*leak - gain*gain);

    mult_yy = (leak*leak + gain*gain)*flux_I;
    mult_yy += (leak*leak - gain*gain)*flux_Q;
    mult_yy += (leak*gain + gain*leak)*flux_U;
    mult_yy += (I*flux_V)* (leak*gain - gain*leak);

    expec_inc_xx = visi_env*mult_xx;
    expec_inc_xy = visi_env*mult_xy;
    expec_inc_yx = visi_env*mult_yx;
    expec_inc_yy = visi_env*mult_yy;

    * expec_re_xx += creal(expec_inc_xx);
    * expec_im_xx += cimag(expec_inc_xx);
    * expec_re_xy += creal(expec_inc_xy);
    * expec_im_xy += cimag(expec_inc_xy);
    * expec_re_yx += creal(expec_inc_yx);
    * expec_im_yx += cimag(expec_inc_yx);
    * expec_re_yy += creal(expec_inc_yy);
    * expec_im_yy += cimag(expec_inc_yy);
  }

  free(expec_flux_I);
  free(expec_flux_Q);
  free(expec_flux_U);
  free(expec_flux_V);
}

//Take input parameters and test whether GPU outputs match expectations
void test_visi_outputs(int num_visis, int num_powers, int num_curves, int num_lists,
                       int num_baselines, int num_freqs, double *extrap_freqs,
                       e_beamtype beamtype,  args_for_testing_t *args_ft,
                       components_t components,
                       e_component_type component_type) {

  double frac_tol;
  //For some reason, AMD doesn't seem as accurate
  if (component_type == POINT) {
    #ifdef DOUBLE_PRECISION
      #ifdef __HIPCC__
        frac_tol = 1e-10;
      #else
        frac_tol = 1e-12;
      #endif
    #else
      frac_tol = 7e-5;
    #endif
  }
  else if (component_type == GAUSSIAN) {
    #ifdef DOUBLE_PRECISION
      #ifdef __HIPCC__
        frac_tol = 1e-10;
      #else
        frac_tol = 1e-11;
      #endif
    #else
      frac_tol = 7e-5;
    #endif
  }
  else {
    #ifdef DOUBLE_PRECISION
      #ifdef __HIPCC__
        frac_tol = 1e-10;
      #else
        frac_tol = 1e-11;
      #endif
    #else
      frac_tol = 5e-3;
    #endif
  }


  double expec_re_xx, expec_im_xx;
  double expec_re_xy, expec_im_xy;
  double expec_re_yx, expec_im_yx;
  double expec_re_yy, expec_im_yy;
  for (int visi = 0; visi < num_visis; visi++) {

      get_expected(visi, num_powers, num_curves, num_lists,
                  num_baselines, num_freqs, extrap_freqs,
                  beamtype, args_ft, components, component_type,
                  &expec_re_xx, &expec_im_xx,
                  &expec_re_xy, &expec_im_xy,
                  &expec_re_yx, &expec_im_yx,
                  &expec_re_yy, &expec_im_yy);

      // printf("%d %.16f\n", visi, expec_re_xx);

      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_xx, expec_re_xx, args_ft->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_xx, expec_im_xx, args_ft->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_xy, expec_re_xy, args_ft->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_xy, expec_im_xy, args_ft->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_yx, expec_re_yx, args_ft->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_yx, expec_im_yx, args_ft->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_yy, expec_re_yy, args_ft->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_yy, expec_im_yy, args_ft->sum_visi_YY_imag[visi]);

  }
}

/*
Test the __global__ code that calculates visibilities for different type
of COMPONENTs
Vary the l,m,n coords but keep all other variables constant
*/
void test_kern_calc_visi_Varylmn(e_beamtype beamtype, e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int n_powers = num_components;
  int n_curves = 0;
  int n_lists = 0;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, n_powers, n_curves, n_lists,
                          num_coeffs, comptype);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  //some beams have leakage, others don't; fill appropriately
  for (size_t visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;

    if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY ) {
      args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    }
    else {
      args_ft->primay_beam_J01[visi] = 0.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 0.0 + I*0.0;
    }
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    components.power_ref_stokesI[comp] = 1.0;
    // components.power_ref_stokesQ[comp] = 0.0;
    // components.power_ref_stokesU[comp] = 0.0;
    // components.power_ref_stokesV[comp] = 0.0;
    components.power_SIs[comp] = 0.0;
    components.power_ref_freqs[comp] = 150e+6;

    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);

    components.power_comp_inds[comp] = comp;

  }
  //Make up some u,v,w values and scale by wavelength in correct order
  setup_uvw_and_freqs(args_ft, num_times, num_freqs, num_baselines);

  // for (int freq_step = 0; freq_step < num_freqs; freq_step++) {
  //   printf("freq %.5e\n",args_ft->extrap_freqs[freq_step] );
  // }

  if (comptype == SHAPELET) {
    //Set the shapelet u,v,w same as the measurement equation one (this is not
    //true in reality but works fine for testing)

    // for (int comp = 0; comp < num_components; comp++) {
    //   for (int visi = 0; visi < num_visis; visi++) {
    //     args_ft->u_shapes[comp*num_visis + visi] = args_ft->us[visi];
    //     args_ft->v_shapes[comp*num_visis + visi] = args_ft->vs[visi];
    //     args_ft->w_shapes[comp*num_visis + visi] = args_ft->ws[visi];
    //   }
    // }

    int count = 0;

    for (int comp_step = 0; comp_step < num_components; comp_step++) {
      for ( int time_step = 0; time_step < num_times; time_step++ ) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {
          args_ft->u_shapes[count] = ((baseline + 1)*10);
          args_ft->v_shapes[count] = ((baseline + 1)*10);

          count ++;
        }
      }
    }

    //This means we have a single basis function for every component
    for (int coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }
  components.do_QUV = 0;
  test_kern_calc_visi_all(n_powers, n_curves, n_lists, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->extrap_freqs,
          args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, comptype);

  free_args_for_testing( args_ft, components, comptype );

}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and beam gains constant and vary the fluxes
Test works for all primary beam types
*/
void test_kern_calc_visi_VarylmnVaryFlux(e_beamtype beamtype,
                                         e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int n_powers = 10;
  int n_curves = 10;
  int n_lists = 5;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, n_powers, n_curves, n_lists,
                          num_coeffs, comptype);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;

    if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY ) {
      args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    }
    else {
      args_ft->primay_beam_J01[visi] = 0.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 0.0 + I*0.0;
    }
  }

  for (int pow_ind = 0; pow_ind < n_powers; pow_ind++) {
    components.power_ref_stokesI[pow_ind] = (pow_ind + 1);
    // components.power_ref_stokesQ[pow_ind] = 0.0;
    // components.power_ref_stokesU[pow_ind] = 0.0;
    // components.power_ref_stokesV[pow_ind] = 0.0;
    components.power_SIs[pow_ind] = -0.1*pow_ind;
    components.power_ref_freqs[pow_ind] = 150e+6;
    components.power_comp_inds[pow_ind] = pow_ind;
  }

  // //Just have some premade arrays living in the header and just stick em in
  // //for the CURVED_POWER_LAW and LIST components
  //
  components.curve_ref_stokesI = k_ref_stokesI;
  // components.curve_ref_stokesQ = k_ref_stokesQ;
  // components.curve_ref_stokesU = k_ref_stokesU;
  // components.curve_ref_stokesV = k_ref_stokesV;
  components.curve_SIs = k_ref_curve_SIs;
  components.curve_ref_freqs = k_ref_freqs;
  components.curve_qs = k_ref_qs;
  components.curve_comp_inds = k_curve_comp_inds;
  //
  components.num_list_values = k_num_list_values;
  components.list_start_indexes = k_list_start_indexes;
  components.list_comp_inds = k_list_comp_inds;
  components.list_freqs = k_list_freqs;
  components.list_stokesI = k_list_stokesI;
  components.list_stokesQ = k_list_stokesQ;
  components.list_stokesU = k_list_stokesU;
  components.list_stokesV = k_list_stokesV;

  components.total_num_flux_entires = total_num_flux_entires;


  //GAUSSIAN an SHAPELET models need major/minor, these are ignored when
  //testing POINT
  for (int comp = 0; comp < num_components; comp++) {
    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  setup_uvw_and_freqs(args_ft, num_times, num_freqs, num_baselines);

  if (comptype == SHAPELET) {
  //Set the shapelet u,v,w same as the measurement equation one (this is not
  //true in reality but works fine for testing)

    int count = 0;

    for (int comp_step = 0; comp_step < num_components; comp_step++) {
      for ( int time_step = 0; time_step < num_times; time_step++ ) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {
          args_ft->u_shapes[count] = ((baseline + 1)*10);
          args_ft->v_shapes[count] = ((baseline + 1)*10);

          count ++;
        }
      }
    }

    //This means we have a single basis function for every component
    for (size_t coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }

  //Polarisation models behbeh
  components.n_stokesV_pol_frac = 3;
  components.stokesV_pol_fracs = k_stokesV_pol_fracs;
  components.stokesV_pol_frac_comp_inds = k_stokesV_pol_frac_comp_inds;

  components.n_stokesV_power = 4;
  components.stokesV_power_ref_flux = k_ref_stokesV;
  components.stokesV_power_SIs = k_stokesV_power_SIs;
  components.stokesV_power_comp_inds =k_stokesV_power_comp_inds;
  
  components.n_stokesV_curve = 4;
  components.stokesV_curve_ref_flux = k_ref_stokesV;
  components.stokesV_curve_SIs = k_stokesV_curve_SIs;
  components.stokesV_curve_qs = k_stokesV_qs;
  components.stokesV_curve_comp_inds = k_stokesV_curve_comp_inds;

  components.n_linpol_pol_frac = 5;
  components.linpol_pol_fracs = k_linpol_pol_fracs;
  components.linpol_pol_frac_comp_inds = k_linpol_pol_frac_comp_inds;

  components.n_linpol_power = 2;
  components.linpol_power_ref_flux = k_ref_linpol;
  components.linpol_power_SIs = k_linpol_power_SIs;
  components.linpol_power_comp_inds = k_linpol_power_comp_inds;
  
  components.n_linpol_curve = 2;
  components.linpol_curve_ref_flux = k_ref_linpol;
  components.linpol_curve_SIs = k_linpol_curve_SIs;
  components.linpol_curve_qs = k_linpol_qs;
  components.linpol_curve_comp_inds = k_linpol_curve_comp_inds;

  components.n_linpol_angles = components.n_linpol_pol_frac + components.n_linpol_power + components.n_linpol_curve;
  components.intr_pol_angle = k_intr_pol_angle;
  components.rm_values = k_rms;
  components.linpol_angle_inds = k_linpol_angle_inds;
  components.do_QUV = 1;

  //Run the CUDA code
  test_kern_calc_visi_all(n_powers, n_curves, n_lists, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->extrap_freqs,
          args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}

/*
Test the __device__ code that updates the summed visibilities by grabbing the
correct beam gain and mesurement equation, multiplying and summing onto the visi
Here we keep the component visibilities and fluxes constant and vary the beam gains
Test works for all primary beam types
*/
void test_kern_calc_visi_VarylmnVaryBeam(e_beamtype beamtype,
                                         e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int n_powers = num_components;
  int n_curves = 0;
  int n_lists = 0;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, n_powers, n_curves, n_lists,
                          num_coeffs, comptype);


  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (size_t beam = 0; beam < num_beam_values; beam++) {
    args_ft->primay_beam_J00[beam] = beam + I*0.0;
    args_ft->primay_beam_J11[beam] = beam + I*0.0;

    if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY ) {
      args_ft->primay_beam_J01[beam] = beam + I*0.0;
      args_ft->primay_beam_J10[beam] = beam + I*0.0;
    }
    else {
      args_ft->primay_beam_J01[beam] = 0.0 + I*0.0;
      args_ft->primay_beam_J10[beam] = 0.0 + I*0.0;
    }
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    components.power_ref_stokesI[comp] = 1.0;
    // components.power_ref_stokesQ[comp] = 0.0;
    // components.power_ref_stokesU[comp] = 0.0;
    // components.power_ref_stokesV[comp] = 0.0;
    components.power_SIs[comp] = 0.0;
    components.power_ref_freqs[comp] = 150e+6;

    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);

    components.power_comp_inds[comp] = comp;

  }

  //Make up some u,v,w values and scale by wavelength in correct order
  setup_uvw_and_freqs(args_ft, num_times, num_freqs, num_baselines);

  if (comptype == SHAPELET) {
    //Set the shapelet u,v,w same as the measurement equation one (this is not
    //true in reality but works fine for testing)

    int count = 0;

    for (int comp_step = 0; comp_step < num_components; comp_step++) {
      for ( int time_step = 0; time_step < num_times; time_step++ ) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {
          args_ft->u_shapes[count] = ((baseline + 1)*10);
          args_ft->v_shapes[count] = ((baseline + 1)*10);

          count ++;
        }
      }
    }

    //This means we have a single basis function for every component
    for (int coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }
  components.do_QUV = 0;

  //Run the CUDA code
  test_kern_calc_visi_all(n_powers, n_curves, n_lists, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->extrap_freqs,
          args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}




void test_kern_calc_visi_VarylmnVaryPAMajMin(e_beamtype beamtype,
                                             e_component_type comptype) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int n_powers = num_components;
  int n_curves = 0;
  int n_lists = 0;
  int num_coeffs = num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, n_powers, n_curves, n_lists,
                          num_coeffs, comptype);

  //Setup l,m args that span a decent chunk of sky
  create_lmn(components);

  int num_beam_values = num_freqs*num_times*num_components;

  //Stick the gains to one everywhere
  for (int visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;

    if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY ) {
      args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    }
    else {
      args_ft->primay_beam_J01[visi] = 0.0 + I*0.0;
      args_ft->primay_beam_J10[visi] = 0.0 + I*0.0;
    }
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    components.power_ref_stokesI[comp] = 1.0;
    // components.power_ref_stokesQ[comp] = 0.0;
    // components.power_ref_stokesU[comp] = 0.0;
    // components.power_ref_stokesV[comp] = 0.0;
    components.power_SIs[comp] = 0.0;
    components.power_ref_freqs[comp] = 150e+6;

    components.power_comp_inds[comp] = comp;

    //Set major,minor to 3 arcmins
    components.pas[comp] = (comp + 1)*DD2R;
    components.majors[comp] = (comp + 1)*(DD2R / 60.0);
    components.minors[comp] = (comp + 2)*(DD2R / 60.0);
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  setup_uvw_and_freqs(args_ft, num_times, num_freqs, num_baselines);

  if (comptype == SHAPELET) {
    //Set the shapelet u,v,w same as the measurement equation one (this is not
    //true in reality but works fine for testing)

    int count = 0;

    for (int comp_step = 0; comp_step < num_components; comp_step++) {
      for ( int time_step = 0; time_step < num_times; time_step++ ) {
        for (int baseline = 0; baseline < num_baselines; baseline++) {
          args_ft->u_shapes[count] = ((baseline + 1)*10);
          args_ft->v_shapes[count] = ((baseline + 1)*10);

          count ++;
        }
      }
    }

    //This means we have a single basis function for every component
    for (int coeff = 0; coeff < num_coeffs; coeff++) {
      components.n1s[coeff] = 0.0;
      components.n2s[coeff] = 0.0;
      components.shape_coeffs[coeff] = 1.0;
      components.param_indexes[coeff] = coeff;
    }
  }
  components.do_QUV = 0;
  //
  test_kern_calc_visi_all(n_powers, n_curves, n_lists, num_baselines, num_coeffs,
          num_freqs, num_visis, num_times, beamtype, comptype,
          components, args_ft->extrap_freqs,
          args_ft->us, args_ft->vs, args_ft->ws,
          args_ft->u_shapes, args_ft->v_shapes,
          args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
          args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
          args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
          args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
          args_ft->allsteps_wavelengths, args_ft->sbf,
          args_ft->primay_beam_J00, args_ft->primay_beam_J01,
          args_ft->primay_beam_J10, args_ft->primay_beam_J11);
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}
