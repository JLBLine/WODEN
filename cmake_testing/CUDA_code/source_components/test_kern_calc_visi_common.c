#include "test_kern_calc_visi_common.h"

user_precision_t k_ref_stokesI[] = {28.0850901732358480,47.6022001610990344,13.1446050508080035,92.6177510299835234,62.9457456456704634,32.6574321930856470,9.4798732914494561,57.4176088590129510,6.6942051856086575,7.1798666254149817};

user_precision_t k_ref_stokesQ[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

user_precision_t k_ref_stokesU[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

user_precision_t k_ref_stokesV[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};


double k_ref_freqs[] = {5.0000000000000000e+07,7.7777777777777776e+07,1.0555555555555555e+08,1.3333333333333333e+08,1.6111111111111110e+08,1.8888888888888890e+08,2.1666666666666666e+08,2.4444444444444442e+08,2.7222222222222221e+08,3.0000000000000000e+08};

user_precision_t k_ref_curve_SIs[] = {55.7708105478831513,33.6649010803953175,28.5228944413592060,-0.7428462271528363,76.1828870445264812,16.5299548910859144,0.7473460208037787,47.8928376930352115,23.7844860951525696,15.0469139365674849};

user_precision_t k_ref_qs[] = {-1.4753683681601231,-0.9053846979507472,-0.7524507745093478,0.0195263350708341,-1.9929605977985849,-0.4328242525295414,-0.0195514444489788,-1.2658870691051862,-0.6227932660710171,-0.4007292504199662};

int k_num_list_values[] = {29,9,19,31,21};

int k_list_start_indexes[] = {0,29,38,57,88};

double k_list_freqs[] = {5.0000000000000000e+07,5.8928571428571433e+07,6.7857142857142866e+07,7.6785714285714284e+07,8.5714285714285716e+07,9.4642857142857149e+07,1.0357142857142857e+08,1.1250000000000000e+08,1.2142857142857143e+08,1.3035714285714287e+08,1.3928571428571430e+08,1.4821428571428573e+08,1.5714285714285713e+08,1.6607142857142860e+08,1.7500000000000000e+08,1.8392857142857143e+08,1.9285714285714287e+08,2.0178571428571430e+08,2.1071428571428573e+08,2.1964285714285716e+08,2.2857142857142860e+08,2.3750000000000000e+08,2.4642857142857143e+08,2.5535714285714287e+08,2.6428571428571430e+08,2.7321428571428573e+08,2.8214285714285719e+08,2.9107142857142860e+08,3.0000000000000000e+08,5.0000000000000000e+07,8.1250000000000000e+07,1.1250000000000000e+08,1.4375000000000000e+08,1.7500000000000000e+08,2.0625000000000000e+08,2.3750000000000000e+08,2.6875000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.3888888888888888e+07,7.7777777777777776e+07,9.1666666666666657e+07,1.0555555555555555e+08,1.1944444444444445e+08,1.3333333333333333e+08,1.4722222222222221e+08,1.6111111111111110e+08,1.7500000000000000e+08,1.8888888888888890e+08,2.0277777777777776e+08,2.1666666666666666e+08,2.3055555555555555e+08,2.4444444444444442e+08,2.5833333333333331e+08,2.7222222222222221e+08,2.8611111111111110e+08,3.0000000000000000e+08,5.0000000000000000e+07,5.8333333333333336e+07,6.6666666666666664e+07,7.5000000000000000e+07,8.3333333333333328e+07,9.1666666666666657e+07,1.0000000000000000e+08,1.0833333333333333e+08,1.1666666666666666e+08,1.2500000000000000e+08,1.3333333333333333e+08,1.4166666666666666e+08,1.5000000000000000e+08,1.5833333333333331e+08,1.6666666666666666e+08,1.7500000000000000e+08,1.8333333333333331e+08,1.9166666666666666e+08,2.0000000000000000e+08,2.0833333333333331e+08,2.1666666666666666e+08,2.2500000000000000e+08,2.3333333333333331e+08,2.4166666666666666e+08,2.5000000000000000e+08,2.5833333333333331e+08,2.6666666666666666e+08,2.7500000000000000e+08,2.8333333333333331e+08,2.9166666666666663e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.2500000000000000e+07,7.5000000000000000e+07,8.7500000000000000e+07,1.0000000000000000e+08,1.1250000000000000e+08,1.2500000000000000e+08,1.3750000000000000e+08,1.5000000000000000e+08,1.6250000000000000e+08,1.7500000000000000e+08,1.8750000000000000e+08,2.0000000000000000e+08,2.1250000000000000e+08,2.2500000000000000e+08,2.3750000000000000e+08,2.5000000000000000e+08,2.6250000000000000e+08,2.7500000000000000e+08,2.8750000000000000e+08,3.0000000000000000e+08};

user_precision_t k_list_stokesI[] = {-1.4845370362285797,-2.2780003826346773,2.0303307674092155,6.6335503574423029,1.0729428554986278,2.2039535928664433,-1.4675282484851406,-0.5928747437746811,7.4072927386772101,7.3464653995738320,0.6690887768437230,2.7155999834918418,6.8285391451154371,-0.8977250968026249,-1.7375185331625760,-1.5036064887785736,8.4975327694349936,-3.6445325205920369,-0.3071962645113064,9.0646511710783511,7.5985122703618959,-2.1417576121034991,-3.1476151830412240,9.3557419924686904,0.3268964680866624,-1.3246146495669904,-4.5066688900232013,-4.2992408155284076,8.7323366317905968,0.4822075959835921,4.5498413600282213,-0.0206434022602453,0.3063026095628043,4.2221111590128331,-0.4447156553240106,5.0195752779618754,0.6280617036800518,-1.5980898708618039,7.1970916501953290,-2.5314675447335229,0.6250474295187942,4.4012917942692074,-0.6728521282671798,5.8741265300500345,3.0664661591418021,-3.3427309998308430,2.6381373186249988,9.6736462307284210,9.2967099033434568,7.6913779157763855,3.5535304298916905,8.1434113750043249,-3.5845066004927868,-0.2424782144328921,3.8065869257996248,-0.5128255663374786,1.6619652793848942,1.8778215445098132,5.2621916301664413,2.0850783909460535,2.8218483046981859,1.0674480815894558,0.3955860441339878,2.2940858083516771,-3.1812633020459051,0.0104220541904256,-1.6763448314343332,8.9686745512421133,2.5508087825081134,3.3723898911254615,2.3957909711618361,9.4727317907747519,2.0643816382963749,8.5678468694387160,-1.0994164729731182,-3.6057419690801980,-3.1897027218156571,4.3290715529177461,7.2006557938731053,5.7919089053311765,8.0739309200452478,2.8172546046374949,8.1553782416704923,2.3487009092268103,7.0929958489892844,6.6623781813219569,1.0541075806542661,9.7525317609459723,-2.5724008632828790,-4.0613592480613390,6.2657882283634621,9.0131286454391493,9.9034897511126818,4.0371749069316021,5.5926748595404021,1.0850378279279873,4.2481153505004361,4.6297122636204229,2.0237867016228748,-0.7529873044172444,8.7090922019704138,6.2704486145533771,6.1166300620556697,3.0113902925759817,-1.8889785948873472,2.4693593706517989,8.9633797107176818,-2.4956858613516353,3.7287062212609747};

user_precision_t k_list_stokesQ[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

user_precision_t k_list_stokesU[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

user_precision_t k_list_stokesV[] = {0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000};

int k_curve_comp_inds[] = {10,11,12,13,14,15,16,17,18,19};
// int k_curve_comp_inds[] = {15,16,17,18,19,20,21,22,23,24};

int k_list_comp_inds[] = {20,21,22,23,24};

int total_num_flux_entires = 109;


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
  components->power_ref_stokesQ = malloc(n_powers*sizeof(user_precision_t));
  components->power_ref_stokesU = malloc(n_powers*sizeof(user_precision_t));
  components->power_ref_stokesV = malloc(n_powers*sizeof(user_precision_t));
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
  free( components.power_ref_stokesQ );
  free( components.power_ref_stokesU );
  free( components.power_ref_stokesV );
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
                  double * expec_re, double * expec_im) {
  double expec_re_inc, expec_im_inc;
  double temp; //, visi_freq;
  double flux_extrap, xx_gain; //flux_ratio, ,
  int time_ind, freq_ind, beam_ind;
  * expec_re = 0.0;
  * expec_im = 0.0;

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
    flux_extrap = expec_flux_I[num_freqs*comp + freq_ind];

    // printf("%.5f %.5f %.5f\n",components.power_ref_freqs[comp],
    //                           components.power_SIs[comp],
    //                           components.power_ref_stokesI[comp] );

    if (beamtype == NO_BEAM) {
      xx_gain = 1.0;
    } else {
      //In tests below have set gains purely real, and all pols equal, so
      //only need to call one pol and square it
      xx_gain = creal(args_ft->primay_beam_J00[beam_ind]);
      xx_gain = xx_gain*xx_gain;
    }

    expec_re_inc = expec_re_inc*flux_extrap*xx_gain;
    expec_im_inc = expec_im_inc*flux_extrap*xx_gain;

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
    * expec_re += (expec_re_inc*creal(visi_env) - expec_im_inc*cimag(visi_env));
    * expec_im += (expec_re_inc*cimag(visi_env) + expec_im_inc*creal(visi_env));
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
  if (component_type == POINT) {
    #ifdef DOUBLE_PRECISION
      frac_tol = 1e-13;
    #else
      frac_tol = 7e-5;
    #endif
  }
  else if (component_type == GAUSSIAN) {
    #ifdef DOUBLE_PRECISION
      frac_tol = 1e-12;
    #else
      frac_tol = 7e-5;
    #endif
  }
  else {
    #ifdef DOUBLE_PRECISION
      frac_tol = 1e-12;
    #else
      frac_tol = 5e-3;
    #endif
  }


  double expec_re, expec_im;
  for (int visi = 0; visi < num_visis; visi++) {

      get_expected(visi, num_powers, num_curves, num_lists,
                  num_baselines, num_freqs, extrap_freqs,
                  beamtype, args_ft, components, component_type,
                  &expec_re, &expec_im);

    // printf("%.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
    //       args_ft->sum_visi_XX_real[visi], args_ft->sum_visi_XX_imag[visi],
    //       args_ft->sum_visi_XY_real[visi], args_ft->sum_visi_XY_imag[visi],
    //       args_ft->sum_visi_YX_real[visi], args_ft->sum_visi_YX_imag[visi],
    //       args_ft->sum_visi_YY_real[visi], args_ft->sum_visi_YY_imag[visi]);
    //
    if (beamtype == FEE_BEAM || beamtype == FEE_BEAM_INTERP || beamtype == MWA_ANALY ) {
      //MWA beam has cross pols which double everything when using just Stokes I
      //and setting the cross pols to 1.0 as well as the gains
      //Also means cross-pols are non-zero

      // printf("%.16f %.16f\n",2*expec_re, args_ft->sum_visi_XX_real[visi]);

      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_re, 2.0*expec_re, args_ft->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(2*frac_tol*expec_im, 2.0*expec_im, args_ft->sum_visi_YY_imag[visi]);
    }
    else {

      // printf("%.8f %.8f\n", expec_re, args_ft->sum_visi_XX_real[visi]);
      // printf("%.8f %.8f\n", expec_im, args_ft->sum_visi_XX_imag[visi]);

      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re, expec_re, args_ft->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im, expec_im, args_ft->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol, 0.0, args_ft->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re, expec_re, args_ft->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im, expec_im, args_ft->sum_visi_YY_imag[visi]);

    }
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
  for (size_t visi = 0; visi < num_beam_values; visi++) {
    args_ft->primay_beam_J00[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    components.power_ref_stokesI[comp] = 1.0;
    components.power_ref_stokesQ[comp] = 0.0;
    components.power_ref_stokesU[comp] = 0.0;
    components.power_ref_stokesV[comp] = 0.0;
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
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  for (int pow_ind = 0; pow_ind < n_powers; pow_ind++) {
    components.power_ref_stokesI[pow_ind] = (pow_ind + 1);
    components.power_ref_stokesQ[pow_ind] = 0.0;
    components.power_ref_stokesU[pow_ind] = 0.0;
    components.power_ref_stokesV[pow_ind] = 0.0;
    components.power_SIs[pow_ind] = -0.1*pow_ind;
    components.power_ref_freqs[pow_ind] = 150e+6;
    components.power_comp_inds[pow_ind] = pow_ind;
  }

  // //Just have some premade arrays living in the header and just stick em in
  // //for the CURVED_POWER_LAW and LIST components
  //
  components.curve_ref_stokesI = k_ref_stokesI;
  components.curve_ref_stokesQ = k_ref_stokesQ;
  components.curve_ref_stokesU = k_ref_stokesU;
  components.curve_ref_stokesV = k_ref_stokesV;
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
    args_ft->primay_beam_J01[beam] = beam + I*0.0;
    args_ft->primay_beam_J10[beam] = beam + I*0.0;
    args_ft->primay_beam_J11[beam] = beam + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    components.power_ref_stokesI[comp] = 1.0;
    components.power_ref_stokesQ[comp] = 0.0;
    components.power_ref_stokesU[comp] = 0.0;
    components.power_ref_stokesV[comp] = 0.0;
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
    args_ft->primay_beam_J01[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J10[visi] = 1.0 + I*0.0;
    args_ft->primay_beam_J11[visi] = 1.0 + I*0.0;
  }

  //Just stick Stokes I to 1.0, SI to zero, and reference freqs to 150MHz
  for (int comp = 0; comp < num_components; comp++) {
    components.power_ref_stokesI[comp] = 1.0;
    components.power_ref_stokesQ[comp] = 0.0;
    components.power_ref_stokesU[comp] = 0.0;
    components.power_ref_stokesV[comp] = 0.0;
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
