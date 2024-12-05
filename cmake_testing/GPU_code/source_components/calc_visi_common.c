#include "calc_visi_common.h"

int k_n_stokesV_pol_frac = 3;
int k_n_stokesV_power = 4;
int k_n_stokesV_curve = 4;
int k_n_stokesV_list = 3;
int k_n_linpol_pol_frac = 5;
int k_n_linpol_power = 2;
int k_n_linpol_curve = 2;
int k_n_linpol_list = 2;
int k_n_linpol_p_list = 3;
int k_n_stokesV_list_flux_entries = 56;
int k_n_stokesQ_list_flux_entries = 20;
int k_n_stokesU_list_flux_entries = 21;
int k_n_linpol_p_list_flux_entries = 60;
user_precision_t k_ref_stokesI[] = {65.8709372597089384,54.7884483174221160,97.5061615843035980,61.6594428464970150,99.5297177363038230,66.9357417400148904,45.0212002120857093,35.7359617314428277,5.9735382808915194,92.2230329259136141};

double k_ref_freqs[] = {5.0000000000000000e+07,7.7777777777777776e+07,1.0555555555555555e+08,1.3333333333333333e+08,1.6111111111111110e+08,1.8888888888888890e+08,2.1666666666666666e+08,2.4444444444444442e+08,2.7222222222222221e+08,3.0000000000000000e+08};

user_precision_t k_ref_curve_SIs[] = {0.7834417047479292,-0.4855866984455952,0.7190725095506663,-0.6764787207262060,0.8770054913828549,-0.5122915335891907,-0.3734335809307321,0.6013225788437182,0.3992912015429384,0.3594219043826412};

user_precision_t k_ref_qs[] = {-0.4764692401215893,-1.8399065653084528,-0.0811322256939060,-1.2718847245617342,-0.6345384442984319,-1.8177562920236174,-1.2563413339395475,-1.2927411814359573,-0.7323838368949018,-0.8809858973081504};

int k_num_list_values[] = {22,12,21,31,31};

int k_list_start_indexes[] = {0,22,34,55,86};

double k_list_freqs[] = {5.0000000000000000e+07,6.1904761904761903e+07,7.3809523809523806e+07,8.5714285714285716e+07,9.7619047619047612e+07,1.0952380952380952e+08,1.2142857142857143e+08,1.3333333333333333e+08,1.4523809523809522e+08,1.5714285714285713e+08,1.6904761904761904e+08,1.8095238095238096e+08,1.9285714285714287e+08,2.0476190476190478e+08,2.1666666666666666e+08,2.2857142857142857e+08,2.4047619047619048e+08,2.5238095238095239e+08,2.6428571428571430e+08,2.7619047619047618e+08,2.8809523809523809e+08,3.0000000000000000e+08,5.0000000000000000e+07,7.2727272727272719e+07,9.5454545454545453e+07,1.1818181818181819e+08,1.4090909090909091e+08,1.6363636363636363e+08,1.8636363636363637e+08,2.0909090909090909e+08,2.3181818181818181e+08,2.5454545454545453e+08,2.7727272727272725e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.2500000000000000e+07,7.5000000000000000e+07,8.7500000000000000e+07,1.0000000000000000e+08,1.1250000000000000e+08,1.2500000000000000e+08,1.3750000000000000e+08,1.5000000000000000e+08,1.6250000000000000e+08,1.7500000000000000e+08,1.8750000000000000e+08,2.0000000000000000e+08,2.1250000000000000e+08,2.2500000000000000e+08,2.3750000000000000e+08,2.5000000000000000e+08,2.6250000000000000e+08,2.7500000000000000e+08,2.8750000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,5.8333333333333336e+07,6.6666666666666664e+07,7.5000000000000000e+07,8.3333333333333328e+07,9.1666666666666657e+07,1.0000000000000000e+08,1.0833333333333333e+08,1.1666666666666666e+08,1.2500000000000000e+08,1.3333333333333333e+08,1.4166666666666666e+08,1.5000000000000000e+08,1.5833333333333331e+08,1.6666666666666666e+08,1.7500000000000000e+08,1.8333333333333331e+08,1.9166666666666666e+08,2.0000000000000000e+08,2.0833333333333331e+08,2.1666666666666666e+08,2.2500000000000000e+08,2.3333333333333331e+08,2.4166666666666666e+08,2.5000000000000000e+08,2.5833333333333331e+08,2.6666666666666666e+08,2.7500000000000000e+08,2.8333333333333331e+08,2.9166666666666663e+08,3.0000000000000000e+08,5.0000000000000000e+07,5.8333333333333336e+07,6.6666666666666664e+07,7.5000000000000000e+07,8.3333333333333328e+07,9.1666666666666657e+07,1.0000000000000000e+08,1.0833333333333333e+08,1.1666666666666666e+08,1.2500000000000000e+08,1.3333333333333333e+08,1.4166666666666666e+08,1.5000000000000000e+08,1.5833333333333331e+08,1.6666666666666666e+08,1.7500000000000000e+08,1.8333333333333331e+08,1.9166666666666666e+08,2.0000000000000000e+08,2.0833333333333331e+08,2.1666666666666666e+08,2.2500000000000000e+08,2.3333333333333331e+08,2.4166666666666666e+08,2.5000000000000000e+08,2.5833333333333331e+08,2.6666666666666666e+08,2.7500000000000000e+08,2.8333333333333331e+08,2.9166666666666663e+08,3.0000000000000000e+08};

user_precision_t k_list_stokesI[] = {0.5919678991675268,1.8228498423748145,5.9405365939697949,1.5625445665470608,4.4084460689458762,9.8435358082632156,4.9391519986156496,-0.6229649227834282,3.5975357440649454,5.3150343885848068,8.5044494441924634,6.0035924620933638,-0.0486189891697234,8.4217583807665193,2.9345973906064202,2.8259859398945721,2.0913063201663831,2.9394929407690848,9.8642299124563593,4.9723647991259288,1.8884774633685799,8.8789583175181956,1.7581768636464705,3.7799313203528921,3.3510132497117562,3.9920751513153778,3.9186088981374105,3.4235566519712624,4.8840585143225299,6.7348361972361381,9.5598614078282385,-0.9858182988499411,9.5639426627951121,9.8213519113121048,2.8806968470309089,4.4422993947468177,4.6378468238054085,3.7880187150402245,4.2516627802273614,6.3677485706884172,8.3020559665977469,2.2557806568867744,-0.8304521052162194,0.1710632534141872,8.6833585789401511,6.7672554816642307,4.4779905434678069,8.2977661851982596,4.8179846938306996,4.4967630484526673,2.4618779631114549,8.0196319264493390,7.5520746109939338,4.0339870155196458,1.0308084101113169,2.0151686874541159,0.0374451142713135,4.3151200178405507,1.6506266741365225,4.3707111028715255,6.6375615020087153,7.1950507965954653,7.4209426862366978,3.0601838336018625,4.8373567316396366,9.9545783987905878,9.4435861043013922,4.4661004669567630,-0.6493501294418613,8.6396843257040636,4.5214413663464654,4.6100144706988031,8.5385740413910529,9.8823275410754441,5.4926166984334426,4.2754452094640598,1.8026930926116811,8.2082771834151771,5.5518712358007027,7.5878071276961148,3.1051310317432312,4.5092688129763632,6.4469260171724265,0.0289984191454642,9.4607413270708776,3.5358824632678374,-0.6102275300014381,2.4399380852001569,2.7512173712287535,8.1224619377481062,-0.2850801896786408,3.1305286232301173,0.8779882606189866,6.7603753315074382,8.0695190741068910,5.1134736252840343,9.9208465403723274,2.9287630149663597,5.8191320751887883,5.9264801562503600,7.6904162365624966,9.9147061056543091,4.1115201815822449,8.7307712481894644,1.1606179093239182,3.7808186072741403,4.4534097644250963,0.9304279423448416,3.4083382810881693,6.5776090074052904,2.1824403782051061,5.8583499625625182,-0.2806297866981025,4.5987069652623749,4.7168754879015822,1.4755864701835053,4.7098082807398782};

int k_curve_comp_inds[] = {10,11,12,13,14,15,16,17,18,19};

int k_list_comp_inds[] = {20,21,22,23,24};

user_precision_t k_stokesV_pol_fracs[] = {0.3727880063483699,0.2722657235369972,0.4536716460387014};

user_precision_t k_ref_stokesV[] = {4.5237369205498208,0.1211208591499693,4.3482893917986116,5.6126003078959874};

user_precision_t k_stokesV_power_SIs[] = {-0.2308043049018835,-0.0800271065470730,0.1387604933787128,-0.8609710769138863};

user_precision_t k_stokesV_qs[] = {-1.0171064164267634,-1.4649381924994342,-1.1992770450377468,-1.9726498902802543};

user_precision_t k_stokesV_curve_SIs[] = {0.2358907864290547,-0.8774455219171808,-0.1169270840449850,-0.2301849313278712};

user_precision_t k_linpol_pol_fracs[] = {0.1639618312305009,-0.7808049675850461,-0.4915756058257663,0.3006811790785389,-0.5046672154237468};

user_precision_t k_ref_linpol[] = {2.7784514381487733,8.3883664824101434};

user_precision_t k_linpol_power_SIs[] = {-0.4193819205566838,-1.0701268561606823};

user_precision_t k_linpol_qs[] = {0.4707351104750854,-1.0942856432833969};

user_precision_t k_linpol_curve_SIs[] = {-0.3915238172990896,0.7816921730430624};

user_precision_t k_intr_pol_angle[] = {1.0700001852292080,4.1596409972184469,0.8576590829591255,3.2958018165048886,2.8089690585343687,5.0975565519404622,2.7264431644107967,4.1795871066003567,4.0302374403335612,0.0272017852272141,5.5381500238193846,0.1820141161926313};

user_precision_t k_rms[] = {26.2642581283263787,32.4051819659594713,59.5627617193057191,21.8789629791023188,32.8222127609839660,73.1691791872270727,77.9414899000087900,40.5007682988345010,13.5857831818523422,52.7269561282646890,0.9821578083106353,46.2036940423003273};

int k_stokesV_pol_frac_comp_inds[] = {6,18,12};

int k_stokesV_power_comp_inds[] = {24,5,3,7};

int k_stokesV_curve_comp_inds[] = {8,0,20,9};

int k_linpol_pol_frac_comp_inds[] = {18,19,2,8,5};

int k_linpol_power_comp_inds[] = {11,24};

int k_linpol_curve_comp_inds[] = {0,9};

int k_linpol_angle_inds[] = {18,19,2,8,5,11,24,0,9,13,15,16};

user_precision_t k_stokesV_list_ref_flux[] = {8.4815709208509364,3.5641547764196551,9.5152457575955367,4.2448576137525267,1.4717749387433523,4.7181635099461490,6.2769316220933566,0.1957631796900456,2.3885233088179061,2.1090957833421022,7.1894055875333542,2.4023026083484451,1.4772078135115403,8.0425923476370507,8.6030937420875109,9.3806627338474708,4.1155911234713338,-0.5941362102501802,8.9688634036075783,-0.0065137488875423,8.8875701259348787,-0.3035087226305109,-0.9338549285210060,2.7221214243960961,8.2609958023771473,1.7417150170027194,0.2556007027068437,5.2586243737818972,9.5448551640616817,2.7528896601974298,1.5035962485431287,-0.7324896825336026,8.4214014315307608,1.8861935328338806,9.9993871108549346,7.1218714964533216,5.8992373859151943,3.2475160962399006,9.0058292928296364,-0.3017802739768610,4.2774383883223051,7.5865813324664302,6.8837394651481238,7.3229829090642973,1.3352976759151858,2.5873047898522343,2.8332896661347098,5.7857895496383884,8.7571004080360009,4.0850143136841108,7.6809023546508417,5.1457544673365057,4.8702122386440072,7.2712126163437709,9.8447667337486315,4.9651327314606881};

double k_stokesV_list_ref_freqs[] = {5.0000000000000000e+07,6.4705882352941178e+07,7.9411764705882356e+07,9.4117647058823526e+07,1.0882352941176471e+08,1.2352941176470588e+08,1.3823529411764705e+08,1.5294117647058824e+08,1.6764705882352942e+08,1.8235294117647058e+08,1.9705882352941176e+08,2.1176470588235295e+08,2.2647058823529410e+08,2.4117647058823529e+08,2.5588235294117647e+08,2.7058823529411763e+08,2.8529411764705884e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.1363636363636360e+07,7.2727272727272719e+07,8.4090909090909094e+07,9.5454545454545453e+07,1.0681818181818181e+08,1.1818181818181819e+08,1.2954545454545455e+08,1.4090909090909091e+08,1.5227272727272725e+08,1.6363636363636363e+08,1.7500000000000000e+08,1.8636363636363637e+08,1.9772727272727272e+08,2.0909090909090909e+08,2.2045454545454544e+08,2.3181818181818181e+08,2.4318181818181819e+08,2.5454545454545453e+08,2.6590909090909091e+08,2.7727272727272725e+08,2.8863636363636363e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.7857142857142866e+07,8.5714285714285716e+07,1.0357142857142857e+08,1.2142857142857143e+08,1.3928571428571430e+08,1.5714285714285713e+08,1.7500000000000000e+08,1.9285714285714287e+08,2.1071428571428573e+08,2.2857142857142860e+08,2.4642857142857143e+08,2.6428571428571430e+08,2.8214285714285719e+08,3.0000000000000000e+08};

int k_stokesV_num_list_values[] = {18,23,15};

int k_stokesV_list_start_indexes[] = {0,18,41};

int k_stokesV_list_comp_inds[] = {23,14,19};

user_precision_t k_stokesQ_list_ref_flux[] = {5.8194199397274318,7.0428476593734803,9.1937946790702547,-0.4355619806412979,6.9323228879429264,9.2542128509020110,1.4111839299335562,-0.7483475987338879,1.2796103567034165,6.6569861559533532,9.9760677306309287,6.4004852753837653,5.3640786580085882,9.8315043112271781,6.7465896366656031,8.6620778428596097,6.4174734365375246,5.0707910047713458,7.4508080961542174,4.7054100239358991};

double k_stokesQ_list_ref_freqs[] = {5.0000000000000000e+07,1.0000000000000000e+08,1.5000000000000000e+08,2.0000000000000000e+08,2.5000000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.9230769230769232e+07,8.8461538461538464e+07,1.0769230769230770e+08,1.2692307692307693e+08,1.4615384615384614e+08,1.6538461538461539e+08,1.8461538461538464e+08,2.0384615384615386e+08,2.2307692307692307e+08,2.4230769230769232e+08,2.6153846153846157e+08,2.8076923076923078e+08,3.0000000000000000e+08};

int k_stokesQ_num_list_values[] = {6,14};

int k_stokesQ_list_start_indexes[] = {0,6};

int k_stokesQ_list_comp_inds[] = {6,17};

user_precision_t k_stokesU_list_ref_flux[] = {9.3146864418986368,9.6635846792127609,-0.9017777596629057,0.6376766098326492,9.8852477970008970,-0.2894888330189717,3.7875523480491147,6.0250582516722941,8.8687190348966354,5.6900227187552641,0.6258232628761622,7.1483926309234036,9.8581391559871641,8.2258786397012393,-0.7111866374449451,8.0247290216725080,1.5912514809750316,6.6749949155066872,9.4959719421033562,4.0685348272837487,0.7299894724561886};

double k_stokesU_list_ref_freqs[] = {5.0000000000000000e+07,7.5000000000000000e+07,1.0000000000000000e+08,1.2500000000000000e+08,1.5000000000000000e+08,1.7500000000000000e+08,2.0000000000000000e+08,2.2500000000000000e+08,2.5000000000000000e+08,2.7500000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,7.7777777777777776e+07,1.0555555555555555e+08,1.3333333333333333e+08,1.6111111111111110e+08,1.8888888888888890e+08,2.1666666666666666e+08,2.4444444444444442e+08,2.7222222222222221e+08,3.0000000000000000e+08};

int k_stokesU_num_list_values[] = {11,10};

int k_stokesU_list_start_indexes[] = {0,11};

int k_stokesU_list_comp_inds[] = {6,17};

user_precision_t k_linpol_p_list_ref_flux[] = {0.6896785847245828,0.8345757291325908,9.8893391654864526,4.4368804994205391,-0.0236045318815440,8.3927380891490522,4.8261284349987434,7.1123839330083083,5.3651157374468941,2.8024558896473994,9.8111688798894168,9.4786521481399664,9.6020595446228132,5.8282598595021220,4.2113500280615117,-0.2822592960869456,6.3962240140124091,6.8853698294932757,6.4851417479246960,2.1543882381554584,8.3141148973802270,5.4822784913322051,2.1984463981734965,-0.2730193774003482,7.4025062269995541,6.6110485321542818,4.1380035253406495,6.9108286075099823,5.5687587416613198,8.6713290825318250,4.8869660603857445,1.2099295899417020,5.1458676670739454,7.5285705389528772,8.6699326074478194,4.0161898493389838,8.1600090095467017,0.7947157851768085,6.4543288863408570,8.3363491982909572,2.0897973508446599,1.3557427821430368,7.0775048837908177,9.1559759369935261,1.9040569339253626,4.0881247560428671,6.9825988278315849,4.8339758812148430,4.7388076963443826,0.6458809200418552,5.7087451055548009,7.5958855056004051,2.7170930097857324,4.3727976757461553,0.9786174372456617,5.1633225917555663,7.6566209398541378,6.5043301280542405,7.9668840831400338,9.8772734205625312};

double k_linpol_p_list_ref_freqs[] = {5.0000000000000000e+07,6.3157894736842103e+07,7.6315789473684207e+07,8.9473684210526317e+07,1.0263157894736841e+08,1.1578947368421052e+08,1.2894736842105263e+08,1.4210526315789473e+08,1.5526315789473683e+08,1.6842105263157895e+08,1.8157894736842105e+08,1.9473684210526314e+08,2.0789473684210527e+08,2.2105263157894737e+08,2.3421052631578946e+08,2.4736842105263159e+08,2.6052631578947368e+08,2.7368421052631581e+08,2.8684210526315790e+08,3.0000000000000000e+08,5.0000000000000000e+07,7.7777777777777776e+07,1.0555555555555555e+08,1.3333333333333333e+08,1.6111111111111110e+08,1.8888888888888890e+08,2.1666666666666666e+08,2.4444444444444442e+08,2.7222222222222221e+08,3.0000000000000000e+08,5.0000000000000000e+07,5.8620689655172415e+07,6.7241379310344830e+07,7.5862068965517238e+07,8.4482758620689660e+07,9.3103448275862068e+07,1.0172413793103448e+08,1.1034482758620688e+08,1.1896551724137931e+08,1.2758620689655171e+08,1.3620689655172414e+08,1.4482758620689654e+08,1.5344827586206895e+08,1.6206896551724136e+08,1.7068965517241377e+08,1.7931034482758620e+08,1.8793103448275861e+08,1.9655172413793102e+08,2.0517241379310343e+08,2.1379310344827586e+08,2.2241379310344827e+08,2.3103448275862068e+08,2.3965517241379309e+08,2.4827586206896549e+08,2.5689655172413790e+08,2.6551724137931034e+08,2.7413793103448272e+08,2.8275862068965518e+08,2.9137931034482753e+08,3.0000000000000000e+08};

int k_linpol_p_num_list_values[] = {20,10,30};

int k_linpol_p_list_start_indexes[] = {0,20,30};

int k_linpol_p_list_comp_inds[] = {13,15,16};

int total_num_flux_entires = 117;






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
  components->n_stokesV_list = 0;
  components->n_linpol_list = 0;
  components->n_linpol_p_list = 0;

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

  // for (int i = 0; i < num_components; i++)
  // {
  //   printf("comp %d: %.5f %.5f %.5f %.5f\n", i, expec_flux_I[num_freqs*i], expec_flux_Q[num_freqs*i], expec_flux_U[num_freqs*i], expec_flux_V[num_freqs*i]);
  // }
  

  for (int comp = 0; comp < num_components; comp++) {
    temp = 2*M_PI*( args_ft->us[visi]*components.ls[comp] + args_ft->vs[visi]*components.ms[comp] + args_ft->ws[visi]*(components.ns[comp]-1) );
    // sincos(temp, &(expec_im_inc), &(expec_re_inc));

    expec_im_inc = sin(temp);
    expec_re_inc = cos(temp);

    //Do some indexing so we can call up correct instrumental gain
    time_ind = (int)floorf( (float)visi / ((float)num_baselines * (float)num_freqs));
    freq_ind = (int)floorf( ((float)visi - ((float)time_ind*(float)num_baselines * (float)num_freqs)) / (float)num_baselines);
    beam_ind = num_freqs*time_ind*num_components + (num_components*freq_ind) + comp;

    flux_I = expec_flux_I[num_freqs*comp + freq_ind];
    flux_Q = expec_flux_Q[num_freqs*comp + freq_ind];
    flux_U = expec_flux_U[num_freqs*comp + freq_ind];
    flux_V = expec_flux_V[num_freqs*comp + freq_ind];

    if (flux_I > 1000 || flux_Q > 1000 || flux_U > 1000 || flux_V > 1000) {
      printf("What have we got %d %f %f %f %f\n", comp, flux_I, flux_Q, flux_U, flux_V);
      
    }

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
  double abs_tol;
  //For some reason, AMD doesn't seem as accurate
  if (component_type == POINT) {
    #ifdef DOUBLE_PRECISION
      #ifdef __HIPCC__
        frac_tol = 1e-10;
        abs_tol = 1e-10;
      #else
        frac_tol = 1e-12;
        abs_tol = 1e-11;
      #endif
    #else
      frac_tol = 7e-5;
      abs_tol = 7e-5;
    #endif
  }
  else if (component_type == GAUSSIAN) {
    #ifdef DOUBLE_PRECISION
      #ifdef __HIPCC__
        frac_tol = 1e-10;
        abs_tol = 1e-10;
      #else
        frac_tol = 1e-11;
        abs_tol = 1e-11;
      #endif
    #else
      frac_tol = 7e-5;
      abs_tol = 7e-5;
    #endif
  }
  else {
    #ifdef DOUBLE_PRECISION
      #ifdef __HIPCC__
        frac_tol = 1e-10;
        abs_tol = 1e-10;
      #else
        frac_tol = 1e-11;
        abs_tol = 1e-11;
      #endif
    #else
      frac_tol = 9e-5;
      abs_tol = 7e-5;
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

      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_xx, expec_re_xx, args_ft->sum_visi_XX_real[visi]);
      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_xx, expec_im_xx, args_ft->sum_visi_XX_imag[visi]);
      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_xy, expec_re_xy, args_ft->sum_visi_XY_real[visi]);
      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_xy, expec_im_xy, args_ft->sum_visi_XY_imag[visi]);
      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_yx, expec_re_yx, args_ft->sum_visi_YX_real[visi]);
      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_yx, expec_im_yx, args_ft->sum_visi_YX_imag[visi]);
      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_yy, expec_re_yy, args_ft->sum_visi_YY_real[visi]);
      // TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_yy, expec_im_yy, args_ft->sum_visi_YY_imag[visi]);

      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_xx + abs_tol, expec_re_xx, args_ft->sum_visi_XX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_xx + abs_tol, expec_im_xx, args_ft->sum_visi_XX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_xy + abs_tol, expec_re_xy, args_ft->sum_visi_XY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_xy + abs_tol, expec_im_xy, args_ft->sum_visi_XY_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_yx + abs_tol, expec_re_yx, args_ft->sum_visi_YX_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_yx + abs_tol, expec_im_yx, args_ft->sum_visi_YX_imag[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_re_yy + abs_tol, expec_re_yy, args_ft->sum_visi_YY_real[visi]);
      TEST_ASSERT_DOUBLE_WITHIN(frac_tol*expec_im_yy + abs_tol, expec_im_yy, args_ft->sum_visi_YY_imag[visi]);

  }
}


void do_gpu_calc_visi(int n_powers, int n_curves, int n_lists,
          int num_baselines, int num_shape_coeffs,
          int num_freqs, int num_cross, int num_times,
          e_beamtype beamtype,  e_component_type comptype,
          components_t *components,
          source_t *chunked_source, 
          args_for_testing_t *args_ft){

  woden_settings_t *woden_settings = malloc(sizeof(woden_settings_t));
  woden_settings->do_autos = 0;
  woden_settings->num_freqs = num_freqs;
  woden_settings->num_time_steps = num_times;
  woden_settings->num_baselines = num_baselines;
  woden_settings->num_cross = num_cross;
  woden_settings->use_dipamps = 0;

  // d_beam_gains_t *d_beam_gains = malloc(sizeof(d_beam_gains_t));
  // visibility_set_t *d_visibility_set = NULL;
  

  source_t *d_chunked_source = copy_chunked_source_to_GPU(chunked_source);
  double *d_freqs = malloc_freqs_gpu(num_freqs, args_ft->extrap_freqs);

  components_t *d_components;
  int num_components;

  if (comptype == POINT) {
    d_components = &d_chunked_source->point_components;
    num_components = d_chunked_source->n_points;
  }
  else if (comptype == GAUSSIAN) {
    d_components = &d_chunked_source->gauss_components;
    num_components = d_chunked_source->n_gauss;
  }
  else {
    d_components = &d_chunked_source->shape_components;
    num_components = d_chunked_source->n_shapes;
  }

  malloc_extrapolated_flux_arrays_gpu(d_components, num_components, num_freqs);
  int do_gpu = 1;
  extrapolate_Stokes(d_chunked_source, d_freqs, num_freqs, comptype, do_gpu);

  test_kern_calc_visi_all(n_powers, n_curves, n_lists, num_baselines, num_shape_coeffs,
            num_freqs, num_cross, num_times, beamtype, comptype, components,
            d_chunked_source, d_freqs,
            args_ft->us, args_ft->vs, args_ft->ws,
            args_ft->u_shapes, args_ft->v_shapes,
            args_ft->sum_visi_XX_real, args_ft->sum_visi_XX_imag,
            args_ft->sum_visi_XY_real, args_ft->sum_visi_XY_imag,
            args_ft->sum_visi_YX_real, args_ft->sum_visi_YX_imag,
            args_ft->sum_visi_YY_real, args_ft->sum_visi_YY_imag,
            args_ft->allsteps_wavelengths, args_ft->sbf,
            args_ft->primay_beam_J00, args_ft->primay_beam_J01,
            args_ft->primay_beam_J10, args_ft->primay_beam_J11);

  free_d_components(d_chunked_source, comptype);
  free_extrapolated_flux_arrays(d_components);
  free_freqs_gpu(d_freqs);

}

// source_t * put_components_into_source(components_t components,
//                                       e_component_type comptype,
//                                       int n_powers, int n_curves, int n_lists,
//                                       int num_shape_coeffs) {

//   source_t *chunked_source = (source_t *)malloc(sizeof(source_t));
//   if (comptype == POINT) {

//     chunked_source->point_components = components;
//     chunked_source->n_points = n_powers + n_curves + n_lists;
//     chunked_source->n_point_powers = n_powers;
//     chunked_source->n_point_curves = n_curves;
//     chunked_source->n_point_lists = n_lists;

//     chunked_source->n_gauss = 0;
//     chunked_source->n_gauss_lists = 0;
//     chunked_source->n_gauss_powers = 0;
//     chunked_source->n_gauss_curves = 0;
//     chunked_source->n_shapes = 0;
//     chunked_source->n_shape_lists = 0;
//     chunked_source->n_shape_powers = 0;
//     chunked_source->n_shape_curves = 0;
//     chunked_source->n_shape_coeffs = 0;

//   }
//   else if (comptype == GAUSSIAN) {

//     chunked_source->gauss_components = components;
//     chunked_source->n_gauss = n_powers + n_curves + n_lists;
//     chunked_source->n_gauss_powers = n_powers;
//     chunked_source->n_gauss_curves = n_curves;
//     chunked_source->n_gauss_lists = n_lists;

//     chunked_source->n_points = 0;
//     chunked_source->n_point_lists = 0;
//     chunked_source->n_point_powers = 0;
//     chunked_source->n_point_curves = 0;
//     chunked_source->n_shapes = 0;
//     chunked_source->n_shape_lists = 0;
//     chunked_source->n_shape_powers = 0;
//     chunked_source->n_shape_curves = 0;
//     chunked_source->n_shape_coeffs = 0;
//   }
//   else if (comptype == SHAPELET) {
//     chunked_source->shape_components = components;
//     chunked_source->n_shapes = n_powers + n_curves + n_lists;
//     chunked_source->n_shape_powers = n_powers;
//     chunked_source->n_shape_curves = n_curves;
//     chunked_source->n_shape_lists = n_lists;
//     chunked_source->n_shape_coeffs = num_shape_coeffs;

//     chunked_source->n_points = 0;
//     chunked_source->n_point_lists = 0;
//     chunked_source->n_point_powers = 0;
//     chunked_source->n_point_curves = 0;
//     chunked_source->n_gauss = 0;
//     chunked_source->n_gauss_lists = 0;
//     chunked_source->n_gauss_powers = 0;
//     chunked_source->n_gauss_curves = 0;
//   }

//   return chunked_source;
// }

/*
Test the __global__ code that calculates visibilities for different type
of COMPONENTs
Vary the l,m,n coords but keep all other variables constant
*/
void test_calc_visi_Varylmn(e_beamtype beamtype, e_component_type comptype,
                            int do_gpu) {

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

  source_t *chunked_source = put_components_into_source(components, comptype,
                                                        n_powers, n_curves,
                                                        n_lists, num_coeffs);

  if (do_gpu == 1){

    do_gpu_calc_visi(n_powers, n_curves, n_lists,
          num_baselines, num_coeffs,
          num_freqs, num_visis, num_times,
          beamtype, comptype, &components,
          chunked_source,
          args_ft);

  }

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
void test_calc_visi_VarylmnVaryFlux(e_beamtype beamtype, e_component_type comptype,
                                    int do_gpu) {

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
    components.power_SIs[pow_ind] = -0.1*pow_ind;
    components.power_ref_freqs[pow_ind] = 150e+6;
    components.power_comp_inds[pow_ind] = pow_ind;
  }

  // //Just have some premade arrays living in the header and just stick em in
  // //for the CURVED_POWER_LAW and LIST components
  //
  components.curve_ref_stokesI = k_ref_stokesI;
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


  components.n_stokesV_list_flux_entries = k_n_stokesV_list_flux_entries;
  components.n_stokesQ_list_flux_entries = k_n_stokesQ_list_flux_entries;
  components.n_stokesU_list_flux_entries = k_n_stokesU_list_flux_entries;
  components.n_linpol_p_list_flux_entries = k_n_linpol_p_list_flux_entries;
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
  components.n_stokesV_pol_frac = k_n_stokesV_pol_frac;
  components.stokesV_pol_fracs = k_stokesV_pol_fracs;
  components.stokesV_pol_frac_comp_inds = k_stokesV_pol_frac_comp_inds;

  components.n_stokesV_power = k_n_stokesV_power;
  components.stokesV_power_ref_flux = k_ref_stokesV;
  components.stokesV_power_SIs = k_stokesV_power_SIs;
  components.stokesV_power_comp_inds =k_stokesV_power_comp_inds;
  
  components.n_stokesV_curve = k_n_stokesV_curve;
  components.stokesV_curve_ref_flux = k_ref_stokesV;
  components.stokesV_curve_SIs = k_stokesV_curve_SIs;
  components.stokesV_curve_qs = k_stokesV_qs;
  components.stokesV_curve_comp_inds = k_stokesV_curve_comp_inds;

  components.n_linpol_pol_frac = k_n_linpol_pol_frac;
  components.linpol_pol_fracs = k_linpol_pol_fracs;
  components.linpol_pol_frac_comp_inds = k_linpol_pol_frac_comp_inds;

  components.n_linpol_power = k_n_linpol_power;
  components.linpol_power_ref_flux = k_ref_linpol;
  components.linpol_power_SIs = k_linpol_power_SIs;
  components.linpol_power_comp_inds = k_linpol_power_comp_inds;
  
  components.n_linpol_curve = k_n_linpol_curve;
  components.linpol_curve_ref_flux = k_ref_linpol;
  components.linpol_curve_SIs = k_linpol_curve_SIs;
  components.linpol_curve_qs = k_linpol_qs;
  components.linpol_curve_comp_inds = k_linpol_curve_comp_inds;

  components.n_stokesV_list = k_n_stokesV_list;
  components.stokesV_list_ref_flux = k_stokesV_list_ref_flux;
  components.stokesV_list_ref_freqs = k_stokesV_list_ref_freqs;
  components.stokesV_num_list_values = k_stokesV_num_list_values;
  components.stokesV_list_start_indexes = k_stokesV_list_start_indexes;
  components.stokesV_list_comp_inds = k_stokesV_list_comp_inds;

  components.n_linpol_list = k_n_linpol_list;
  components.stokesQ_list_ref_flux = k_stokesQ_list_ref_flux;
  components.stokesQ_list_ref_freqs = k_stokesQ_list_ref_freqs;
  components.stokesQ_num_list_values = k_stokesQ_num_list_values;
  components.stokesQ_list_start_indexes = k_stokesQ_list_start_indexes;
  components.stokesQ_list_comp_inds = k_stokesQ_list_comp_inds;
  components.stokesU_list_ref_flux = k_stokesU_list_ref_flux;
  components.stokesU_list_ref_freqs = k_stokesU_list_ref_freqs;
  components.stokesU_num_list_values = k_stokesU_num_list_values;
  components.stokesU_list_start_indexes = k_stokesU_list_start_indexes;
  components.stokesU_list_comp_inds = k_stokesU_list_comp_inds;

  components.n_linpol_p_list = k_n_linpol_p_list;
  components.linpol_p_list_ref_flux = k_linpol_p_list_ref_flux;
  components.linpol_p_list_ref_freqs = k_linpol_p_list_ref_freqs;
  components.linpol_p_num_list_values = k_linpol_p_num_list_values;
  components.linpol_p_list_start_indexes = k_linpol_p_list_start_indexes;
  components.linpol_p_list_comp_inds = k_linpol_p_list_comp_inds;

  components.n_linpol_angles = components.n_linpol_pol_frac + components.n_linpol_power + components.n_linpol_curve + components.n_linpol_p_list;
  components.intr_pol_angle = k_intr_pol_angle;
  components.rm_values = k_rms;
  components.linpol_angle_inds = k_linpol_angle_inds;
  components.do_QUV = 1;

  source_t *chunked_source = put_components_into_source(components, comptype,
                                                        n_powers, n_curves,
                                                        n_lists, num_coeffs);

  if (do_gpu == 1){

    do_gpu_calc_visi(n_powers, n_curves, n_lists,
          num_baselines, num_coeffs,
          num_freqs, num_visis, num_times,
          beamtype, comptype, &components,
          chunked_source,
          args_ft);

  }
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
void test_calc_visi_VarylmnVaryBeam(e_beamtype beamtype, e_component_type comptype,
                                    int do_gpu) {

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

  source_t *chunked_source = put_components_into_source(components, comptype,
                                                        n_powers, n_curves,
                                                        n_lists, num_coeffs);

  if (do_gpu == 1){

    do_gpu_calc_visi(n_powers, n_curves, n_lists,
          num_baselines, num_coeffs,
          num_freqs, num_visis, num_times,
          beamtype, comptype, &components,
          chunked_source,
          args_ft);

  }
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}




void test_calc_visi_VarylmnVaryPAMajMin(e_beamtype beamtype, e_component_type comptype,
                                        int do_gpu) {

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

  source_t *chunked_source = put_components_into_source(components, comptype,
                                                        n_powers, n_curves,
                                                        n_lists, num_coeffs);

  if (do_gpu == 1){

    do_gpu_calc_visi(n_powers, n_curves, n_lists,
          num_baselines, num_coeffs,
          num_freqs, num_visis, num_times,
          beamtype, comptype, &components,
          chunked_source,
          args_ft);

  }
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, comptype);

  free_args_for_testing( args_ft, components, comptype );
}

//This test varies the shapelet coeff params
void test_calc_visi_shape_VarylmnMultipleCoeff(int beamtype, int do_gpu) {

  int num_baselines = 10.0;
  int num_times = 5.0;
  int num_freqs = 3.0;

  int num_visis = num_baselines*num_times*num_freqs;

  int num_components = 25;
  int n_powers = num_components;
  int n_curves = 0;
  int n_lists = 0;
  int num_coeffs_per_component = 3;

  int num_coeffs = num_coeffs_per_component*num_components;

  //Container for many arrays to feed the GPU
  args_for_testing_t *args_ft = malloc(sizeof(args_for_testing_t));
  //Component information
  components_t components;
  //Allocate memory
  malloc_args_for_testing(args_ft, &components, num_baselines, num_times,
                          num_freqs, num_components, n_powers, n_curves, n_lists,
                          num_coeffs, SHAPELET);

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

    //Set major,minor to 3 arcmins
    components.pas[comp] = 0.0;
    components.majors[comp] = 3.0*(DD2R / 60.0);
    components.minors[comp] = 3.0*(DD2R / 60.0);

    components.power_comp_inds[comp] = comp;
  }

  //Make up some u,v,w values and scale by wavelength in correct order
  setup_uvw_and_freqs(args_ft, num_times, num_freqs, num_baselines);

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

  //Stick a number of coeffs in per component
  user_precision_t sign;
  count = 0;
  for (int comp = 0; comp < num_components; comp++) {
    for (int coeff = 0; coeff < num_coeffs_per_component; coeff++) {

      if (count % 2 == 0) {
        sign = 1.0;
      } else {
        sign = -1.0;
      }

      components.n1s[count] = count;
      components.n2s[count] = count + 1;
      components.shape_coeffs[count] = sign*1e-3*(coeff + 1);
      components.param_indexes[count] = comp;

      count ++;
    }
  }
  components.do_QUV = 0;

  e_component_type comptype = SHAPELET;

  source_t *chunked_source = put_components_into_source(components, comptype,
                                                        n_powers, n_curves,
                                                        n_lists, num_coeffs);

  if (do_gpu == 1){

    do_gpu_calc_visi(n_powers, n_curves, n_lists,
          num_baselines, num_coeffs,
          num_freqs, num_visis, num_times,
          beamtype, comptype, &components,
          chunked_source,
          args_ft);

  }
  //
  // //Check all results are within 0.1% of expected value
  // // double frac_tol = 1e-3;
  test_visi_outputs(num_visis, n_powers, n_curves, n_lists,
                    num_baselines, num_freqs, args_ft->extrap_freqs,
                    beamtype, args_ft, components, SHAPELET);

  free_args_for_testing( args_ft, components, SHAPELET );
}