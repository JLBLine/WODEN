#include "woden_precision_defs.h"

int num_powers = 6;
int num_curves = 6;
int num_lists = 6;
int num_extrap_freqs = 25;
user_precision_t ref_stokesI[] = {5.0179772520326402,21.7233756369023290,82.2103892979358761,13.2978162103577713,99.1283179729780812,54.4223450765814860};
user_precision_t ref_stokesQ[] = {98.0736625153917601,10.0424902486137260,54.8358073648647490,86.4288064069258581,20.6079017356986149,15.9193471780000468};
user_precision_t ref_stokesU[] = {65.5048270666696766,85.2912544824885543,74.4989584335259849,99.6087707648812142,77.0206604297566741,42.2871399696728432};
user_precision_t ref_stokesV[] = {90.2256259763843360,92.5075374231781353,28.0053834234872134,7.1744788905378849,98.9061604540046346,27.5073258526496645};
double extrap_freqs[] = {5.0000000000000000e+07,6.0416666666666664e+07,7.0833333333333328e+07,8.1250000000000000e+07,9.1666666666666657e+07,1.0208333333333333e+08,1.1250000000000000e+08,1.2291666666666666e+08,1.3333333333333333e+08,1.4375000000000000e+08,1.5416666666666666e+08,1.6458333333333331e+08,1.7500000000000000e+08,1.8541666666666666e+08,1.9583333333333331e+08,2.0625000000000000e+08,2.1666666666666666e+08,2.2708333333333331e+08,2.3750000000000000e+08,2.4791666666666666e+08,2.5833333333333331e+08,2.6875000000000000e+08,2.7916666666666663e+08,2.8958333333333331e+08,3.0000000000000000e+08};
double ref_freqs[] = {5.0000000000000000e+07,1.0000000000000000e+08,1.5000000000000000e+08,2.0000000000000000e+08,2.5000000000000000e+08,3.0000000000000000e+08};
user_precision_t ref_power_SIs[] = {-0.1022213328764028,-0.5603205193197343,-1.0306542637199814,-0.0882522570164057,0.4444977423552279,-1.1954750626111341};
user_precision_t ref_curve_SIs[] = {36.5944497812515337,55.3879959282168173,49.5959293153752014,29.5432628638802619,46.8648888586917778,11.8518790482611927};
user_precision_t ref_qs[] = {-0.9770475382131205,-1.4936985186850849,-1.2990771564856129,-0.7946912886162745,-1.2313969572622310,-0.3122362036366431};
int num_list_values[] = {11,8,5,23,21,5};
int list_start_indexes[] = {0,11,19,24,47,68};
double list_freqs[] = {5.0000000000000000e+07,7.5000000000000000e+07,1.0000000000000000e+08,1.2500000000000000e+08,1.5000000000000000e+08,1.7500000000000000e+08,2.0000000000000000e+08,2.2500000000000000e+08,2.5000000000000000e+08,2.7500000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,8.5714285714285716e+07,1.2142857142857143e+08,1.5714285714285713e+08,1.9285714285714287e+08,2.2857142857142860e+08,2.6428571428571430e+08,3.0000000000000000e+08,5.0000000000000000e+07,1.1250000000000000e+08,1.7500000000000000e+08,2.3750000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.1363636363636360e+07,7.2727272727272719e+07,8.4090909090909094e+07,9.5454545454545453e+07,1.0681818181818181e+08,1.1818181818181819e+08,1.2954545454545455e+08,1.4090909090909091e+08,1.5227272727272725e+08,1.6363636363636363e+08,1.7500000000000000e+08,1.8636363636363637e+08,1.9772727272727272e+08,2.0909090909090909e+08,2.2045454545454544e+08,2.3181818181818181e+08,2.4318181818181819e+08,2.5454545454545453e+08,2.6590909090909091e+08,2.7727272727272725e+08,2.8863636363636363e+08,3.0000000000000000e+08,5.0000000000000000e+07,6.2500000000000000e+07,7.5000000000000000e+07,8.7500000000000000e+07,1.0000000000000000e+08,1.1250000000000000e+08,1.2500000000000000e+08,1.3750000000000000e+08,1.5000000000000000e+08,1.6250000000000000e+08,1.7500000000000000e+08,1.8750000000000000e+08,2.0000000000000000e+08,2.1250000000000000e+08,2.2500000000000000e+08,2.3750000000000000e+08,2.5000000000000000e+08,2.6250000000000000e+08,2.7500000000000000e+08,2.8750000000000000e+08,3.0000000000000000e+08,5.0000000000000000e+07,1.1250000000000000e+08,1.7500000000000000e+08,2.3750000000000000e+08,3.0000000000000000e+08};
user_precision_t list_stokesI[] = {-1.8549529265775067,-3.9146507921992937,-4.4764012402425921,9.3600416205224981,-1.7012471614906346,8.4011848681654904,9.9889314112354484,9.6281164394528460,-3.0903934832508413,3.8111338207483278,4.8932999885235304,-2.7171329103975213,-2.1617886475103081,5.0442735387167694,0.9043630209790345,-4.8348123500841025,6.3374367439529244,5.6070859466615861,2.9300980981860132,0.9952063527081751,4.9260406557371628,-0.5549443035474502,3.9539627020484556,-0.7990320616730493,-1.1697145708651813,3.7618578783085397,-0.6741937500734183,5.4019243809650490,2.2414784733953983,7.7902457478980480,-4.1064011347182223,4.9584376992743682,4.4201512691895495,2.2220323980681069,-0.6871756542480822,-3.7079952325402039,9.5194285396986569,2.0008394293323679,3.1285956297268491,-3.8736429424824914,0.8960397079698943,5.9848375542456846,-0.9946295352614190,6.1078168121869272,4.0995527266753431,7.4016132113030757,1.3793442387145207,1.7841063184638255,2.3409575119150876,4.8744045869740589,0.6013909817416438,-1.2465348253029118,-3.8179482435649059,5.3782910502322032,4.0234111947026925,8.8756965539709061,1.8199376997102981,1.2292479805234038,4.9973755619985596,9.0035199629633684,-2.5926042665326330,3.1976165374205578,1.5867750719898419,9.3901217271314454,1.0730176075327540,-2.2782419258915332,3.9956497181961392,-2.7116615481619881,8.3499890336995097,-3.7242962710774070,0.3034882539154999,7.0188131130816043,7.3976234467145101};
user_precision_t list_stokesQ[] = {-4.7822840232832018,8.2730486450783918,5.8318297651195685,3.9821773296995584,6.4790620236773115,0.7476316870510864,9.1542454768245296,0.7109167209660914,-2.1489579976331248,0.1044310874668790,3.4409280779554052,9.8102944696061165,9.0682389874886571,1.3051221908380013,8.9607197005714685,3.8184097595737025,9.6623467468046371,-2.6682126452396036,-1.5949924172479690,0.1573659386105799,8.3885180301977087,6.3694911646350523,-1.6442316774411481,3.6917087124979240,-0.9096062692032483,6.7686466488446335,-4.7418116592198247,3.4500833052428366,7.3341746769672600,9.8141079761126306,8.2100629413165187,6.6545682687751793,3.8221297608906681,4.6376122920987406,4.2006024058575022,2.8778548763404173,3.3873045281773457,4.8161042147507516,4.4740880561866838,5.6155331372919601,9.0201654289760445,-1.9557063396026795,0.0321925945428907,7.1381754762706269,1.5321063594148683,3.3629492766171438,-3.6975707426406625,2.6980703501222711,1.3630337613237051,7.0685553318323020,6.3490902219013901,9.7032785476199059,2.9577996090871546,-2.9083764116494937,7.2253454061679001,3.6642499594833371,2.6376866299704806,6.7320070447456821,2.6858728994599055,-3.9799605775728164,6.3413493552586218,-4.6402654954349867,5.8888999849235351,-3.1599996605677187,2.7484910190718645,8.7880652454800448,-1.0526651965630824,-2.5561343655354407,5.2451597942800099,-3.3282607560366739,6.7017125146034449,-2.3719992743459626,1.6586968738764210};
user_precision_t list_stokesU[] = {7.3715348094993303,1.2531343585464425,2.7443892443449824,1.0418164215999290,-3.7517940491737347,4.7334105053683597,1.6249819384983635,-4.0729062569967009,6.1527408889986965,4.6725828764202770,4.1699096643498716,1.4846594152147139,-0.0078429497509900,0.4380556245755880,-4.4288811567337731,4.2595102371836795,-1.0346368710698268,7.3886233639465129,4.0254529824575300,5.6459197838626558,-3.5283763841330544,6.1480845955430432,-2.0492625727127893,5.9848571401735473,3.4317047218058825,0.6783810734764169,5.5166160188445357,5.0044258706748508,5.9149365242108782,4.9019493958668416,6.1953980341138752,0.9407315314685958,-3.5247985613521955,8.6518600934964525,0.1043837350628598,-4.7291291029742775,-1.2078224127345285,0.7888286783103586,4.0880880003405391,-4.2079960207057931,1.4887579313724721,7.2163489252292798,5.4235380328748768,3.2038679114385662,0.5908040427944954,-0.0038235938018696,-1.3630376074508392,-4.9065530614043587,7.9561689442210088,3.9477729831424195,2.6085734689282045,8.8028481158576195,-0.6958646651127456,-2.1343705685341643,4.0437369547653024,3.9188808880173553,5.4401708809639171,-4.6904453709564358,-1.1398405564878491,-1.5738394855723130,9.0661341220046943,1.0355495203078391,-0.6765430890357642,-0.5495843129188920,4.6386723348110159,5.0518865051215620,5.8608830689077944,-2.8886421360891630,2.0796794188225087,-2.4847576473111448,-0.0182437068873229,0.9798951947009282,8.6026311776012356};
user_precision_t list_stokesV[] = {9.6779525654360263,9.7964497618843893,8.2316635714398334,1.2801834921363566,3.3152249271884635,9.6468971605099085,1.1765725749508080,-4.5763595842097535,9.0457219563896718,3.2474146476836943,-1.4715592497009227,8.1020996091048652,9.9331393752485653,8.4808485110023089,-0.4830788558716224,0.1727638573408798,8.6274780381179870,1.1835607370636598,4.3043503660149725,2.4563766924755654,7.2535314087063814,-3.5117001706130839,9.6568999938781310,1.3775281729512461,-3.3137602114892100,3.6735724993874950,7.4327920297000443,-3.4184830852876851,9.7671840072578000,0.1086048002491786,0.7207332485471616,-2.3322165802336698,-2.9765776025954809,-3.4878741654605880,6.7778887457302393,-0.6358341514877255,3.9915450747627599,-4.3549102137209363,7.1427194907506415,9.9900271970821617,4.0905813734066854,5.3525209822488016,1.5027135009276007,4.2213345270935143,-0.8422549951053862,4.4702598185889855,2.9834180682538110,8.1234707878957479,7.6735803157371549,6.0268828298551753,-3.7919198601949375,9.0606973584596862,4.4072514545053991,-1.6768281264776470,3.3322677827514422,0.0414361397751097,-0.1675207443727960,9.2832711159128269,1.1930281576324582,4.0563077212476202,1.8544866954748844,-1.7214641787842395,2.1269776132351490,3.8822961894981685,9.5280300892136722,4.5861624275783015,3.4825902757741183,-4.3889978670939875,5.1997771241435444,8.0428437289893360,7.7125223386116097,-0.2286726247015745,6.6738877719948970};