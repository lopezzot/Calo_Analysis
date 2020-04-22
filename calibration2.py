
def calibscin(Barrel_VectorSignals):
	"""function to apply calibrations at all the towers (right and left is the same)
	for scin signals"""
	#scincalibrations = [415.46911408615017, 415.49032413530045, 414.2257295618793, 413.52149204021396, 411.5992702498399, 410.64900009885315, 410.05084897192495, 410.2933477391751, 411.3409180068102, 412.31065478145797, 410.66077419019877, 411.16318351136084, 411.44777420544165, 411.9325417788443, 411.37997153809016, 411.388755835253, 411.8810225437739, 411.58323678675737, 411.40360405115445, 411.4363536068476, 411.1434892175137, 411.0311309144043, 411.45278182281913, 411.295529214398, 410.66722584949844, 411.95872061830596, 411.5134821865462, 411.8234427877747, 411.0133769672026, 410.49154689961557, 411.8200149562403, 411.4650824673985, 410.2795836068218, 410.99784486503626, 411.3719137475675, 411.1496647172716, 411.10361391486657, 411.01421880697376, 410.597406149097, 410.8673343110665, 410.4613021375798, 410.9360845492963, 410.9528787632441, 410.8574880829973, 411.03315409936096, 410.9379153502092, 411.01901911573697, 410.80356948239154, 411.0444263559407, 410.84313941730574, 410.94072813786545, 410.204939014593, 410.9738410909085, 410.3391884191118, 410.6721258934903, 410.51534130232534, 409.62060911658847, 408.77108688926836, 409.99777980214225, 409.5992379372193, 409.83732253845403, 410.210413105926, 409.1679944501218, 409.7662734551919, 408.11740676906703, 408.3924355876955, 408.03209385093214, 407.6136434771539, 407.49508269173134, 407.10864535201216, 406.49601594991765, 407.97586047362637, 407.6622587467224, 407.70767503648653, 408.04075414921584]*36	
	#scincorrections = [411.41483745459146, 411.52181636082577, 410.25098543903624, 409.76481640097643, 407.8648162057763, 407.0433073672767, 406.4280147748465, 406.6202788727341, 407.6662499733772, 408.6662149874733, 407.04700222023894, 407.5469707080265, 407.8482602981088, 408.34396984724685, 407.8356135111425, 407.9006104165864, 408.4005126486694, 408.08624713560107, 408.0027574234192, 408.03512549563607, 407.7646802314756, 407.66902911324036, 408.13219631572156, 408.0856386283693, 407.38021535979243, 408.76602535322274, 408.3354792012611, 408.7493505933448, 407.97452518970755, 407.4711792144778, 408.78371452969293, 408.59189054866147, 407.3678692910586, 408.1980036194848, 408.6716381355258, 408.4772494588298, 408.48623697063414, 408.4283387274703, 408.0399995983036, 408.459433126158, 408.0009678326388, 408.4220840256861, 408.401124384891, 408.22498572575756, 408.2946388759221, 408.16241045483764, 408.2162208231712, 407.86038369360386, 408.10072483636185, 407.8308977000335, 407.8739032108822, 407.05402921668883, 407.81948812803313, 407.0379965972104, 407.2807418021031, 407.2472061094793, 406.2461313455679, 405.3115718891361, 406.46924320089, 406.04385491758575, 406.23962673351735, 406.485581819831, 405.37184790415176, 405.9200835792759, 404.174107208703, 404.3861928079304, 403.94869689782337, 403.3611910496827, 403.06973154521165, 402.5554299506423, 401.6688569907378, 402.8571098242905, 402.03193211770736, 401.6356578529606, 401.29776188137134]*36
	#containment = [0.9922257350540894, 0.9924029080295286, 0.9925314933889849, 0.9923710947914168, 0.9925005601234895, 0.9924571682663081, 0.9924736876499678, 0.9923107627997272, 0.9923717693883665, 0.992462552656521, 0.9924100067589308, 0.9924641871673259, 0.9923359607471853, 0.9924773952914967, 0.9924178141593918, 0.9924772080326696, 0.9922415412669103, 0.992379952512389, 0.9924161038441175, 0.9923936205889852, 0.9924236863162547, 0.9924316563576252, 0.9923568814774528, 0.9924405399594359, 0.9924611132742223, 0.9924930212297, 0.9923853174427201, 0.9923924066611056, 0.9924806854716944, 0.992390053683815, 0.9923893436872515, 0.9924596118186407, 0.9924440389467973, 0.9923388445611183, 0.9924589198963458, 0.9924056855957695, 0.9923602278246142, 0.9924218711993472, 0.9922791098871626, 0.9924370426348486, 0.9923948825562914, 0.9924726323036805, 0.9924027753718472, 0.9924520521026018, 0.9924006741483085, 0.9923441256883894, 0.9925201121485623, 0.9924546188776329, 0.9924782088159386, 0.9923947843091397, 0.9924698150227844, 0.9922952793358721, 0.9924636652757733, 0.9924582419298134, 0.9924127350793948, 0.9924450958333008, 0.9924155137751113, 0.992411906630523, 0.9924154871234933, 0.992445097600654, 0.9924117857125934, 0.9924140566552296, 0.9924301803483659, 0.9924252926230659, 0.992427658103939, 0.9924930313657724, 0.9925539612633278, 0.9924001191478654, 0.9924655317037001, 0.9924423192903769, 0.9925128053546892, 0.9925077413435955, 0.9923661833622193, 0.9920001139857957, 0.9762479738306956]
	s_cont = [408.21638950554075, 408.3954472740771, 407.1870232421094, 406.63875945884087, 404.8060585388971, 403.97304819147996, 403.3691105878475, 403.49367909804056, 404.55647780600043, 405.58591491094637, 403.9575182245898, 404.4757730162475, 404.72249522199195, 405.272159576985, 404.74332809708255, 404.83205898107536, 405.23195412471205, 404.9766105533868, 404.9085068798063, 404.9314555180952, 404.67532710488985, 404.58364980855805, 405.012793566413, 405.0007315500301, 404.30902206187204, 405.6974274788762, 405.2261341502687, 405.63975175649347, 404.90683641527, 404.37034541526305, 405.67260217215875, 405.5109490861691, 404.2898135363692, 405.07073526391474, 405.58981257625425, 405.3751447994642, 405.36549518339785, 405.3332161707569, 404.88956759976287, 405.37027184803094, 404.8980725551248, 405.34774082392767, 405.2984093045488, 405.14372480308344, 405.19187487160525, 405.03757034167137, 405.16280927227615, 404.7829216539207, 405.03107640207867, 404.7292557576276, 404.8025372723253, 403.9177916263665, 404.7460239584375, 403.96821450150077, 404.1905949169899, 404.1704924951662, 403.16496315846314, 402.2360298379118, 403.3863719919289, 402.9762332238292, 403.15699339382735, 403.4020052256797, 402.3032561236677, 402.8453577277423, 401.11356268338346, 401.3504783424065, 400.94087925309395, 400.29569405733, 400.0328154316862, 399.5130445431503, 398.66148407548866, 399.83880015591535, 398.96289406538807, 398.42261837089694, 391.76612693948175]*36

	#sequence is: scincalibrations are calibrations estimated at first by firing electrons 
	#containment is containment tower by tower
	#scincorrections are the corrections to apply to get correct calibration constant to reconstruct energy deposited (s vector)
	#s_cont is s vector corrected by containment to reconstrct primary electron energy
	
	scincalibrations = [0.1]+s_cont
	if len(Barrel_VectorSignals) != len(scincalibrations):
		print "wrong calibration length s"+str(len(Barrel_VectorSignals))+" "+str(len(scincalibrations))
		quit()
	
	calib_Barrel_VectorSignals = [Barrel_VectorSignals[counter]*(1/(entry)) for counter, entry in enumerate(scincalibrations)]
	return calib_Barrel_VectorSignals

def calibcher(Barrel_VectorSignalsCher):
	"""function to apply calibrations at all the towers (right and left is the same)
	for cher signals"""
	#chercalibrations = [106.37048257816659, 106.13146567275112, 105.86590903712592, 105.81634772611199, 105.70865247970764, 105.64190045805486, 105.6639555914862, 105.5283776283045, 105.65723353008754, 105.7715133071762, 105.5658607357618, 105.61417315615921, 105.70536383822864, 105.77758372693368, 105.77043752701637, 105.5506909020155, 105.70623524436972, 105.65209417795693, 105.73892093109146, 105.61237691382905, 105.61860853169827, 105.63069081304641, 105.63120977003204, 105.72883741949272, 105.64351076067801, 105.51191943482567, 105.75579206781183, 105.6255099492349, 105.52676483799216, 105.50729849902562, 105.49691196235294, 105.51232976499537, 105.49294739086051, 105.46743366546217, 105.39944256662419, 105.43093327017503, 105.42218126102728, 105.33353578857226, 105.26156500953066, 105.19128846086149, 105.31605046745744, 105.35904116180494, 105.31714780361818, 105.34284763793717, 105.31700630345172, 105.32123783068911, 105.3888770522753, 105.37422214386339, 105.52604928728498, 105.46278011845477, 105.28372058763158, 105.4206155499903, 105.44174663935398, 105.56551159861326, 105.43779145304754, 105.54146571740893, 105.55442493693165, 105.67687272835073, 105.64965481223994, 105.65043713995333, 105.39546198539618, 105.5315247334327, 105.46869729570527, 105.45532655577286, 105.45938285691861, 105.49087149606235, 105.44596875109438, 105.52868954952874, 105.53867506889421, 105.69754494826158, 105.86874369528533, 105.70806060408836, 105.80369193129457, 106.1513324795608, 106.22235486676682]*36
	#chercorrections = [103.89550278429044, 103.70085240913916, 103.47143663620487, 103.40757852713733, 103.32416049481843, 103.2595513611303, 103.2771449523405, 103.14870022496955, 103.26744048386126, 103.40758019116339, 103.21664103253485, 103.25622797189547, 103.34575627566997, 103.44938982761327, 103.45740559874831, 103.25964142067909, 103.36696644454749, 103.34908426908248, 103.46274313918914, 103.34676995410369, 103.38553010427495, 103.41998384923339, 103.4223873054996, 103.54772829312716, 103.47303521689206, 103.38099295102664, 103.64838254840568, 103.52882284094234, 103.47351928457604, 103.45766349134928, 103.4805975522724, 103.53959817237843, 103.54212246362432, 103.52461109374697, 103.49938164982743, 103.56077151211231, 103.57780504091401, 103.53601233216357, 103.49962719536036, 103.47251867627905, 103.57132307678117, 103.59517284094395, 103.49923289942933, 103.51277318075437, 103.43104568143158, 103.40570108524098, 103.44047928171713, 103.39003361389342, 103.5289699793597, 103.41222913634279, 103.2043558306036, 103.30902132105307, 103.29948932213777, 103.37716306057568, 103.24111809153652, 103.31134211867192, 103.25381251814082, 103.35638287866391, 103.34984997105141, 103.24469534151174, 102.97376718686942, 103.06332233060792, 102.9776969458658, 102.9126043760532, 102.87127991738039, 102.82959630902663, 102.72936950682644, 102.74096868887399, 102.625968873716, 102.70935760664186, 102.73558010993771, 102.45477830507473, 102.33168129960764, 102.49598661880496, 102.19939271116233]*36
	c_cont = [103.08779161895677, 102.91302749597065, 102.69865952763615, 102.61869191270468, 102.54928716539662, 102.48068194031679, 102.49984890080964, 102.35556540203991, 102.47969263317724, 102.6281510005559, 102.43322742473204, 102.47810836409134, 102.55371034296142, 102.67118096060427, 102.67297232291142, 102.48284061965019, 102.5649981010228, 102.56155933915096, 102.67809243921879, 102.56067521092992, 102.60224889784466, 102.63726587197354, 102.63191774143888, 102.76496337880408, 102.6929637252195, 102.60491403169074, 102.85913301772406, 102.741217657914, 102.69546934772463, 102.67035622618218, 102.69304228926421, 102.75886941001674, 102.75976221892324, 102.731492956408, 102.7188845221274, 102.77429845330465, 102.78649420797491, 102.75140309520445, 102.70051794706535, 102.68996042906552, 102.78365100098196, 102.8153738834064, 102.71292597825087, 102.73146416207084, 102.6450394621172, 102.61404003462839, 102.66675609739092, 102.60991640602225, 102.750246685674, 102.62575682868824, 102.42720794074478, 102.51305416968992, 102.52098979376447, 102.59751750679058, 102.45780037787654, 102.53083482963227, 102.47068539942974, 102.5721049950492, 102.56599170316093, 102.46469174495641, 102.19238017547394, 102.28148980648412, 102.19817435184497, 102.1330715125064, 102.09230341456059, 102.05765775486448, 101.9644426420847, 101.96014956820567, 101.85273676485993, 101.93311307596035, 101.96637882465569, 101.68716060542853, 101.55050000833062, 101.67603040894112, 99.77195006099979]*36

	chercalibrations = [0.1]+c_cont
	if len(Barrel_VectorSignalsCher) != len(chercalibrations):
		print "wrong calibration length \n"
		quit()

	calib_Barrel_VectorSignalsCher = [Barrel_VectorSignalsCher[counter]*(1/(entry)) for counter, entry in enumerate(chercalibrations)]
	return calib_Barrel_VectorSignalsCher