sun_oe = {
    'N' : (20.552775873722602, 0.2136832880120676, -3.580978347278808, 174.97389301955806),
    'i' : (0.001101954782290859, 0.21369711394557486, -5.275796592097575, 0.00015489408592731231, 3.621325113811721e-07),
    'w' : (-20.630936449162093, 0.21368386463476155, -3.5858347681354386, 288.0401174521317),
    'a' : (0.0008991876561645148, 0.21276978132643024, 6.440250190274306, 1.000002237215457),
    'e' : (0.000624492722824836, 0.2299825556661901, 6.272769979941623, 0.01670036421395417),
    'M' : [ 9.85645883e-01, -6.84436431e+03]
    }

mercury_oe = {
    'N' : [2.5943137649076465e-33, -1.6458405468796637e-28, 4.1886116153189024e-24, -4.986022118098174e-20, 1.425664896833767e-16, 3.579902179195779e-12, -5.08222878798609e-08, 0.0003084130928174047, -0.9418197727279626, 1229.9985299396456],
    'i' : [-5.377291695654706e-36, 4.985297704367898e-31, -1.7171150986776112e-26, 2.7243783637920533e-22, -1.432580698949224e-18, -1.7317540599834052e-14, 3.4635949171992256e-10, -2.5036145236034457e-06, 0.00875602156565588, -5.3234684354500335],
    'w' : [7.266185018994313e-06, 29.130043647166488],
    'a' : (-7.302881780107657e-07, 0.015618668472226156, -0.16626523478656402, 0.3870983517945791),
    'e' : [-9.795222209711759e-79, 4.948344530013818e-74, -7.165646271347012e-70, -1.1487378401636948e-66, 6.912912897843176e-62, 4.067964777551289e-58, -5.014514038358126e-54, -8.163511441235564e-50, -8.754412696402857e-47, 8.607607039568674e-42, 8.339252306647572e-38, -3.0418308726218236e-34, -1.2193128129962863e-29, -5.338302719089814e-26, 1.1024320681178665e-21, 1.1731761292564323e-17, -9.765034347403027e-14, -1.4031494864886622e-09, 1.9017897264846455e-05, -0.08182346678242457, 125.15061900135669],
    'M' : [ 4.09233458e+00, -2.97113448e+04]
    }

"""
venus_oe = {
    'N' : ,
    'i' : ,
    'w' : ,
    'a' : ,
    'e' : ,
    'M' : 2
    }
mars_oe = {
    'N' : ,
    'i' : ,
    'w' : ,
    'a' : ,
    'e' : ,
    'M' : 2
    }


Fitting < N >
=============
mercury: polyfit(9) => 

venus: polyfit(9) => [8.883608712132133e-33, -5.6784062266386335e-28, 1.4587320685450877e-23, -1.759389614126735e-19, 5.272692894172828e-16, 1.2541053620125329e-11, -1.8200822892364454e-07, 0.0011218340122082697, -3.4777646043631467, 4507.8193294953935]

earth: : sin => A, w, p, c = (20.552775873722602, 0.2136832880120676, -3.580978347278808, 174.97389301955806)

mars: polyfit(9) => [-1.4975526200274625e-32, 9.320215280863113e-28, -2.3297249402101756e-23, 2.72213457640512e-19, -7.398179831242255e-16, -1.964590336281937e-11, 2.7406208795274616e-07, -0.001649340232889223, 5.011598248964769, -6221.330742340181]

jupiter: polyfit(9)  ZAIF => [-3.9126723651894454e-33, 2.545404735881869e-28, -6.662511597738446e-24, 8.218648240866026e-20, -2.6341532185757007e-16, -5.800086268786816e-12, 8.680503269243621e-08, -0.000545212498827264, 1.717713797579897, -2120.748293944222]

saturn: polyfit(9) + sin (raje be tarkibeshoon fekr shavad)
felan polyfit(9) => [2.3773816078678274e-32, -1.5117086306563767e-27, 3.8613024695068167e-23, -4.623258134068388e-19, 1.3520049165925946e-15, 3.307962858502426e-11, -4.7465986506670604e-07, 0.0029027291020599826, -8.930828598738573, 11404.707504343633]

uranus: polyfit(9)  ZAIF => [7.899427519273183e-34, -6.817088781653662e-29, 2.278277330067026e-24, -3.5832671617281494e-20, 1.9043231302967451e-16, 2.255924803400395e-12, -4.6024922011574896e-08, 0.0003378036094103372, -1.1995235218981393, 1788.3740274570002]

neptun: polyfit(9) + sin (raje be tarkibeshoon fekr shavad)
felan polyfit(9) => [-2.2082474073603896e-32, 1.3640095288667195e-27, -3.3808512658213116e-23, 3.9090801393282146e-19, -1.0259972318391186e-15, -2.8317482020125772e-11, 3.8999199075437446e-07, -0.0023302443253615064, 7.045351475431107, -8656.072890706693]


Fitting < i >
=============
mercury: polyfit(9) => 

venus: polyfit(9) => [1.4064221204108598e-34, -8.521753182907468e-30, 2.0605705954325883e-25, -2.2962854986367836e-21, 5.114135925769306e-18, 1.6976569276153283e-13, -2.191256592521741e-09, 1.2437865744731987e-05, -0.03561566031280168, 45.19853435473836]

earth : sin_trend => A, w, p, c, m = (0.001101954782290859, 0.21369711394557486, -5.275796592097575, 0.00015489408592731231, 3.621325113811721e-07)
A * np.sin(w*x + p) + c + m*x

mars: polyfit(15) => [1.171082541673678e-56, -6.338745899831259e-52, 1.150153523361735e-47, -4.0254348399981885e-44, -9.503105111836479e-40, 3.800172652156616e-36, 1.0362454927085267e-31, -1.4355677831526594e-28, -1.145970368538569e-23, -8.777558325428042e-21, 1.2506985066914e-15, 7.986872775884659e-14, -1.4189843705798802e-07, 0.001207227680474841, -4.228619762057916, 5645.854121145258]

jupiter: polyfit(9)  ZAIF => [-5.50978269373824e-35, 3.4289524099446574e-30, -8.532565895792388e-26, 9.851714039633446e-22, -2.510067820400572e-18, -7.186139291827545e-14, 9.711153163557227e-10, -5.683263750540564e-06, 0.01672833308414061, -18.86016082715479]

saturn: polyfit(9) ZAIF => [7.139547770569372e-34, -4.482524945297787e-29, 1.126793710186633e-24, -1.3189365249962781e-20, 3.544949649985849e-17, 9.552257603628643e-13, -1.3199208659456419e-08, 7.856736276550627e-05, -0.23543876261153462, 292.13896950148387]
uranus: polyfit(9) ZAIF => [2.7135196267729225e-34, -1.6926268122936698e-29, 4.226201057505193e-25, -4.906008879808542e-21, 1.2801638381137727e-17, 3.566311014037476e-13, -4.870049167944941e-09, 2.8753274749679968e-05, -0.08550801953627063, 105.18270771617652]
neptun: polyfit(9) ZAIF => [-4.666138640598417e-34, 2.897349913668551e-29, -7.231445544755431e-25, 8.449031646523512e-21, -2.3137647768477404e-17, -6.083585360498975e-13, 8.536098347815676e-09, -5.172067308383399e-05, 0.15856660274207288, -198.89793203544457]


Fitting < w >
=============
mercury: polyfit(1) => 
venus: polyfit(20) => [3.885724519806261e-73, -1.6805644141761093e-68, 1.929439129881774e-64, 7.812718352683747e-61, -1.5921442042476676e-56, -1.4821551680656297e-52, 6.610922236196345e-49, 2.0023949729735412e-44, 8.917924029651088e-41, -1.3918385119037133e-36, -2.099452830930911e-32, -1.7767295899649145e-29, 2.147006898584575e-24, 1.657208323216468e-20, -1.3988074666227706e-16, -2.418704223238284e-12, 1.0269141962488598e-08, 0.00026047171826611663, -2.7316389611311918, 10215.801090948815, -13901043.735505505]
mars: polyfit(20) => [2.1500017719911625e-73, -9.679463108192685e-69, 1.1902054301100992e-64, 3.896062287935365e-61, -1.0353500574437038e-56, -8.39755261335942e-53, 5.422334379176946e-49, 1.2752020643953513e-44, 4.247340066257932e-41, -1.0386134594843835e-36, -1.326903081098734e-32, 7.970929399412635e-30, 1.5448437549615103e-24, 9.854915301263418e-21, -1.1612138735341645e-16, -1.6381951714190799e-12, 9.339325245242289e-09, 0.00018333771964727157, -2.133901868843252, 8496.051724752697, -12205886.626958573]
jupiter:  polyfit(9) ZAIF => [1.1269343458884313e-31, -7.291541333318378e-27, 1.9040372569979172e-22, -2.352530042469386e-18, 7.672553471365106e-15, 1.651918807396791e-10, -2.5046175108791134e-06, 0.01592301325955861, -50.87891359973493, 67158.25516429092]
saturn: polyfit(9) => [4.055222987049113e-31, -2.5695230677633402e-26, 6.516906182277206e-22, -7.698561172392914e-18, 2.118123253394868e-14, 5.567664460323397e-10, -7.743198365595809e-06, 0.046134694561964215, -137.9087744644863, 169024.27988371646]
uranus: polyfit(9) => [-1.7434576965114308e-31, 8.474645945365895e-27, -1.4376471167490101e-22, 6.452042696099745e-19, 8.179268688250957e-15, -7.995715964526355e-11, -3.814997484286599e-07, 0.008226741096739536, -41.04459615492501, 70334.01373715977]
neptun: polyfit(9) => [1.8012364389512415e-30, -1.1190744959252465e-25, 2.798465116404526e-21, -3.2817781564957556e-17, 9.123787577033405e-14, 2.3583206957456476e-09, -3.330820912722121e-05, 0.20271219517866695, -624.0508766329033, 793427.4885507486]


Fitting < a >
=============
Sin is the best fit for <a> of ALL planets. (A, w, p, c)

mercury = 
venus = (-4.398726985167627e-06, 0.021535750705983677, 7.38795406021702, 0.7233297518634283)
mars = (6.99973029538498e-05, 0.015286987960212791, 1.1327289006004275, 1.5236769001329489)
jupiter = (0.0007842230849913618, 0.001335336582944745, 11.201102489430493, 5.203064043289066)
saturn = (0.03270543872547703, 0.0009886647253105937, 6.415190020182706, 9.549749229808656)
uranus = (0.06841859575729173, 0.0014300561641803187, 1.2665566674296087, 19.234972726608998)
neptune = (0.14897430159014052, 0.001394513189679602, 2.6802584975116233, 30.15126475563609)

sin_func = lambda t: A * np.sin(w*t + p) + c

Fitting < e >
=============
mercury : polyfit(20) => 
venus : polyfit(20) => [-1.0381777086897013e-76, 4.655154633928052e-72, -5.689854175005972e-68, -1.8959279377799044e-64, 4.9285125355206404e-60, 4.0370676752167875e-56, -2.5458926933730127e-52, -6.076432248109977e-48, -2.064776467667157e-44, 4.905821239612119e-40, 6.32215778587972e-36, -3.328436567169255e-33, -7.314701325241423e-28, -4.705581172143132e-24, 5.474258549178694e-20, 7.777640567573351e-16, -4.398865291026659e-12, -8.695678283345249e-08, 0.0010099689841294244, -4.01899261968326, 5774.214226892281]
mars : polyfit(20) => [-1.3623100489685777e-76, 5.72277706676604e-72, -6.2238227176251686e-68, -2.9201654825163845e-64, 4.902832678884259e-60, 5.097945284843599e-56, -1.553298324318741e-52, -6.27873222653364e-48, -3.404826583011182e-44, 3.717711290256498e-40, 6.621429636478276e-36, 1.3657554030478848e-32, -5.983579650005481e-28, -5.47995537432202e-24, 3.2615330431158544e-20, 7.16836889334241e-16, -2.0678425662730782e-12, -7.438769342568587e-08, 0.0006947387223147139, -2.392561773988341, 3001.329122260517]
jupiter: polyfit(5) ZAIF => [1.7676070455155284e-20, -8.098581027462522e-16, 1.4745778474496384e-11, -1.3330352822954084e-07, 0.0005978825782425681, -1.0150303791616306]
saturn : polyfit(5) ZAIF => [1.0368080003684683e-20, -2.452298391091555e-16, 1.3543453418128338e-13, 3.836301675648586e-08, -0.0003469947023894754, 0.9685756181402778]
uranus : polyfit(5) ZAIF => [1.4253164639927094e-19, -6.47519304948677e-15, 1.1629406374485125e-10, -1.0310553067570207e-06, 0.0045086344678756315, -7.728070475847356]
neptune : polyfit(5) ZAIF => [-1.3638940844223904e-19, 6.3211333605446095e-15, -1.156524081388617e-10, 1.0425637568273344e-06, -0.004624577518673766, 8.078252052987942]

Fitting < M >
=============
mercury: polyfit(1) => [ 4.09233458e+00, -2.97113448e+04]
venus: polyfit(1)   => [ 1.60212135e+00, -1.14719551e+04]
earth: polyfit(1)   => [ 9.85645883e-01, -6.84436431e+03]
mars: polyfit(1)    => [ 5.24018817e-01, -3.58138004e+03]
jupiter: polyfit(1) => [ 8.30064336e-02, -3.39025711e+02]
saturn: polyfit(3)  => [ 8.25642457e-10, -2.33880342e-05,  2.50880835e-01, -7.04167266e+02]
uranus: polyfit(5)  => [-1.37442630e-16,  6.87312937e-12, -1.35804129e-07,  1.32377262e-03, -6.34818463e+00,  1.21847232e+04]
neptun: polyfit(5)  => [-9.79520106e-18,  1.02238328e-12, -2.43680314e-08,  2.30922525e-04, -9.12446413e-01,  1.46442642e+03]
"""