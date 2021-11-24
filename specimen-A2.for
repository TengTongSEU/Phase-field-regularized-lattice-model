*     *************************************************************************************************
*
*     Common variable module, which defines some common variables 
*     and allocates memory
*
*     *************************************************************************************************
*      
      module vars_module
          parameter (NumEle= 10000)
          real*8,save :: allD(NumEle),allH(NumEle)
          real*8,save :: allArea(NumEle),allFt(NumEle)
          real*8,save :: allsigLocX(NumEle),allsigLocY(NumEle)
          real*8,save :: allstressLocX(NumEle),allstressLocY(NumEle)
          real*8,save :: allLE11(NumEle),allLE22(NumEle)
      end module

      
*     *************************************************************************************************
*
*     VUMAT for Visulization 
*
*     *************************************************************************************************     
*      
      subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      use vars_module
      include 'vaba_param.inc'
C
C All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(*), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),
     4  stretchOld(*), defgradOld(*),
     5  fieldOld(*), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(*),
     8  stretchNew(*), defgradNew(*), fieldNew(*),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
      parameter( zero = 0., one = 1., two = 2., three = 3.,
     1  third = one/three, half = .5, twoThirds = two/three,
     2  threeHalfs = 1.5 )
      
      integer jElemUid
*     ================================================================================================     *
*     Statev(1)    allsigLocX
*     Statev(2)    allsigLocY
*     Statev(3)    allstressLocX
*     Statev(4)    allstressLocY
*     Statev(8)    element identification number                                                           *
*     ================================================================================================     *
*     Parameter initilization
*
      e = 1.0d0
*
      do 100 i = 1,nblock
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         stress calculation
          sig1   = stressOld(i,1) + e*strainInc(i,1)
          stressNew(i,1) = sig1 
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         element number identification
          jElemUid      =  stateOld(i,8)
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         area input

          if (jElemUid .eq.	1	) area=	0.796592985
          if (jElemUid .eq.	2	) area=	1.072580516
          if (jElemUid .eq.	3	) area=	0.353314948
          if (jElemUid .eq.	4	) area=	2.347072153
          if (jElemUid .eq.	5	) area=	2.147966543
          if (jElemUid .eq.	6	) area=	0.451549783
          if (jElemUid .eq.	7	) area=	2.213739955
          if (jElemUid .eq.	8	) area=	2.308975817
          if (jElemUid .eq.	9	) area=	0.724919
          if (jElemUid .eq.	10	) area=	1.503885273
          if (jElemUid .eq.	11	) area=	1.56876778
          if (jElemUid .eq.	12	) area=	1.344942973
          if (jElemUid .eq.	13	) area=	1.641663527
          if (jElemUid .eq.	14	) area=	1.81829696
          if (jElemUid .eq.	15	) area=	1.910774489
          if (jElemUid .eq.	16	) area=	2.529642502
          if (jElemUid .eq.	17	) area=	2.726847367
          if (jElemUid .eq.	18	) area=	0.668224552
          if (jElemUid .eq.	19	) area=	0.989338633
          if (jElemUid .eq.	20	) area=	2.907061835
          if (jElemUid .eq.	21	) area=	0.435474748
          if (jElemUid .eq.	22	) area=	2.515573778
          if (jElemUid .eq.	23	) area=	2.056447427
          if (jElemUid .eq.	24	) area=	2.206553689
          if (jElemUid .eq.	25	) area=	0.156694711
          if (jElemUid .eq.	26	) area=	0.43201801
          if (jElemUid .eq.	27	) area=	1.15616832
          if (jElemUid .eq.	28	) area=	2.363435102
          if (jElemUid .eq.	29	) area=	2.142076181
          if (jElemUid .eq.	30	) area=	0.497794729
          if (jElemUid .eq.	31	) area=	2.819998146
          if (jElemUid .eq.	32	) area=	1.607307658
          if (jElemUid .eq.	33	) area=	2.084356873
          if (jElemUid .eq.	34	) area=	0.453102469
          if (jElemUid .eq.	35	) area=	1.485164526
          if (jElemUid .eq.	36	) area=	1.692151071
          if (jElemUid .eq.	37	) area=	2.140587836
          if (jElemUid .eq.	38	) area=	0.065668384
          if (jElemUid .eq.	39	) area=	2.013406651
          if (jElemUid .eq.	40	) area=	1.540836346
          if (jElemUid .eq.	41	) area=	0.147106337
          if (jElemUid .eq.	42	) area=	2.332901531
          if (jElemUid .eq.	43	) area=	2.038757805
          if (jElemUid .eq.	44	) area=	0.663706122
          if (jElemUid .eq.	45	) area=	0.673327103
          if (jElemUid .eq.	46	) area=	0.965891064
          if (jElemUid .eq.	47	) area=	2.373430459
          if (jElemUid .eq.	48	) area=	2.022635323
          if (jElemUid .eq.	49	) area=	0.943132455
          if (jElemUid .eq.	50	) area=	1.691598548
          if (jElemUid .eq.	51	) area=	1.005856147
          if (jElemUid .eq.	52	) area=	2.086246448
          if (jElemUid .eq.	53	) area=	0.582942114
          if (jElemUid .eq.	54	) area=	1.63905418
          if (jElemUid .eq.	55	) area=	0.846479045
          if (jElemUid .eq.	56	) area=	1.666357891
          if (jElemUid .eq.	57	) area=	1.860040776
          if (jElemUid .eq.	58	) area=	1.794211242
          if (jElemUid .eq.	59	) area=	2.203637546
          if (jElemUid .eq.	60	) area=	1.171571123
          if (jElemUid .eq.	61	) area=	0.166892978
          if (jElemUid .eq.	62	) area=	1.322628151
          if (jElemUid .eq.	63	) area=	2.097001828
          if (jElemUid .eq.	64	) area=	1.873410273
          if (jElemUid .eq.	65	) area=	1.820780178
          if (jElemUid .eq.	66	) area=	1.713734509
          if (jElemUid .eq.	67	) area=	0.522541117
          if (jElemUid .eq.	68	) area=	2.067543521
          if (jElemUid .eq.	69	) area=	1.631727842
          if (jElemUid .eq.	70	) area=	1.980820373
          if (jElemUid .eq.	71	) area=	1.169794323
          if (jElemUid .eq.	72	) area=	1.078360144
          if (jElemUid .eq.	73	) area=	2.205419421
          if (jElemUid .eq.	74	) area=	2.415643946
          if (jElemUid .eq.	75	) area=	0.749134676
          if (jElemUid .eq.	76	) area=	1.334832748
          if (jElemUid .eq.	77	) area=	2.136840987
          if (jElemUid .eq.	78	) area=	0.914510543
          if (jElemUid .eq.	79	) area=	0.987013572
          if (jElemUid .eq.	80	) area=	2.122290712
          if (jElemUid .eq.	81	) area=	2.246958478
          if (jElemUid .eq.	82	) area=	2.155549477
          if (jElemUid .eq.	83	) area=	1.348869546
          if (jElemUid .eq.	84	) area=	1.9876796
          if (jElemUid .eq.	85	) area=	1.620537851
          if (jElemUid .eq.	86	) area=	2.579195236
          if (jElemUid .eq.	87	) area=	2.905408566
          if (jElemUid .eq.	88	) area=	0.557701209
          if (jElemUid .eq.	89	) area=	1.780954346
          if (jElemUid .eq.	90	) area=	2.316782029
          if (jElemUid .eq.	91	) area=	0.138308137
          if (jElemUid .eq.	92	) area=	0.327047302
          if (jElemUid .eq.	93	) area=	2.398525884
          if (jElemUid .eq.	94	) area=	1.521952984
          if (jElemUid .eq.	95	) area=	2.501079739
          if (jElemUid .eq.	96	) area=	2.445128358
          if (jElemUid .eq.	97	) area=	1.991753945
          if (jElemUid .eq.	98	) area=	2.228844682
          if (jElemUid .eq.	99	) area=	1.834616107
          if (jElemUid .eq.	100	) area=	1.405780402
          if (jElemUid .eq.	101	) area=	0.479831495
          if (jElemUid .eq.	102	) area=	2.380065511
          if (jElemUid .eq.	103	) area=	1.964064629
          if (jElemUid .eq.	104	) area=	2.049547415
          if (jElemUid .eq.	105	) area=	1.856034962
          if (jElemUid .eq.	106	) area=	0.368584469
          if (jElemUid .eq.	107	) area=	0.996617794
          if (jElemUid .eq.	108	) area=	1.071443193
          if (jElemUid .eq.	109	) area=	1.469240853
          if (jElemUid .eq.	110	) area=	2.291013357
          if (jElemUid .eq.	111	) area=	2.162502853
          if (jElemUid .eq.	112	) area=	0.823536927
          if (jElemUid .eq.	113	) area=	1.659657843
          if (jElemUid .eq.	114	) area=	2.325904187
          if (jElemUid .eq.	115	) area=	1.814499501
          if (jElemUid .eq.	116	) area=	0.862200399
          if (jElemUid .eq.	117	) area=	0.875104896
          if (jElemUid .eq.	118	) area=	1.789265722
          if (jElemUid .eq.	119	) area=	2.078550969
          if (jElemUid .eq.	120	) area=	1.861921424
          if (jElemUid .eq.	121	) area=	2.149892215
          if (jElemUid .eq.	122	) area=	1.761052691
          if (jElemUid .eq.	123	) area=	1.216269081
          if (jElemUid .eq.	124	) area=	1.274217492
          if (jElemUid .eq.	125	) area=	1.177858499
          if (jElemUid .eq.	126	) area=	2.85063735
          if (jElemUid .eq.	127	) area=	2.087767593
          if (jElemUid .eq.	128	) area=	0.984081937
          if (jElemUid .eq.	129	) area=	1.37400615
          if (jElemUid .eq.	130	) area=	0.001981754
          if (jElemUid .eq.	131	) area=	0.328063399
          if (jElemUid .eq.	132	) area=	1.657322061
          if (jElemUid .eq.	133	) area=	2.546745281
          if (jElemUid .eq.	134	) area=	1.932482551
          if (jElemUid .eq.	135	) area=	2.587362968
          if (jElemUid .eq.	136	) area=	1.831581718
          if (jElemUid .eq.	137	) area=	2.303022152
          if (jElemUid .eq.	138	) area=	1.49342593
          if (jElemUid .eq.	139	) area=	1.670986095
          if (jElemUid .eq.	140	) area=	0.86196464
          if (jElemUid .eq.	141	) area=	1.943894041
          if (jElemUid .eq.	142	) area=	2.118011934
          if (jElemUid .eq.	143	) area=	1.035888818
          if (jElemUid .eq.	144	) area=	2.278501852
          if (jElemUid .eq.	145	) area=	2.737740397
          if (jElemUid .eq.	146	) area=	0.779082132
          if (jElemUid .eq.	147	) area=	1.946807816
          if (jElemUid .eq.	148	) area=	1.665716916
          if (jElemUid .eq.	149	) area=	1.066400475
          if (jElemUid .eq.	150	) area=	1.496853978
          if (jElemUid .eq.	151	) area=	1.713167708
          if (jElemUid .eq.	152	) area=	1.292006515
          if (jElemUid .eq.	153	) area=	1.596965827
          if (jElemUid .eq.	154	) area=	0.210325447
          if (jElemUid .eq.	155	) area=	2.946547478
          if (jElemUid .eq.	156	) area=	0.221120758
          if (jElemUid .eq.	157	) area=	2.942411777
          if (jElemUid .eq.	158	) area=	3.071852461
          if (jElemUid .eq.	159	) area=	1.030934593
          if (jElemUid .eq.	160	) area=	2.449367722
          if (jElemUid .eq.	161	) area=	2.110838833
          if (jElemUid .eq.	162	) area=	1.798615053
          if (jElemUid .eq.	163	) area=	2.439286322
          if (jElemUid .eq.	164	) area=	1.0173886
          if (jElemUid .eq.	165	) area=	0.165581952
          if (jElemUid .eq.	166	) area=	1.964173678
          if (jElemUid .eq.	167	) area=	1.796688424
          if (jElemUid .eq.	168	) area=	2.02980588
          if (jElemUid .eq.	169	) area=	2.669615561
          if (jElemUid .eq.	170	) area=	0.999082789
          if (jElemUid .eq.	171	) area=	1.017562479
          if (jElemUid .eq.	172	) area=	1.594793346
          if (jElemUid .eq.	173	) area=	1.960804943
          if (jElemUid .eq.	174	) area=	2.212927811
          if (jElemUid .eq.	175	) area=	0.899781376
          if (jElemUid .eq.	176	) area=	1.224560176
          if (jElemUid .eq.	177	) area=	1.98519297
          if (jElemUid .eq.	178	) area=	1.523671061
          if (jElemUid .eq.	179	) area=	2.163055538
          if (jElemUid .eq.	180	) area=	0.500240372
          if (jElemUid .eq.	181	) area=	1.230092686
          if (jElemUid .eq.	182	) area=	2.074965418
          if (jElemUid .eq.	183	) area=	1.448147624
          if (jElemUid .eq.	184	) area=	1.943725168
          if (jElemUid .eq.	185	) area=	1.552234061
          if (jElemUid .eq.	186	) area=	1.862148241
          if (jElemUid .eq.	187	) area=	1.608058641
          if (jElemUid .eq.	188	) area=	1.984259549
          if (jElemUid .eq.	189	) area=	1.982392063
          if (jElemUid .eq.	190	) area=	1.049193745
          if (jElemUid .eq.	191	) area=	1.305986258
          if (jElemUid .eq.	192	) area=	1.323325014
          if (jElemUid .eq.	193	) area=	1.703154093
          if (jElemUid .eq.	194	) area=	1.188891704
          if (jElemUid .eq.	195	) area=	1.908963389
          if (jElemUid .eq.	196	) area=	1.582404293
          if (jElemUid .eq.	197	) area=	0.369377992
          if (jElemUid .eq.	198	) area=	2.51730228
          if (jElemUid .eq.	199	) area=	2.632091973
          if (jElemUid .eq.	200	) area=	1.771002886
          if (jElemUid .eq.	201	) area=	1.338205754
          if (jElemUid .eq.	202	) area=	1.782187561
          if (jElemUid .eq.	203	) area=	1.554626473
          if (jElemUid .eq.	204	) area=	2.678511737
          if (jElemUid .eq.	205	) area=	1.164023528
          if (jElemUid .eq.	206	) area=	0.631661456
          if (jElemUid .eq.	207	) area=	2.37982426
          if (jElemUid .eq.	208	) area=	1.652374243
          if (jElemUid .eq.	209	) area=	0.740933268
          if (jElemUid .eq.	210	) area=	1.001418869
          if (jElemUid .eq.	211	) area=	2.765816837
          if (jElemUid .eq.	212	) area=	2.367236624
          if (jElemUid .eq.	213	) area=	1.227178337
          if (jElemUid .eq.	214	) area=	2.06344996
          if (jElemUid .eq.	215	) area=	2.196788013
          if (jElemUid .eq.	216	) area=	1.436250902
          if (jElemUid .eq.	217	) area=	0.889720224
          if (jElemUid .eq.	218	) area=	2.001811447
          if (jElemUid .eq.	219	) area=	0.045938993
          if (jElemUid .eq.	220	) area=	2.494614302
          if (jElemUid .eq.	221	) area=	1.763875228
          if (jElemUid .eq.	222	) area=	0.738400164
          if (jElemUid .eq.	223	) area=	1.932823302
          if (jElemUid .eq.	224	) area=	2.116599614
          if (jElemUid .eq.	225	) area=	1.339970106
          if (jElemUid .eq.	226	) area=	0.070059347
          if (jElemUid .eq.	227	) area=	1.755530477
          if (jElemUid .eq.	228	) area=	1.774289457
          if (jElemUid .eq.	229	) area=	1.286200593
          if (jElemUid .eq.	230	) area=	0.763395362
          if (jElemUid .eq.	231	) area=	2.603122422
          if (jElemUid .eq.	232	) area=	0.245258077
          if (jElemUid .eq.	233	) area=	2.128012408
          if (jElemUid .eq.	234	) area=	1.667284302
          if (jElemUid .eq.	235	) area=	0.665374663
          if (jElemUid .eq.	236	) area=	1.11857718
          if (jElemUid .eq.	237	) area=	2.393678023
          if (jElemUid .eq.	238	) area=	1.943130209
          if (jElemUid .eq.	239	) area=	2.107346848
          if (jElemUid .eq.	240	) area=	2.359454804
          if (jElemUid .eq.	241	) area=	2.082254857
          if (jElemUid .eq.	242	) area=	1.567777889
          if (jElemUid .eq.	243	) area=	2.039648103
          if (jElemUid .eq.	244	) area=	3.285756095
          if (jElemUid .eq.	245	) area=	1.377034147
          if (jElemUid .eq.	246	) area=	2.462414712
          if (jElemUid .eq.	247	) area=	1.998675304
          if (jElemUid .eq.	248	) area=	1.894543183
          if (jElemUid .eq.	249	) area=	0.255111243
          if (jElemUid .eq.	250	) area=	2.385865083
          if (jElemUid .eq.	251	) area=	0.360703479
          if (jElemUid .eq.	252	) area=	3.22690358
          if (jElemUid .eq.	253	) area=	1.971590883
          if (jElemUid .eq.	254	) area=	2.246291401
          if (jElemUid .eq.	255	) area=	2.924567708
          if (jElemUid .eq.	256	) area=	1.336403088
          if (jElemUid .eq.	257	) area=	1.784033823
          if (jElemUid .eq.	258	) area=	2.325664302
          if (jElemUid .eq.	259	) area=	1.413302919
          if (jElemUid .eq.	260	) area=	2.237851445
          if (jElemUid .eq.	261	) area=	0.716731707
          if (jElemUid .eq.	262	) area=	0.607229524
          if (jElemUid .eq.	263	) area=	1.662097792
          if (jElemUid .eq.	264	) area=	2.103500433
          if (jElemUid .eq.	265	) area=	1.922892249
          if (jElemUid .eq.	266	) area=	2.315106468
          if (jElemUid .eq.	267	) area=	1.215660323
          if (jElemUid .eq.	268	) area=	3.33302207
          if (jElemUid .eq.	269	) area=	1.094288313
          if (jElemUid .eq.	270	) area=	0.986535078
          if (jElemUid .eq.	271	) area=	1.448039112
          if (jElemUid .eq.	272	) area=	1.576621944
          if (jElemUid .eq.	273	) area=	1.635572614
          if (jElemUid .eq.	274	) area=	1.482784171
          if (jElemUid .eq.	275	) area=	1.269578929
          if (jElemUid .eq.	276	) area=	2.394155619
          if (jElemUid .eq.	277	) area=	2.934309552
          if (jElemUid .eq.	278	) area=	0.374479976
          if (jElemUid .eq.	279	) area=	1.058308193
          if (jElemUid .eq.	280	) area=	0.146006181
          if (jElemUid .eq.	281	) area=	2.8660233
          if (jElemUid .eq.	282	) area=	2.766318178
          if (jElemUid .eq.	283	) area=	2.237134642
          if (jElemUid .eq.	284	) area=	2.589914443
          if (jElemUid .eq.	285	) area=	2.380251073
          if (jElemUid .eq.	286	) area=	1.30167132
          if (jElemUid .eq.	287	) area=	1.061020717
          if (jElemUid .eq.	288	) area=	0.801778366
          if (jElemUid .eq.	289	) area=	1.765228237
          if (jElemUid .eq.	290	) area=	2.117632495
          if (jElemUid .eq.	291	) area=	2.325729063
          if (jElemUid .eq.	292	) area=	1.406992202
          if (jElemUid .eq.	293	) area=	1.006163868
          if (jElemUid .eq.	294	) area=	0.069089718
          if (jElemUid .eq.	295	) area=	2.046801039
          if (jElemUid .eq.	296	) area=	2.061174438
          if (jElemUid .eq.	297	) area=	2.140106375
          if (jElemUid .eq.	298	) area=	0.86058806
          if (jElemUid .eq.	299	) area=	1.909801991
          if (jElemUid .eq.	300	) area=	1.942284226
          if (jElemUid .eq.	301	) area=	1.35964735
          if (jElemUid .eq.	302	) area=	2.356591929
          if (jElemUid .eq.	303	) area=	1.046291453
          if (jElemUid .eq.	304	) area=	1.512071999
          if (jElemUid .eq.	305	) area=	2.202999117
          if (jElemUid .eq.	306	) area=	1.124666439
          if (jElemUid .eq.	307	) area=	2.638113419
          if (jElemUid .eq.	308	) area=	0.77471928
          if (jElemUid .eq.	309	) area=	1.829766018
          if (jElemUid .eq.	310	) area=	1.21708484
          if (jElemUid .eq.	311	) area=	1.110115759
          if (jElemUid .eq.	312	) area=	1.920675849
          if (jElemUid .eq.	313	) area=	2.277969943
          if (jElemUid .eq.	314	) area=	2.09938103
          if (jElemUid .eq.	315	) area=	0.998109352
          if (jElemUid .eq.	316	) area=	2.016753655
          if (jElemUid .eq.	317	) area=	2.865951252
          if (jElemUid .eq.	318	) area=	1.190684405
          if (jElemUid .eq.	319	) area=	1.447308618
          if (jElemUid .eq.	320	) area=	1.438486204
          if (jElemUid .eq.	321	) area=	2.263141311
          if (jElemUid .eq.	322	) area=	1.517577588
          if (jElemUid .eq.	323	) area=	1.596819682
          if (jElemUid .eq.	324	) area=	0.994105427
          if (jElemUid .eq.	325	) area=	1.06234516
          if (jElemUid .eq.	326	) area=	2.181677833
          if (jElemUid .eq.	327	) area=	1.643431061
          if (jElemUid .eq.	328	) area=	1.535280597
          if (jElemUid .eq.	329	) area=	1.62417649
          if (jElemUid .eq.	330	) area=	2.037481189
          if (jElemUid .eq.	331	) area=	2.189462456
          if (jElemUid .eq.	332	) area=	1.923459097
          if (jElemUid .eq.	333	) area=	1.836051845
          if (jElemUid .eq.	334	) area=	1.707862311
          if (jElemUid .eq.	335	) area=	2.610404717
          if (jElemUid .eq.	336	) area=	2.362924499
          if (jElemUid .eq.	337	) area=	1.775247887
          if (jElemUid .eq.	338	) area=	1.487942451
          if (jElemUid .eq.	339	) area=	2.220093837
          if (jElemUid .eq.	340	) area=	1.227812034
          if (jElemUid .eq.	341	) area=	2.798248859
          if (jElemUid .eq.	342	) area=	1.322880693
          if (jElemUid .eq.	343	) area=	1.838567671
          if (jElemUid .eq.	344	) area=	1.94513405
          if (jElemUid .eq.	345	) area=	2.322471218
          if (jElemUid .eq.	346	) area=	1.269337164
          if (jElemUid .eq.	347	) area=	1.699117044
          if (jElemUid .eq.	348	) area=	2.207376219
          if (jElemUid .eq.	349	) area=	2.183201432
          if (jElemUid .eq.	350	) area=	0.063433466
          if (jElemUid .eq.	351	) area=	1.948116922
          if (jElemUid .eq.	352	) area=	2.21534515
          if (jElemUid .eq.	353	) area=	0.597697513
          if (jElemUid .eq.	354	) area=	2.955070645
          if (jElemUid .eq.	355	) area=	1.424155057
          if (jElemUid .eq.	356	) area=	2.014390453
          if (jElemUid .eq.	357	) area=	2.257891355
          if (jElemUid .eq.	358	) area=	1.692561338
          if (jElemUid .eq.	359	) area=	1.287229251
          if (jElemUid .eq.	360	) area=	1.511483088
          if (jElemUid .eq.	361	) area=	2.147869872
          if (jElemUid .eq.	362	) area=	2.336213356
          if (jElemUid .eq.	363	) area=	1.934063807
          if (jElemUid .eq.	364	) area=	1.294177357
          if (jElemUid .eq.	365	) area=	2.149390792
          if (jElemUid .eq.	366	) area=	1.732488235
          if (jElemUid .eq.	367	) area=	1.625133906
          if (jElemUid .eq.	368	) area=	1.266112392
          if (jElemUid .eq.	369	) area=	1.861909268
          if (jElemUid .eq.	370	) area=	1.78836535
          if (jElemUid .eq.	371	) area=	1.079252612
          if (jElemUid .eq.	372	) area=	1.88410436
          if (jElemUid .eq.	373	) area=	1.710349359
          if (jElemUid .eq.	374	) area=	1.472774132
          if (jElemUid .eq.	375	) area=	2.854564773
          if (jElemUid .eq.	376	) area=	1.75966006
          if (jElemUid .eq.	377	) area=	0.69958446
          if (jElemUid .eq.	378	) area=	0.23045644
          if (jElemUid .eq.	379	) area=	2.51015125
          if (jElemUid .eq.	380	) area=	0.828145727
          if (jElemUid .eq.	381	) area=	1.912001416
          if (jElemUid .eq.	382	) area=	1.799894632
          if (jElemUid .eq.	383	) area=	2.120284408
          if (jElemUid .eq.	384	) area=	1.050295159
          if (jElemUid .eq.	385	) area=	0.909723225
          if (jElemUid .eq.	386	) area=	0.745563757
          if (jElemUid .eq.	387	) area=	2.582848178
          if (jElemUid .eq.	388	) area=	1.99362005
          if (jElemUid .eq.	389	) area=	1.445674913
          if (jElemUid .eq.	390	) area=	2.526932816
          if (jElemUid .eq.	391	) area=	0.608186825
          if (jElemUid .eq.	392	) area=	2.643074205
          if (jElemUid .eq.	393	) area=	2.714219518
          if (jElemUid .eq.	394	) area=	1.494338426
          if (jElemUid .eq.	395	) area=	0.239161373
          if (jElemUid .eq.	396	) area=	2.824444017
          if (jElemUid .eq.	397	) area=	1.384312799
          if (jElemUid .eq.	398	) area=	2.086542847
          if (jElemUid .eq.	399	) area=	1.349254741
          if (jElemUid .eq.	400	) area=	1.477934285
          if (jElemUid .eq.	401	) area=	2.050533692
          if (jElemUid .eq.	402	) area=	2.418635176
          if (jElemUid .eq.	403	) area=	1.305547747
          if (jElemUid .eq.	404	) area=	1.827549357
          if (jElemUid .eq.	405	) area=	1.364906014
          if (jElemUid .eq.	406	) area=	0.427565127
          if (jElemUid .eq.	407	) area=	2.003785438
          if (jElemUid .eq.	408	) area=	2.357671269
          if (jElemUid .eq.	409	) area=	0.898156859
          if (jElemUid .eq.	410	) area=	0.975468499
          if (jElemUid .eq.	411	) area=	1.754774963
          if (jElemUid .eq.	412	) area=	1.898716159
          if (jElemUid .eq.	413	) area=	2.242522286
          if (jElemUid .eq.	414	) area=	0.705361278
          if (jElemUid .eq.	415	) area=	1.752817443
          if (jElemUid .eq.	416	) area=	2.250497556
          if (jElemUid .eq.	417	) area=	0.926745084
          if (jElemUid .eq.	418	) area=	0.60672722
          if (jElemUid .eq.	419	) area=	1.823766722
          if (jElemUid .eq.	420	) area=	2.409657814
          if (jElemUid .eq.	421	) area=	1.889870158
          if (jElemUid .eq.	422	) area=	2.458533576
          if (jElemUid .eq.	423	) area=	1.788924349
          if (jElemUid .eq.	424	) area=	1.823520328
          if (jElemUid .eq.	425	) area=	2.657854256
          if (jElemUid .eq.	426	) area=	2.302928483
          if (jElemUid .eq.	427	) area=	2.494715
          if (jElemUid .eq.	428	) area=	2.430910382
          if (jElemUid .eq.	429	) area=	0.960329605
          if (jElemUid .eq.	430	) area=	2.064080994
          if (jElemUid .eq.	431	) area=	0.71501105
          if (jElemUid .eq.	432	) area=	1.43132501
          if (jElemUid .eq.	433	) area=	1.739775773
          if (jElemUid .eq.	434	) area=	1.714779267
          if (jElemUid .eq.	435	) area=	1.72471497
          if (jElemUid .eq.	436	) area=	1.835948464
          if (jElemUid .eq.	437	) area=	2.071867053
          if (jElemUid .eq.	438	) area=	1.563156955
          if (jElemUid .eq.	439	) area=	1.899687988
          if (jElemUid .eq.	440	) area=	0.807662118
          if (jElemUid .eq.	441	) area=	1.706418015
          if (jElemUid .eq.	442	) area=	1.647392053
          if (jElemUid .eq.	443	) area=	1.45502739
          if (jElemUid .eq.	444	) area=	1.30298835
          if (jElemUid .eq.	445	) area=	2.452113609
          if (jElemUid .eq.	446	) area=	1.80598073
          if (jElemUid .eq.	447	) area=	1.209661734
          if (jElemUid .eq.	448	) area=	0.506010165
          if (jElemUid .eq.	449	) area=	1.880018855
          if (jElemUid .eq.	450	) area=	1.602954451
          if (jElemUid .eq.	451	) area=	1.877823922
          if (jElemUid .eq.	452	) area=	2.080597895
          if (jElemUid .eq.	453	) area=	2.389822585
          if (jElemUid .eq.	454	) area=	1.534939966
          if (jElemUid .eq.	455	) area=	1.149531983
          if (jElemUid .eq.	456	) area=	0.059544313
          if (jElemUid .eq.	457	) area=	2.548447682
          if (jElemUid .eq.	458	) area=	2.541370188
          if (jElemUid .eq.	459	) area=	1.969244444
          if (jElemUid .eq.	460	) area=	1.42010364
          if (jElemUid .eq.	461	) area=	2.079872074
          if (jElemUid .eq.	462	) area=	1.874112267
          if (jElemUid .eq.	463	) area=	2.762349351
          if (jElemUid .eq.	464	) area=	1.562588869
          if (jElemUid .eq.	465	) area=	1.743394055
          if (jElemUid .eq.	466	) area=	1.70295233
          if (jElemUid .eq.	467	) area=	2.425161125
          if (jElemUid .eq.	468	) area=	1.602618086
          if (jElemUid .eq.	469	) area=	0.609501523
          if (jElemUid .eq.	470	) area=	2.661253964
          if (jElemUid .eq.	471	) area=	1.978673629
          if (jElemUid .eq.	472	) area=	2.306675357
          if (jElemUid .eq.	473	) area=	1.000128579
          if (jElemUid .eq.	474	) area=	2.66472971
          if (jElemUid .eq.	475	) area=	1.233327292
          if (jElemUid .eq.	476	) area=	1.917465622
          if (jElemUid .eq.	477	) area=	1.834583149
          if (jElemUid .eq.	478	) area=	1.597854548
          if (jElemUid .eq.	479	) area=	1.978031765
          if (jElemUid .eq.	480	) area=	2.07586506
          if (jElemUid .eq.	481	) area=	2.859596166
          if (jElemUid .eq.	482	) area=	2.360842567
          if (jElemUid .eq.	483	) area=	0.083423271
          if (jElemUid .eq.	484	) area=	1.066772967
          if (jElemUid .eq.	485	) area=	1.219916905
          if (jElemUid .eq.	486	) area=	2.509727768
          if (jElemUid .eq.	487	) area=	2.77795342
          if (jElemUid .eq.	488	) area=	2.209061172
          if (jElemUid .eq.	489	) area=	1.812316655
          if (jElemUid .eq.	490	) area=	1.831690153
          if (jElemUid .eq.	491	) area=	2.086107402
          if (jElemUid .eq.	492	) area=	1.426485507
          if (jElemUid .eq.	493	) area=	0.86198685
          if (jElemUid .eq.	494	) area=	1.438444492
          if (jElemUid .eq.	495	) area=	2.537335291
          if (jElemUid .eq.	496	) area=	2.375981922
          if (jElemUid .eq.	497	) area=	1.492036913
          if (jElemUid .eq.	498	) area=	1.657691411
          if (jElemUid .eq.	499	) area=	1.19698999
          if (jElemUid .eq.	500	) area=	1.801137513
          if (jElemUid .eq.	501	) area=	1.627187484
          if (jElemUid .eq.	502	) area=	1.621134455
          if (jElemUid .eq.	503	) area=	0.977319116
          if (jElemUid .eq.	504	) area=	1.615247717
          if (jElemUid .eq.	505	) area=	2.48361493
          if (jElemUid .eq.	506	) area=	1.584905913
          if (jElemUid .eq.	507	) area=	2.295408142
          if (jElemUid .eq.	508	) area=	1.280990912
          if (jElemUid .eq.	509	) area=	1.006928878
          if (jElemUid .eq.	510	) area=	2.547721093
          if (jElemUid .eq.	511	) area=	2.121674147
          if (jElemUid .eq.	512	) area=	2.521474625
          if (jElemUid .eq.	513	) area=	0.855937488
          if (jElemUid .eq.	514	) area=	1.376921784
          if (jElemUid .eq.	515	) area=	1.816497066
          if (jElemUid .eq.	516	) area=	1.813071452
          if (jElemUid .eq.	517	) area=	1.375826403
          if (jElemUid .eq.	518	) area=	1.427401359
          if (jElemUid .eq.	519	) area=	1.71369869
          if (jElemUid .eq.	520	) area=	1.345299087
          if (jElemUid .eq.	521	) area=	2.275484693
          if (jElemUid .eq.	522	) area=	1.585613775
          if (jElemUid .eq.	523	) area=	2.257114376
          if (jElemUid .eq.	524	) area=	1.91133238
          if (jElemUid .eq.	525	) area=	2.254264315
          if (jElemUid .eq.	526	) area=	0.119706067
          if (jElemUid .eq.	527	) area=	1.582049255
          if (jElemUid .eq.	528	) area=	2.021800061
          if (jElemUid .eq.	529	) area=	1.385572473
          if (jElemUid .eq.	530	) area=	1.570021339
          if (jElemUid .eq.	531	) area=	2.399642672
          if (jElemUid .eq.	532	) area=	1.770532752
          if (jElemUid .eq.	533	) area=	0.013321428
          if (jElemUid .eq.	534	) area=	1.585254992
          if (jElemUid .eq.	535	) area=	1.366431556
          if (jElemUid .eq.	536	) area=	1.728310133
          if (jElemUid .eq.	537	) area=	1.949739943
          if (jElemUid .eq.	538	) area=	0.030033325
          if (jElemUid .eq.	539	) area=	3.376102848
          if (jElemUid .eq.	540	) area=	0.115196579
          if (jElemUid .eq.	541	) area=	2.95946101
          if (jElemUid .eq.	542	) area=	1.129559434
          if (jElemUid .eq.	543	) area=	2.707057118
          if (jElemUid .eq.	544	) area=	0.769513458
          if (jElemUid .eq.	545	) area=	2.870024635
          if (jElemUid .eq.	546	) area=	0.736721457
          if (jElemUid .eq.	547	) area=	2.156016269
          if (jElemUid .eq.	548	) area=	2.192720947
          if (jElemUid .eq.	549	) area=	0.870545645
          if (jElemUid .eq.	550	) area=	1.880281521
          if (jElemUid .eq.	551	) area=	1.997386161
          if (jElemUid .eq.	552	) area=	2.717137809
          if (jElemUid .eq.	553	) area=	1.517729463
          if (jElemUid .eq.	554	) area=	0.901560675
          if (jElemUid .eq.	555	) area=	2.254618368
          if (jElemUid .eq.	556	) area=	1.04378256
          if (jElemUid .eq.	557	) area=	1.824096511
          if (jElemUid .eq.	558	) area=	1.381733391
          if (jElemUid .eq.	559	) area=	2.031291924
          if (jElemUid .eq.	560	) area=	1.504753261
          if (jElemUid .eq.	561	) area=	1.043354882
          if (jElemUid .eq.	562	) area=	0.947961303
          if (jElemUid .eq.	563	) area=	1.044595566
          if (jElemUid .eq.	564	) area=	0.769125443
          if (jElemUid .eq.	565	) area=	3.192556668
          if (jElemUid .eq.	566	) area=	2.292187843
          if (jElemUid .eq.	567	) area=	1.788809965
          if (jElemUid .eq.	568	) area=	0.731718088
          if (jElemUid .eq.	569	) area=	0.975424987
          if (jElemUid .eq.	570	) area=	1.283489409
          if (jElemUid .eq.	571	) area=	2.739359561
          if (jElemUid .eq.	572	) area=	1.489883026
          if (jElemUid .eq.	573	) area=	0.174865207
          if (jElemUid .eq.	574	) area=	2.159657422
          if (jElemUid .eq.	575	) area=	2.140154791
          if (jElemUid .eq.	576	) area=	1.852698682
          if (jElemUid .eq.	577	) area=	1.489126517
          if (jElemUid .eq.	578	) area=	2.238268185
          if (jElemUid .eq.	579	) area=	0.324621276
          if (jElemUid .eq.	580	) area=	1.843560403
          if (jElemUid .eq.	581	) area=	0.355914162
          if (jElemUid .eq.	582	) area=	0.260569191
          if (jElemUid .eq.	583	) area=	1.786079006
          if (jElemUid .eq.	584	) area=	2.229597261
          if (jElemUid .eq.	585	) area=	2.23237153
          if (jElemUid .eq.	586	) area=	2.441739712
          if (jElemUid .eq.	587	) area=	0.829903231
          if (jElemUid .eq.	588	) area=	2.201036618
          if (jElemUid .eq.	589	) area=	1.699723651
          if (jElemUid .eq.	590	) area=	2.665678877
          if (jElemUid .eq.	591	) area=	2.703222849
          if (jElemUid .eq.	592	) area=	1.921549468
          if (jElemUid .eq.	593	) area=	2.048779423
          if (jElemUid .eq.	594	) area=	1.706876209
          if (jElemUid .eq.	595	) area=	2.43130866
          if (jElemUid .eq.	596	) area=	1.783541348
          if (jElemUid .eq.	597	) area=	2.241993189
          if (jElemUid .eq.	598	) area=	0.870175913
          if (jElemUid .eq.	599	) area=	1.517824102
          if (jElemUid .eq.	600	) area=	2.361267822
          if (jElemUid .eq.	601	) area=	0.48429055
          if (jElemUid .eq.	602	) area=	1.42983337
          if (jElemUid .eq.	603	) area=	2.767307429
          if (jElemUid .eq.	604	) area=	1.824772696
          if (jElemUid .eq.	605	) area=	1.030547449
          if (jElemUid .eq.	606	) area=	1.236081392
          if (jElemUid .eq.	607	) area=	1.663218271
          if (jElemUid .eq.	608	) area=	2.369425078
          if (jElemUid .eq.	609	) area=	1.36481774
          if (jElemUid .eq.	610	) area=	2.446203427
          if (jElemUid .eq.	611	) area=	1.629096028
          if (jElemUid .eq.	612	) area=	3.22701118
          if (jElemUid .eq.	613	) area=	0.410080786
          if (jElemUid .eq.	614	) area=	0.664507227
          if (jElemUid .eq.	615	) area=	0.350830719
          if (jElemUid .eq.	616	) area=	2.315106243
          if (jElemUid .eq.	617	) area=	0.174126177
          if (jElemUid .eq.	618	) area=	2.154871438
          if (jElemUid .eq.	619	) area=	2.364180904
          if (jElemUid .eq.	620	) area=	1.813186637
          if (jElemUid .eq.	621	) area=	1.666561549
          if (jElemUid .eq.	622	) area=	1.641989863
          if (jElemUid .eq.	623	) area=	1.854761893
          if (jElemUid .eq.	624	) area=	2.627427882
          if (jElemUid .eq.	625	) area=	1.133889042
          if (jElemUid .eq.	626	) area=	1.805353475
          if (jElemUid .eq.	627	) area=	0.894680664
          if (jElemUid .eq.	628	) area=	1.737603662
          if (jElemUid .eq.	629	) area=	1.489450188
          if (jElemUid .eq.	630	) area=	1.52248662
          if (jElemUid .eq.	631	) area=	2.454358289
          if (jElemUid .eq.	632	) area=	0.523701459
          if (jElemUid .eq.	633	) area=	2.615395187
          if (jElemUid .eq.	634	) area=	2.082340144
          if (jElemUid .eq.	635	) area=	1.014352007
          if (jElemUid .eq.	636	) area=	1.573762109
          if (jElemUid .eq.	637	) area=	2.397495252
          if (jElemUid .eq.	638	) area=	2.174301192
          if (jElemUid .eq.	639	) area=	1.197852707
          if (jElemUid .eq.	640	) area=	2.164811238
          if (jElemUid .eq.	641	) area=	2.036813324
          if (jElemUid .eq.	642	) area=	1.607400473
          if (jElemUid .eq.	643	) area=	1.948281371
          if (jElemUid .eq.	644	) area=	1.463807455
          if (jElemUid .eq.	645	) area=	1.686037173
          if (jElemUid .eq.	646	) area=	1.732644322
          if (jElemUid .eq.	647	) area=	2.472664006
          if (jElemUid .eq.	648	) area=	2.692775032
          if (jElemUid .eq.	649	) area=	1.51370733
          if (jElemUid .eq.	650	) area=	0.605462665
          if (jElemUid .eq.	651	) area=	1.761482306
          if (jElemUid .eq.	652	) area=	2.321927458
          if (jElemUid .eq.	653	) area=	2.27993898
          if (jElemUid .eq.	654	) area=	1.063297415
          if (jElemUid .eq.	655	) area=	2.144710974
          if (jElemUid .eq.	656	) area=	0.710790542
          if (jElemUid .eq.	657	) area=	2.012958993
          if (jElemUid .eq.	658	) area=	2.82480324
          if (jElemUid .eq.	659	) area=	1.978337517
          if (jElemUid .eq.	660	) area=	1.353711891
          if (jElemUid .eq.	661	) area=	0.756939923
          if (jElemUid .eq.	662	) area=	1.651258903
          if (jElemUid .eq.	663	) area=	2.348824188
          if (jElemUid .eq.	664	) area=	1.126563702
          if (jElemUid .eq.	665	) area=	1.749134406
          if (jElemUid .eq.	666	) area=	2.61046629
          if (jElemUid .eq.	667	) area=	1.134284737
          if (jElemUid .eq.	668	) area=	2.738135542
          if (jElemUid .eq.	669	) area=	0.114348196
          if (jElemUid .eq.	670	) area=	1.880859752
          if (jElemUid .eq.	671	) area=	0.534964344
          if (jElemUid .eq.	672	) area=	0.469486682
          if (jElemUid .eq.	673	) area=	1.861948938
          if (jElemUid .eq.	674	) area=	2.05031659
          if (jElemUid .eq.	675	) area=	1.009934294
          if (jElemUid .eq.	676	) area=	0.720977611
          if (jElemUid .eq.	677	) area=	1.301617862
          if (jElemUid .eq.	678	) area=	2.659719428
          if (jElemUid .eq.	679	) area=	1.476925917
          if (jElemUid .eq.	680	) area=	1.120273873
          if (jElemUid .eq.	681	) area=	1.832415411
          if (jElemUid .eq.	682	) area=	1.338180337
          if (jElemUid .eq.	683	) area=	2.437743459
          if (jElemUid .eq.	684	) area=	1.405585829
          if (jElemUid .eq.	685	) area=	0.861230807
          if (jElemUid .eq.	686	) area=	1.524146068
          if (jElemUid .eq.	687	) area=	0.680565909
          if (jElemUid .eq.	688	) area=	1.82443336
          if (jElemUid .eq.	689	) area=	2.187029759
          if (jElemUid .eq.	690	) area=	1.329983575
          if (jElemUid .eq.	691	) area=	2.221342172
          if (jElemUid .eq.	692	) area=	0.367890989
          if (jElemUid .eq.	693	) area=	2.175951365
          if (jElemUid .eq.	694	) area=	0.810516679
          if (jElemUid .eq.	695	) area=	1.778107727
          if (jElemUid .eq.	696	) area=	2.27395742
          if (jElemUid .eq.	697	) area=	2.268463658
          if (jElemUid .eq.	698	) area=	1.788013442
          if (jElemUid .eq.	699	) area=	1.197333431
          if (jElemUid .eq.	700	) area=	2.353045682
          if (jElemUid .eq.	701	) area=	2.396893007
          if (jElemUid .eq.	702	) area=	0.682544968
          if (jElemUid .eq.	703	) area=	0.92236134
          if (jElemUid .eq.	704	) area=	2.528334794
          if (jElemUid .eq.	705	) area=	2.260357107
          if (jElemUid .eq.	706	) area=	1.441266314
          if (jElemUid .eq.	707	) area=	1.119386286
          if (jElemUid .eq.	708	) area=	1.647515338
          if (jElemUid .eq.	709	) area=	1.157898038
          if (jElemUid .eq.	710	) area=	1.77148087
          if (jElemUid .eq.	711	) area=	1.389202295
          if (jElemUid .eq.	712	) area=	2.124156503
          if (jElemUid .eq.	713	) area=	1.921977332
          if (jElemUid .eq.	714	) area=	1.813103046
          if (jElemUid .eq.	715	) area=	2.533542269
          if (jElemUid .eq.	716	) area=	1.664038269
          if (jElemUid .eq.	717	) area=	0.412452592
          if (jElemUid .eq.	718	) area=	2.001325802
          if (jElemUid .eq.	719	) area=	1.759330989
          if (jElemUid .eq.	720	) area=	0.399215636
          if (jElemUid .eq.	721	) area=	1.451004993
          if (jElemUid .eq.	722	) area=	1.639638204
          if (jElemUid .eq.	723	) area=	1.838819836
          if (jElemUid .eq.	724	) area=	1.902742945
          if (jElemUid .eq.	725	) area=	1.722900244
          if (jElemUid .eq.	726	) area=	1.667802678
          if (jElemUid .eq.	727	) area=	1.638091671
          if (jElemUid .eq.	728	) area=	1.578177777
          if (jElemUid .eq.	729	) area=	1.621235326
          if (jElemUid .eq.	730	) area=	0.804838794
          if (jElemUid .eq.	731	) area=	1.839263725
          if (jElemUid .eq.	732	) area=	1.987834801
          if (jElemUid .eq.	733	) area=	0.894657032
          if (jElemUid .eq.	734	) area=	1.748718236
          if (jElemUid .eq.	735	) area=	1.957213716
          if (jElemUid .eq.	736	) area=	2.042619242
          if (jElemUid .eq.	737	) area=	1.977790513
          if (jElemUid .eq.	738	) area=	2.13276335
          if (jElemUid .eq.	739	) area=	2.117050283
          if (jElemUid .eq.	740	) area=	0.747966949
          if (jElemUid .eq.	741	) area=	2.428335974
          if (jElemUid .eq.	742	) area=	1.163839675
          if (jElemUid .eq.	743	) area=	1.469995121
          if (jElemUid .eq.	744	) area=	1.67133298
          if (jElemUid .eq.	745	) area=	1.764505308
          if (jElemUid .eq.	746	) area=	1.963146221
          if (jElemUid .eq.	747	) area=	1.550450988
          if (jElemUid .eq.	748	) area=	0.291081247
          if (jElemUid .eq.	749	) area=	1.504691653
          if (jElemUid .eq.	750	) area=	2.171360091
          if (jElemUid .eq.	751	) area=	2.177293097
          if (jElemUid .eq.	752	) area=	1.068900319
          if (jElemUid .eq.	753	) area=	0.25487374
          if (jElemUid .eq.	754	) area=	2.242400624
          if (jElemUid .eq.	755	) area=	1.297495529
          if (jElemUid .eq.	756	) area=	0.710558328
          if (jElemUid .eq.	757	) area=	2.55015749
          if (jElemUid .eq.	758	) area=	2.480940211
          if (jElemUid .eq.	759	) area=	0.285599648
          if (jElemUid .eq.	760	) area=	2.08286469
          if (jElemUid .eq.	761	) area=	1.379383106
          if (jElemUid .eq.	762	) area=	2.288184129
          if (jElemUid .eq.	763	) area=	1.708941993
          if (jElemUid .eq.	764	) area=	0.008588589
          if (jElemUid .eq.	765	) area=	2.494547429
          if (jElemUid .eq.	766	) area=	1.005909951
          if (jElemUid .eq.	767	) area=	2.543037127
          if (jElemUid .eq.	768	) area=	0.045259136
          if (jElemUid .eq.	769	) area=	1.918071302
          if (jElemUid .eq.	770	) area=	1.805132478
          if (jElemUid .eq.	771	) area=	1.942134995
          if (jElemUid .eq.	772	) area=	2.441987337
          if (jElemUid .eq.	773	) area=	1.358397794
          if (jElemUid .eq.	774	) area=	0.783885491
          if (jElemUid .eq.	775	) area=	1.384102338
          if (jElemUid .eq.	776	) area=	1.764268571
          if (jElemUid .eq.	777	) area=	1.488774772
          if (jElemUid .eq.	778	) area=	2.751402416
          if (jElemUid .eq.	779	) area=	2.300713725
          if (jElemUid .eq.	780	) area=	0.666920582
          if (jElemUid .eq.	781	) area=	1.034740622
          if (jElemUid .eq.	782	) area=	2.314203068
          if (jElemUid .eq.	783	) area=	2.545336894
          if (jElemUid .eq.	784	) area=	1.995364135
          if (jElemUid .eq.	785	) area=	0.497200842
          if (jElemUid .eq.	786	) area=	0.260184058
          if (jElemUid .eq.	787	) area=	2.766006336
          if (jElemUid .eq.	788	) area=	1.48703524
          if (jElemUid .eq.	789	) area=	2.699414763
          if (jElemUid .eq.	790	) area=	1.362754504
          if (jElemUid .eq.	791	) area=	1.649078183
          if (jElemUid .eq.	792	) area=	2.14059418
          if (jElemUid .eq.	793	) area=	1.896061971
          if (jElemUid .eq.	794	) area=	1.335291496
          if (jElemUid .eq.	795	) area=	1.716438678
          if (jElemUid .eq.	796	) area=	0.629805576
          if (jElemUid .eq.	797	) area=	2.637193252
          if (jElemUid .eq.	798	) area=	2.54976515
          if (jElemUid .eq.	799	) area=	1.404806132
          if (jElemUid .eq.	800	) area=	1.776959967
          if (jElemUid .eq.	801	) area=	1.700518577
          if (jElemUid .eq.	802	) area=	2.923536889
          if (jElemUid .eq.	803	) area=	0.015321478
          if (jElemUid .eq.	804	) area=	1.426428767
          if (jElemUid .eq.	805	) area=	2.400792662
          if (jElemUid .eq.	806	) area=	1.895112742
          if (jElemUid .eq.	807	) area=	2.069830307
          if (jElemUid .eq.	808	) area=	1.598897508
          if (jElemUid .eq.	809	) area=	0.568683463
          if (jElemUid .eq.	810	) area=	0.777748131
          if (jElemUid .eq.	811	) area=	1.309276504
          if (jElemUid .eq.	812	) area=	2.331121194
          if (jElemUid .eq.	813	) area=	1.23862425
          if (jElemUid .eq.	814	) area=	2.080576743
          if (jElemUid .eq.	815	) area=	1.771015807
          if (jElemUid .eq.	816	) area=	1.624192363
          if (jElemUid .eq.	817	) area=	0.762578455
          if (jElemUid .eq.	818	) area=	1.852006003
          if (jElemUid .eq.	819	) area=	1.679032617
          if (jElemUid .eq.	820	) area=	1.082209744
          if (jElemUid .eq.	821	) area=	3.075994628
          if (jElemUid .eq.	822	) area=	2.234534767
          if (jElemUid .eq.	823	) area=	0.608780388
          if (jElemUid .eq.	824	) area=	0.966115244
          if (jElemUid .eq.	825	) area=	1.252888928
          if (jElemUid .eq.	826	) area=	1.290114916
          if (jElemUid .eq.	827	) area=	0.292325947
          if (jElemUid .eq.	828	) area=	2.393363442
          if (jElemUid .eq.	829	) area=	2.162224429
          if (jElemUid .eq.	830	) area=	0.227832085
          if (jElemUid .eq.	831	) area=	0.355931386
          if (jElemUid .eq.	832	) area=	0.804563543
          if (jElemUid .eq.	833	) area=	1.713526289
          if (jElemUid .eq.	834	) area=	2.443822916
          if (jElemUid .eq.	835	) area=	1.532921852
          if (jElemUid .eq.	836	) area=	1.165420224
          if (jElemUid .eq.	837	) area=	1.457919612
          if (jElemUid .eq.	838	) area=	0.281307587
          if (jElemUid .eq.	839	) area=	2.327174373
          if (jElemUid .eq.	840	) area=	1.492530945
          if (jElemUid .eq.	841	) area=	1.317816127
          if (jElemUid .eq.	842	) area=	1.958076766
          if (jElemUid .eq.	843	) area=	0.391185199
          if (jElemUid .eq.	844	) area=	2.00250233
          if (jElemUid .eq.	845	) area=	2.88602471
          if (jElemUid .eq.	846	) area=	0.432028626
          if (jElemUid .eq.	847	) area=	1.122821956
          if (jElemUid .eq.	848	) area=	1.647785149
          if (jElemUid .eq.	849	) area=	1.82427611
          if (jElemUid .eq.	850	) area=	2.685168345
          if (jElemUid .eq.	851	) area=	0.923285693
          if (jElemUid .eq.	852	) area=	1.197764762
          if (jElemUid .eq.	853	) area=	1.740070812
          if (jElemUid .eq.	854	) area=	1.558259026
          if (jElemUid .eq.	855	) area=	1.897456833
          if (jElemUid .eq.	856	) area=	1.774241927
          if (jElemUid .eq.	857	) area=	1.072021703
          if (jElemUid .eq.	858	) area=	2.025282828
          if (jElemUid .eq.	859	) area=	1.134984404
          if (jElemUid .eq.	860	) area=	1.917025135
          if (jElemUid .eq.	861	) area=	1.610891541
          if (jElemUid .eq.	862	) area=	0.640956421
          if (jElemUid .eq.	863	) area=	1.123885575
          if (jElemUid .eq.	864	) area=	1.998419038
          if (jElemUid .eq.	865	) area=	2.529755962
          if (jElemUid .eq.	866	) area=	2.133128025
          if (jElemUid .eq.	867	) area=	2.070676005
          if (jElemUid .eq.	868	) area=	1.318719387
          if (jElemUid .eq.	869	) area=	0.371736417
          if (jElemUid .eq.	870	) area=	2.507858923
          if (jElemUid .eq.	871	) area=	2.06313631
          if (jElemUid .eq.	872	) area=	1.090887976
          if (jElemUid .eq.	873	) area=	1.957118821
          if (jElemUid .eq.	874	) area=	1.808833418
          if (jElemUid .eq.	875	) area=	2.773342598
          if (jElemUid .eq.	876	) area=	0.504693909
          if (jElemUid .eq.	877	) area=	0.721174518
          if (jElemUid .eq.	878	) area=	0.134684789
          if (jElemUid .eq.	879	) area=	2.224560729
          if (jElemUid .eq.	880	) area=	0.99485026
          if (jElemUid .eq.	881	) area=	0.036469731
          if (jElemUid .eq.	882	) area=	2.182794395
          if (jElemUid .eq.	883	) area=	2.170782145
          if (jElemUid .eq.	884	) area=	0.77390769
          if (jElemUid .eq.	885	) area=	0.951039298
          if (jElemUid .eq.	886	) area=	0.833332401
          if (jElemUid .eq.	887	) area=	1.714945418
          if (jElemUid .eq.	888	) area=	2.584324845
          if (jElemUid .eq.	889	) area=	0.71171584
          if (jElemUid .eq.	890	) area=	2.161979153
          if (jElemUid .eq.	891	) area=	1.238071126
          if (jElemUid .eq.	892	) area=	1.000633848
          if (jElemUid .eq.	893	) area=	1.981847058
          if (jElemUid .eq.	894	) area=	1.834786109
          if (jElemUid .eq.	895	) area=	2.552623949
          if (jElemUid .eq.	896	) area=	1.961320834
          if (jElemUid .eq.	897	) area=	2.517748978
          if (jElemUid .eq.	898	) area=	2.112455238
          if (jElemUid .eq.	899	) area=	1.508618637
          if (jElemUid .eq.	900	) area=	1.08468914
          if (jElemUid .eq.	901	) area=	2.010112642
          if (jElemUid .eq.	902	) area=	2.29494847
          if (jElemUid .eq.	903	) area=	2.008491684
          if (jElemUid .eq.	904	) area=	1.740965335
          if (jElemUid .eq.	905	) area=	1.749776728
          if (jElemUid .eq.	906	) area=	2.13610454
          if (jElemUid .eq.	907	) area=	1.231761661
          if (jElemUid .eq.	908	) area=	1.611449886
          if (jElemUid .eq.	909	) area=	0.768477594
          if (jElemUid .eq.	910	) area=	0.853236767
          if (jElemUid .eq.	911	) area=	1.52444032
          if (jElemUid .eq.	912	) area=	1.111195644
          if (jElemUid .eq.	913	) area=	0.381566156
          if (jElemUid .eq.	914	) area=	2.361042028
          if (jElemUid .eq.	915	) area=	2.14914463
          if (jElemUid .eq.	916	) area=	1.631005398
          if (jElemUid .eq.	917	) area=	0.078136303
          if (jElemUid .eq.	918	) area=	1.962583499
          if (jElemUid .eq.	919	) area=	1.764261154
          if (jElemUid .eq.	920	) area=	2.039040352
          if (jElemUid .eq.	921	) area=	1.662063606
          if (jElemUid .eq.	922	) area=	2.202533247
          if (jElemUid .eq.	923	) area=	2.467204137
          if (jElemUid .eq.	924	) area=	1.806104535
          if (jElemUid .eq.	925	) area=	0.540352513
          if (jElemUid .eq.	926	) area=	2.159965655
          if (jElemUid .eq.	927	) area=	1.43700613
          if (jElemUid .eq.	928	) area=	1.827572336
          if (jElemUid .eq.	929	) area=	2.636506063
          if (jElemUid .eq.	930	) area=	1.082626595
          if (jElemUid .eq.	931	) area=	1.153474315
          if (jElemUid .eq.	932	) area=	1.337491444
          if (jElemUid .eq.	933	) area=	1.709977184
          if (jElemUid .eq.	934	) area=	1.947221254
          if (jElemUid .eq.	935	) area=	1.907700668
          if (jElemUid .eq.	936	) area=	0.958067668
          if (jElemUid .eq.	937	) area=	1.519798096
          if (jElemUid .eq.	938	) area=	2.747250087
          if (jElemUid .eq.	939	) area=	2.185622766
          if (jElemUid .eq.	940	) area=	1.512355634
          if (jElemUid .eq.	941	) area=	2.251389129
          if (jElemUid .eq.	942	) area=	0.846381035
          if (jElemUid .eq.	943	) area=	0.912091258
          if (jElemUid .eq.	944	) area=	1.781279336
          if (jElemUid .eq.	945	) area=	2.542350353
          if (jElemUid .eq.	946	) area=	0.038930958
          if (jElemUid .eq.	947	) area=	2.592005621
          if (jElemUid .eq.	948	) area=	1.356941573
          if (jElemUid .eq.	949	) area=	2.20889248
          if (jElemUid .eq.	950	) area=	0.785261743
          if (jElemUid .eq.	951	) area=	1.803884114
          if (jElemUid .eq.	952	) area=	1.985198289
          if (jElemUid .eq.	953	) area=	1.826026957
          if (jElemUid .eq.	954	) area=	2.444629902
          if (jElemUid .eq.	955	) area=	1.817917381
          if (jElemUid .eq.	956	) area=	1.513921346
          if (jElemUid .eq.	957	) area=	1.401581461
          if (jElemUid .eq.	958	) area=	2.146901137
          if (jElemUid .eq.	959	) area=	1.732791903
          if (jElemUid .eq.	960	) area=	1.088790372
          if (jElemUid .eq.	961	) area=	0.843316628
          if (jElemUid .eq.	962	) area=	1.633450587
          if (jElemUid .eq.	963	) area=	1.659222315
          if (jElemUid .eq.	964	) area=	1.455550856
          if (jElemUid .eq.	965	) area=	2.277478744
          if (jElemUid .eq.	966	) area=	1.488674926
          if (jElemUid .eq.	967	) area=	0.897006452
          if (jElemUid .eq.	968	) area=	1.671028952
          if (jElemUid .eq.	969	) area=	1.952757861
          if (jElemUid .eq.	970	) area=	2.047175566
          if (jElemUid .eq.	971	) area=	3.030698673
          if (jElemUid .eq.	972	) area=	1.719069802
          if (jElemUid .eq.	973	) area=	0.042347291
          if (jElemUid .eq.	974	) area=	0.584512984
          if (jElemUid .eq.	975	) area=	2.636436179
          if (jElemUid .eq.	976	) area=	1.622386945
          if (jElemUid .eq.	977	) area=	0.802798144
          if (jElemUid .eq.	978	) area=	1.179923582
          if (jElemUid .eq.	979	) area=	1.80796659
          if (jElemUid .eq.	980	) area=	2.588147735
          if (jElemUid .eq.	981	) area=	1.483431576
          if (jElemUid .eq.	982	) area=	1.881323582
          if (jElemUid .eq.	983	) area=	0.571943153
          if (jElemUid .eq.	984	) area=	1.849287538
          if (jElemUid .eq.	985	) area=	2.385687532
          if (jElemUid .eq.	986	) area=	1.898672656
          if (jElemUid .eq.	987	) area=	2.013813886
          if (jElemUid .eq.	988	) area=	1.465207423
          if (jElemUid .eq.	989	) area=	1.587087945
          if (jElemUid .eq.	990	) area=	1.547755584
          if (jElemUid .eq.	991	) area=	1.235158577
          if (jElemUid .eq.	992	) area=	1.069894712
          if (jElemUid .eq.	993	) area=	1.919019816
          if (jElemUid .eq.	994	) area=	0.479222266
          if (jElemUid .eq.	995	) area=	2.394323239
          if (jElemUid .eq.	996	) area=	1.320002552
          if (jElemUid .eq.	997	) area=	2.175527297
          if (jElemUid .eq.	998	) area=	0.313823175
          if (jElemUid .eq.	999	) area=	2.569909516
          if (jElemUid .eq.	1000	) area=	1.54267399
          if (jElemUid .eq.	1001	) area=	1.840684384
          if (jElemUid .eq.	1002	) area=	2.592330673
          if (jElemUid .eq.	1003	) area=	1.123437846
          if (jElemUid .eq.	1004	) area=	1.657174701
          if (jElemUid .eq.	1005	) area=	1.681710165
          if (jElemUid .eq.	1006	) area=	2.007497894
          if (jElemUid .eq.	1007	) area=	0.296225067
          if (jElemUid .eq.	1008	) area=	3.030961898
          if (jElemUid .eq.	1009	) area=	2.565140003
          if (jElemUid .eq.	1010	) area=	1.958141447
          if (jElemUid .eq.	1011	) area=	2.117838119
          if (jElemUid .eq.	1012	) area=	0.959741771
          if (jElemUid .eq.	1013	) area=	2.3244322
          if (jElemUid .eq.	1014	) area=	0.329707247
          if (jElemUid .eq.	1015	) area=	1.825397824
          if (jElemUid .eq.	1016	) area=	2.400433112
          if (jElemUid .eq.	1017	) area=	2.774443661
          if (jElemUid .eq.	1018	) area=	1.024352949
          if (jElemUid .eq.	1019	) area=	2.427693087
          if (jElemUid .eq.	1020	) area=	0.742285273
          if (jElemUid .eq.	1021	) area=	1.901334519
          if (jElemUid .eq.	1022	) area=	2.040219203
          if (jElemUid .eq.	1023	) area=	1.954258419
          if (jElemUid .eq.	1024	) area=	1.782707454
          if (jElemUid .eq.	1025	) area=	0.955814253
          if (jElemUid .eq.	1026	) area=	0.540823238
          if (jElemUid .eq.	1027	) area=	1.473595666
          if (jElemUid .eq.	1028	) area=	1.749864625
          if (jElemUid .eq.	1029	) area=	1.082717926
          if (jElemUid .eq.	1030	) area=	0.450574245
          if (jElemUid .eq.	1031	) area=	1.95135196
          if (jElemUid .eq.	1032	) area=	1.231878083
          if (jElemUid .eq.	1033	) area=	2.51775938
          if (jElemUid .eq.	1034	) area=	1.407598245
          if (jElemUid .eq.	1035	) area=	2.700260454
          if (jElemUid .eq.	1036	) area=	1.47104133
          if (jElemUid .eq.	1037	) area=	1.691401413
          if (jElemUid .eq.	1038	) area=	1.187248781
          if (jElemUid .eq.	1039	) area=	0.73307307
          if (jElemUid .eq.	1040	) area=	2.284852679
          if (jElemUid .eq.	1041	) area=	0.788408493
          if (jElemUid .eq.	1042	) area=	1.355361575
          if (jElemUid .eq.	1043	) area=	1.095596249
          if (jElemUid .eq.	1044	) area=	2.483286836
          if (jElemUid .eq.	1045	) area=	1.971415071
          if (jElemUid .eq.	1046	) area=	1.401784761
          if (jElemUid .eq.	1047	) area=	2.261555289
          if (jElemUid .eq.	1048	) area=	1.57464602
          if (jElemUid .eq.	1049	) area=	1.704496979
          if (jElemUid .eq.	1050	) area=	2.7360902
          if (jElemUid .eq.	1051	) area=	2.259318583
          if (jElemUid .eq.	1052	) area=	0.288848986
          if (jElemUid .eq.	1053	) area=	2.70001999
          if (jElemUid .eq.	1054	) area=	0.839542558
          if (jElemUid .eq.	1055	) area=	0.179114775
          if (jElemUid .eq.	1056	) area=	2.071885997
          if (jElemUid .eq.	1057	) area=	2.615967002
          if (jElemUid .eq.	1058	) area=	1.545762731
          if (jElemUid .eq.	1059	) area=	1.537029713
          if (jElemUid .eq.	1060	) area=	1.958537421
          if (jElemUid .eq.	1061	) area=	1.676391528
          if (jElemUid .eq.	1062	) area=	2.033758123
          if (jElemUid .eq.	1063	) area=	1.446279332
          if (jElemUid .eq.	1064	) area=	1.425068213
          if (jElemUid .eq.	1065	) area=	1.747253497
          if (jElemUid .eq.	1066	) area=	1.334891701
          if (jElemUid .eq.	1067	) area=	1.598974243
          if (jElemUid .eq.	1068	) area=	1.865806993
          if (jElemUid .eq.	1069	) area=	1.287651471
          if (jElemUid .eq.	1070	) area=	0.728651808
          if (jElemUid .eq.	1071	) area=	2.340228597
          if (jElemUid .eq.	1072	) area=	1.703960039
          if (jElemUid .eq.	1073	) area=	1.631214873
          if (jElemUid .eq.	1074	) area=	1.791455177
          if (jElemUid .eq.	1075	) area=	2.212071074
          if (jElemUid .eq.	1076	) area=	0.856038723
          if (jElemUid .eq.	1077	) area=	2.544034827
          if (jElemUid .eq.	1078	) area=	0.055879362
          if (jElemUid .eq.	1079	) area=	0.885699833
          if (jElemUid .eq.	1080	) area=	1.924480271
          if (jElemUid .eq.	1081	) area=	0.680426112
          if (jElemUid .eq.	1082	) area=	2.12071733
          if (jElemUid .eq.	1083	) area=	1.219281631
          if (jElemUid .eq.	1084	) area=	2.524855051
          if (jElemUid .eq.	1085	) area=	2.39767305
          if (jElemUid .eq.	1086	) area=	1.82841626
          if (jElemUid .eq.	1087	) area=	1.155378774
          if (jElemUid .eq.	1088	) area=	1.633583977
          if (jElemUid .eq.	1089	) area=	0.071132536
          if (jElemUid .eq.	1090	) area=	2.182658448
          if (jElemUid .eq.	1091	) area=	2.652625921
          if (jElemUid .eq.	1092	) area=	1.479994484
          if (jElemUid .eq.	1093	) area=	2.68011702
          if (jElemUid .eq.	1094	) area=	0.639039612
          if (jElemUid .eq.	1095	) area=	2.07823919
          if (jElemUid .eq.	1096	) area=	1.536819235
          if (jElemUid .eq.	1097	) area=	1.413356335
          if (jElemUid .eq.	1098	) area=	0.673381065
          if (jElemUid .eq.	1099	) area=	2.393667249
          if (jElemUid .eq.	1100	) area=	0.390658836
          if (jElemUid .eq.	1101	) area=	0.762047741
          if (jElemUid .eq.	1102	) area=	1.515527128
          if (jElemUid .eq.	1103	) area=	0.314805808
          if (jElemUid .eq.	1104	) area=	2.391990917
          if (jElemUid .eq.	1105	) area=	2.155018393
          if (jElemUid .eq.	1106	) area=	2.994328956
          if (jElemUid .eq.	1107	) area=	2.624647617
          if (jElemUid .eq.	1108	) area=	1.065549344
          if (jElemUid .eq.	1109	) area=	1.224535585
          if (jElemUid .eq.	1110	) area=	2.289028299
          if (jElemUid .eq.	1111	) area=	1.971681743
          if (jElemUid .eq.	1112	) area=	1.157804164
          if (jElemUid .eq.	1113	) area=	1.088771169
          if (jElemUid .eq.	1114	) area=	1.372330181
          if (jElemUid .eq.	1115	) area=	2.601941861
          if (jElemUid .eq.	1116	) area=	1.49069553
          if (jElemUid .eq.	1117	) area=	2.219185333
          if (jElemUid .eq.	1118	) area=	1.221449662
          if (jElemUid .eq.	1119	) area=	2.008251622
          if (jElemUid .eq.	1120	) area=	1.383326692
          if (jElemUid .eq.	1121	) area=	1.12383686
          if (jElemUid .eq.	1122	) area=	2.194617819
          if (jElemUid .eq.	1123	) area=	1.454190611
          if (jElemUid .eq.	1124	) area=	2.376874438
          if (jElemUid .eq.	1125	) area=	2.171846132
          if (jElemUid .eq.	1126	) area=	0.395462231
          if (jElemUid .eq.	1127	) area=	2.340421141
          if (jElemUid .eq.	1128	) area=	1.81907758
          if (jElemUid .eq.	1129	) area=	1.640914361
          if (jElemUid .eq.	1130	) area=	2.016178274
          if (jElemUid .eq.	1131	) area=	2.718730026
          if (jElemUid .eq.	1132	) area=	1.857777633
          if (jElemUid .eq.	1133	) area=	0.699651355
          if (jElemUid .eq.	1134	) area=	1.357397666
          if (jElemUid .eq.	1135	) area=	3.010680877
          if (jElemUid .eq.	1136	) area=	2.380338268
          if (jElemUid .eq.	1137	) area=	0.629611914
          if (jElemUid .eq.	1138	) area=	1.364235056
          if (jElemUid .eq.	1139	) area=	1.149431584
          if (jElemUid .eq.	1140	) area=	2.059419178
          if (jElemUid .eq.	1141	) area=	2.278529405
          if (jElemUid .eq.	1142	) area=	2.356272315
          if (jElemUid .eq.	1143	) area=	2.255247193
          if (jElemUid .eq.	1144	) area=	0.444776491
          if (jElemUid .eq.	1145	) area=	3.080050944
          if (jElemUid .eq.	1146	) area=	2.381048085
          if (jElemUid .eq.	1147	) area=	2.724867432
          if (jElemUid .eq.	1148	) area=	2.237073811
          if (jElemUid .eq.	1149	) area=	2.154306383
          if (jElemUid .eq.	1150	) area=	2.388860368
          if (jElemUid .eq.	1151	) area=	1.12072277
          if (jElemUid .eq.	1152	) area=	0.368345377
          if (jElemUid .eq.	1153	) area=	1.124580764
          if (jElemUid .eq.	1154	) area=	1.622184091
          if (jElemUid .eq.	1155	) area=	0.465453299
          if (jElemUid .eq.	1156	) area=	1.644376547
          if (jElemUid .eq.	1157	) area=	1.67419504
          if (jElemUid .eq.	1158	) area=	0.811153692
          if (jElemUid .eq.	1159	) area=	2.118004205
          if (jElemUid .eq.	1160	) area=	0.191270038
          if (jElemUid .eq.	1161	) area=	1.348950614
          if (jElemUid .eq.	1162	) area=	2.056077234
          if (jElemUid .eq.	1163	) area=	1.702957629
          if (jElemUid .eq.	1164	) area=	2.449272161
          if (jElemUid .eq.	1165	) area=	1.985200464
          if (jElemUid .eq.	1166	) area=	2.061840695
          if (jElemUid .eq.	1167	) area=	1.781297629
          if (jElemUid .eq.	1168	) area=	0.030653677
          if (jElemUid .eq.	1169	) area=	1.370087435
          if (jElemUid .eq.	1170	) area=	1.468187086
          if (jElemUid .eq.	1171	) area=	2.643437956
          if (jElemUid .eq.	1172	) area=	1.939051678
          if (jElemUid .eq.	1173	) area=	1.711221765
          if (jElemUid .eq.	1174	) area=	1.57179496
          if (jElemUid .eq.	1175	) area=	2.040927249
          if (jElemUid .eq.	1176	) area=	1.932493045
          if (jElemUid .eq.	1177	) area=	2.402059584
          if (jElemUid .eq.	1178	) area=	1.069871655
          if (jElemUid .eq.	1179	) area=	2.182599311
          if (jElemUid .eq.	1180	) area=	0.571107371
          if (jElemUid .eq.	1181	) area=	0.257514718
          if (jElemUid .eq.	1182	) area=	2.756421777
          if (jElemUid .eq.	1183	) area=	1.013588086
          if (jElemUid .eq.	1184	) area=	0.025624052
          if (jElemUid .eq.	1185	) area=	1.539035185
          if (jElemUid .eq.	1186	) area=	2.017243804
          if (jElemUid .eq.	1187	) area=	1.703580906
          if (jElemUid .eq.	1188	) area=	1.63321978
          if (jElemUid .eq.	1189	) area=	1.116676337
          if (jElemUid .eq.	1190	) area=	0.481546263
          if (jElemUid .eq.	1191	) area=	2.577190079
          if (jElemUid .eq.	1192	) area=	2.601964021
          if (jElemUid .eq.	1193	) area=	2.562002343
          if (jElemUid .eq.	1194	) area=	1.789440332
          if (jElemUid .eq.	1195	) area=	1.93610398
          if (jElemUid .eq.	1196	) area=	1.947227003
          if (jElemUid .eq.	1197	) area=	1.007391239
          if (jElemUid .eq.	1198	) area=	1.961715091
          if (jElemUid .eq.	1199	) area=	0.006608962
          if (jElemUid .eq.	1200	) area=	2.859634787
          if (jElemUid .eq.	1201	) area=	1.678042859
          if (jElemUid .eq.	1202	) area=	1.90873653
          if (jElemUid .eq.	1203	) area=	0.213120188
          if (jElemUid .eq.	1204	) area=	1.018238772
          if (jElemUid .eq.	1205	) area=	2.547046019
          if (jElemUid .eq.	1206	) area=	1.965240209
          if (jElemUid .eq.	1207	) area=	1.287379227
          if (jElemUid .eq.	1208	) area=	1.881619468
          if (jElemUid .eq.	1209	) area=	0.63675713
          if (jElemUid .eq.	1210	) area=	2.330010147
          if (jElemUid .eq.	1211	) area=	2.540528633
          if (jElemUid .eq.	1212	) area=	2.491912154
          if (jElemUid .eq.	1213	) area=	0.577532551
          if (jElemUid .eq.	1214	) area=	0.581623355
          if (jElemUid .eq.	1215	) area=	2.116200275
          if (jElemUid .eq.	1216	) area=	0.335494263
          if (jElemUid .eq.	1217	) area=	0.494306595
          if (jElemUid .eq.	1218	) area=	1.867481035
          if (jElemUid .eq.	1219	) area=	1.964468597
          if (jElemUid .eq.	1220	) area=	0.58238799
          if (jElemUid .eq.	1221	) area=	0.154275838
          if (jElemUid .eq.	1222	) area=	1.011381249
          if (jElemUid .eq.	1223	) area=	1.894816146
          if (jElemUid .eq.	1224	) area=	2.305106944
          if (jElemUid .eq.	1225	) area=	0.607259243
          if (jElemUid .eq.	1226	) area=	0.472449231
          if (jElemUid .eq.	1227	) area=	0.890771957
          if (jElemUid .eq.	1228	) area=	1.967193671
          if (jElemUid .eq.	1229	) area=	2.191571889
          if (jElemUid .eq.	1230	) area=	1.671322523
          if (jElemUid .eq.	1231	) area=	1.763984682
          if (jElemUid .eq.	1232	) area=	0.969167968
          if (jElemUid .eq.	1233	) area=	1.821198963
          if (jElemUid .eq.	1234	) area=	0.374396524
          if (jElemUid .eq.	1235	) area=	1.996085212
          if (jElemUid .eq.	1236	) area=	0.452096886
          if (jElemUid .eq.	1237	) area=	1.882817656
          if (jElemUid .eq.	1238	) area=	1.061220225
          if (jElemUid .eq.	1239	) area=	1.782023794
          if (jElemUid .eq.	1240	) area=	1.698674883
          if (jElemUid .eq.	1241	) area=	2.175002557
          if (jElemUid .eq.	1242	) area=	2.906230363
          if (jElemUid .eq.	1243	) area=	0.437530677
          if (jElemUid .eq.	1244	) area=	1.512255
          if (jElemUid .eq.	1245	) area=	2.323125207
          if (jElemUid .eq.	1246	) area=	0.477903536
          if (jElemUid .eq.	1247	) area=	1.680291903
          if (jElemUid .eq.	1248	) area=	1.719185003
          if (jElemUid .eq.	1249	) area=	0.985383485
          if (jElemUid .eq.	1250	) area=	0.878226911
          if (jElemUid .eq.	1251	) area=	0.826505197
          if (jElemUid .eq.	1252	) area=	2.089272029
          if (jElemUid .eq.	1253	) area=	2.095913601
          if (jElemUid .eq.	1254	) area=	1.49255348
          if (jElemUid .eq.	1255	) area=	2.011196336
          if (jElemUid .eq.	1256	) area=	0.785404193
          if (jElemUid .eq.	1257	) area=	2.464432022
          if (jElemUid .eq.	1258	) area=	2.016663807
          if (jElemUid .eq.	1259	) area=	2.274244564
          if (jElemUid .eq.	1260	) area=	1.796982684
          if (jElemUid .eq.	1261	) area=	1.996546275
          if (jElemUid .eq.	1262	) area=	0.521535326
          if (jElemUid .eq.	1263	) area=	0.678710416
          if (jElemUid .eq.	1264	) area=	1.756349722
          if (jElemUid .eq.	1265	) area=	1.748412474
          if (jElemUid .eq.	1266	) area=	1.943046459
          if (jElemUid .eq.	1267	) area=	1.052629469
          if (jElemUid .eq.	1268	) area=	2.269259767
          if (jElemUid .eq.	1269	) area=	0.144065184
          if (jElemUid .eq.	1270	) area=	0.996344475
          if (jElemUid .eq.	1271	) area=	1.386273188
          if (jElemUid .eq.	1272	) area=	0.544666098
          if (jElemUid .eq.	1273	) area=	1.065695144
          if (jElemUid .eq.	1274	) area=	1.950025786
          if (jElemUid .eq.	1275	) area=	1.076709506
          if (jElemUid .eq.	1276	) area=	1.848383734
          if (jElemUid .eq.	1277	) area=	1.60520014
          if (jElemUid .eq.	1278	) area=	2.012441181
          if (jElemUid .eq.	1279	) area=	1.413924171
          if (jElemUid .eq.	1280	) area=	1.451214756
          if (jElemUid .eq.	1281	) area=	1.063687655
          if (jElemUid .eq.	1282	) area=	2.298808755
          if (jElemUid .eq.	1283	) area=	1.963078047
          if (jElemUid .eq.	1284	) area=	1.201145235
          if (jElemUid .eq.	1285	) area=	2.796006502
          if (jElemUid .eq.	1286	) area=	2.029109083
          if (jElemUid .eq.	1287	) area=	1.456447165
          if (jElemUid .eq.	1288	) area=	1.886079957
          if (jElemUid .eq.	1289	) area=	2.79800146
          if (jElemUid .eq.	1290	) area=	2.56226258
          if (jElemUid .eq.	1291	) area=	2.084752609
          if (jElemUid .eq.	1292	) area=	1.744529202
          if (jElemUid .eq.	1293	) area=	2.178669683
          if (jElemUid .eq.	1294	) area=	0.145810711
          if (jElemUid .eq.	1295	) area=	1.914993611
          if (jElemUid .eq.	1296	) area=	1.299238845
          if (jElemUid .eq.	1297	) area=	0.852238193
          if (jElemUid .eq.	1298	) area=	0.650632246
          if (jElemUid .eq.	1299	) area=	2.532500992
          if (jElemUid .eq.	1300	) area=	2.44062899
          if (jElemUid .eq.	1301	) area=	1.110244965
          if (jElemUid .eq.	1302	) area=	2.898179084
          if (jElemUid .eq.	1303	) area=	0.006785911
          if (jElemUid .eq.	1304	) area=	1.080185912
          if (jElemUid .eq.	1305	) area=	2.100114909
          if (jElemUid .eq.	1306	) area=	2.273196945
          if (jElemUid .eq.	1307	) area=	2.209709922
          if (jElemUid .eq.	1308	) area=	1.369120844
          if (jElemUid .eq.	1309	) area=	2.263024446
          if (jElemUid .eq.	1310	) area=	1.968990334
          if (jElemUid .eq.	1311	) area=	2.546441536
          if (jElemUid .eq.	1312	) area=	2.488168551
          if (jElemUid .eq.	1313	) area=	0.569710252
          if (jElemUid .eq.	1314	) area=	2.435361826
          if (jElemUid .eq.	1315	) area=	0.627995796
          if (jElemUid .eq.	1316	) area=	1.977597801
          if (jElemUid .eq.	1317	) area=	0.471000485
          if (jElemUid .eq.	1318	) area=	1.968635047
          if (jElemUid .eq.	1319	) area=	1.808832262
          if (jElemUid .eq.	1320	) area=	0.856860569
          if (jElemUid .eq.	1321	) area=	2.451831145
          if (jElemUid .eq.	1322	) area=	0.546959279
          if (jElemUid .eq.	1323	) area=	1.139189874
          if (jElemUid .eq.	1324	) area=	2.477300114
          if (jElemUid .eq.	1325	) area=	2.863263038
          if (jElemUid .eq.	1326	) area=	1.529456457
          if (jElemUid .eq.	1327	) area=	1.499881554
          if (jElemUid .eq.	1328	) area=	1.811217371
          if (jElemUid .eq.	1329	) area=	2.581802813
          if (jElemUid .eq.	1330	) area=	0.941311737
          if (jElemUid .eq.	1331	) area=	1.734330183
          if (jElemUid .eq.	1332	) area=	0.171451869
          if (jElemUid .eq.	1333	) area=	2.278557883
          if (jElemUid .eq.	1334	) area=	3.065729858
          if (jElemUid .eq.	1335	) area=	1.848953618
          if (jElemUid .eq.	1336	) area=	1.064836407
          if (jElemUid .eq.	1337	) area=	0.529351821
          if (jElemUid .eq.	1338	) area=	2.154509434
          if (jElemUid .eq.	1339	) area=	1.073434316
          if (jElemUid .eq.	1340	) area=	1.040904448
          if (jElemUid .eq.	1341	) area=	0.152293844
          if (jElemUid .eq.	1342	) area=	1.600085331
          if (jElemUid .eq.	1343	) area=	0.781481356
          if (jElemUid .eq.	1344	) area=	0.823328377
          if (jElemUid .eq.	1345	) area=	2.871493418
          if (jElemUid .eq.	1346	) area=	0.336911777
          if (jElemUid .eq.	1347	) area=	1.797998921
          if (jElemUid .eq.	1348	) area=	1.813570742
          if (jElemUid .eq.	1349	) area=	2.791069047
          if (jElemUid .eq.	1350	) area=	1.413603685
          if (jElemUid .eq.	1351	) area=	1.918590677
          if (jElemUid .eq.	1352	) area=	0.410100877
          if (jElemUid .eq.	1353	) area=	2.449485705
          if (jElemUid .eq.	1354	) area=	2.178751993
          if (jElemUid .eq.	1355	) area=	1.757644428
          if (jElemUid .eq.	1356	) area=	1.286685595
          if (jElemUid .eq.	1357	) area=	2.022358641
          if (jElemUid .eq.	1358	) area=	2.264452731
          if (jElemUid .eq.	1359	) area=	1.421531249
          if (jElemUid .eq.	1360	) area=	2.489374414
          if (jElemUid .eq.	1361	) area=	1.89671887
          if (jElemUid .eq.	1362	) area=	1.938025521
          if (jElemUid .eq.	1363	) area=	1.411504316
          if (jElemUid .eq.	1364	) area=	2.121476894
          if (jElemUid .eq.	1365	) area=	1.603333527
          if (jElemUid .eq.	1366	) area=	2.530823856
          if (jElemUid .eq.	1367	) area=	2.302592183
          if (jElemUid .eq.	1368	) area=	1.186705139
          if (jElemUid .eq.	1369	) area=	1.61685478
          if (jElemUid .eq.	1370	) area=	0.680050047
          if (jElemUid .eq.	1371	) area=	0.316158928
          if (jElemUid .eq.	1372	) area=	0.932380056
          if (jElemUid .eq.	1373	) area=	0.84117009
          if (jElemUid .eq.	1374	) area=	1.106641112
          if (jElemUid .eq.	1375	) area=	1.884412176
          if (jElemUid .eq.	1376	) area=	0.808727008
          if (jElemUid .eq.	1377	) area=	0.753754584
          if (jElemUid .eq.	1378	) area=	2.721345998
          if (jElemUid .eq.	1379	) area=	0.507459421
          if (jElemUid .eq.	1380	) area=	2.335065308
          if (jElemUid .eq.	1381	) area=	1.527444089
          if (jElemUid .eq.	1382	) area=	2.81225529
          if (jElemUid .eq.	1383	) area=	0.564578144
          if (jElemUid .eq.	1384	) area=	2.874970837
          if (jElemUid .eq.	1385	) area=	0.768700379
          if (jElemUid .eq.	1386	) area=	2.417693543
          if (jElemUid .eq.	1387	) area=	1.480326579
          if (jElemUid .eq.	1388	) area=	2.357678872
          if (jElemUid .eq.	1389	) area=	0.088476743
          if (jElemUid .eq.	1390	) area=	1.094060954
          if (jElemUid .eq.	1391	) area=	1.517592334
          if (jElemUid .eq.	1392	) area=	1.82775681
          if (jElemUid .eq.	1393	) area=	2.092008482
          if (jElemUid .eq.	1394	) area=	0.013590687
          if (jElemUid .eq.	1395	) area=	0.36719448
          if (jElemUid .eq.	1396	) area=	2.730932722
          if (jElemUid .eq.	1397	) area=	2.809984327
          if (jElemUid .eq.	1398	) area=	1.42092692
          if (jElemUid .eq.	1399	) area=	1.997918596
          if (jElemUid .eq.	1400	) area=	2.272682271
          if (jElemUid .eq.	1401	) area=	2.018356723
          if (jElemUid .eq.	1402	) area=	2.236387854
          if (jElemUid .eq.	1403	) area=	1.353599232
          if (jElemUid .eq.	1404	) area=	2.121146248
          if (jElemUid .eq.	1405	) area=	0.27509694
          if (jElemUid .eq.	1406	) area=	2.51904458
          if (jElemUid .eq.	1407	) area=	1.017927822
          if (jElemUid .eq.	1408	) area=	1.006924514
          if (jElemUid .eq.	1409	) area=	2.365039442
          if (jElemUid .eq.	1410	) area=	1.750626862
          if (jElemUid .eq.	1411	) area=	1.420331746
          if (jElemUid .eq.	1412	) area=	1.654435751
          if (jElemUid .eq.	1413	) area=	1.671323407
          if (jElemUid .eq.	1414	) area=	0.701990313
          if (jElemUid .eq.	1415	) area=	2.592738985
          if (jElemUid .eq.	1416	) area=	1.284563264
          if (jElemUid .eq.	1417	) area=	2.931989021
          if (jElemUid .eq.	1418	) area=	2.214098945
          if (jElemUid .eq.	1419	) area=	0.233354388
          if (jElemUid .eq.	1420	) area=	0.855307345
          if (jElemUid .eq.	1421	) area=	1.70214048
          if (jElemUid .eq.	1422	) area=	0.832838496
          if (jElemUid .eq.	1423	) area=	0.907642474
          if (jElemUid .eq.	1424	) area=	2.32545161
          if (jElemUid .eq.	1425	) area=	2.15699774
          if (jElemUid .eq.	1426	) area=	0.972677612
          if (jElemUid .eq.	1427	) area=	1.171398156
          if (jElemUid .eq.	1428	) area=	1.360369302
          if (jElemUid .eq.	1429	) area=	2.267726481
          if (jElemUid .eq.	1430	) area=	2.450953524
          if (jElemUid .eq.	1431	) area=	1.413551251
          if (jElemUid .eq.	1432	) area=	2.031456803
          if (jElemUid .eq.	1433	) area=	1.04644556
          if (jElemUid .eq.	1434	) area=	1.666363463
          if (jElemUid .eq.	1435	) area=	1.739927634
          if (jElemUid .eq.	1436	) area=	0.582678144
          if (jElemUid .eq.	1437	) area=	2.052814051
          if (jElemUid .eq.	1438	) area=	1.176331358
          if (jElemUid .eq.	1439	) area=	1.211550501
          if (jElemUid .eq.	1440	) area=	2.507214435
          if (jElemUid .eq.	1441	) area=	0.506303138
          if (jElemUid .eq.	1442	) area=	0.98553636
          if (jElemUid .eq.	1443	) area=	0.718151953
          if (jElemUid .eq.	1444	) area=	0.017549313
          if (jElemUid .eq.	1445	) area=	2.049695583
          if (jElemUid .eq.	1446	) area=	1.468997625
          if (jElemUid .eq.	1447	) area=	2.03272849
          if (jElemUid .eq.	1448	) area=	2.195725633
          if (jElemUid .eq.	1449	) area=	1.061551956
          if (jElemUid .eq.	1450	) area=	1.59491813
          if (jElemUid .eq.	1451	) area=	1.620733955
          if (jElemUid .eq.	1452	) area=	1.190381219
          if (jElemUid .eq.	1453	) area=	1.314111857
          if (jElemUid .eq.	1454	) area=	1.863142053
          if (jElemUid .eq.	1455	) area=	1.638312048
          if (jElemUid .eq.	1456	) area=	0.416819199
          if (jElemUid .eq.	1457	) area=	0.070177963
          if (jElemUid .eq.	1458	) area=	1.534572416
          if (jElemUid .eq.	1459	) area=	1.443676215
          if (jElemUid .eq.	1460	) area=	0.844320179
          if (jElemUid .eq.	1461	) area=	1.5257872
          if (jElemUid .eq.	1462	) area=	0.076384631
          if (jElemUid .eq.	1463	) area=	0.625844732
          if (jElemUid .eq.	1464	) area=	0.582439024
          if (jElemUid .eq.	1465	) area=	1.97393246
          if (jElemUid .eq.	1466	) area=	0.02035388
          if (jElemUid .eq.	1467	) area=	2.186797851
          if (jElemUid .eq.	1468	) area=	1.212218273
          if (jElemUid .eq.	1469	) area=	0.608103693
          if (jElemUid .eq.	1470	) area=	1.930733998
          if (jElemUid .eq.	1471	) area=	1.936260067
          if (jElemUid .eq.	1472	) area=	1.670052011
          if (jElemUid .eq.	1473	) area=	1.631833335
          if (jElemUid .eq.	1474	) area=	0.243229155
          if (jElemUid .eq.	1475	) area=	1.890424122
          if (jElemUid .eq.	1476	) area=	0.968607389
          if (jElemUid .eq.	1477	) area=	2.039899614
          if (jElemUid .eq.	1478	) area=	0.672535943
          if (jElemUid .eq.	1479	) area=	1.998617392
          if (jElemUid .eq.	1480	) area=	0.533355653
          if (jElemUid .eq.	1481	) area=	1.462222069
          if (jElemUid .eq.	1482	) area=	0.824997194
          if (jElemUid .eq.	1483	) area=	2.173629566
          if (jElemUid .eq.	1484	) area=	1.506597342
          if (jElemUid .eq.	1485	) area=	0.166445485
          if (jElemUid .eq.	1486	) area=	1.302893969
          if (jElemUid .eq.	1487	) area=	1.378638306
          if (jElemUid .eq.	1488	) area=	1.010190644
          if (jElemUid .eq.	1489	) area=	0.315720698
  
          allArea(jElemUid) = area
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         variables transfer
          stateNew(i,1)  = allsigLocX(jElemUid)
          stateNew(i,2)  = allsigLocY(jElemUid)
          stateNew(i,3)  = allstressLocX(jElemUid)
          stateNew(i,4)  = allstressLocY(jElemUid)
          stateNew(i,5)  = area
          stateNew(i,6)  = allLE11(jElemUid)
          stateNew(i,7)  = allLE22(jElemUid)
          stateNew(i,9)  = allFt(jElemUid)
*          
          stateNew(i,14) = allH(jElemUid) 
          stateNew(i,15) = allD(jElemUid)
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
  100 continue
*     ================================================================================================     *
      return
      end

      
*     *************************************************************************************************
*
*     VUEL for the displacement field  
*     Activated DOFs: 1 and 2
*
*     *************************************************************************************************
      subroutine vuel(
     *     nblock,
c          to be defined
     *     rhs,amass,dtimeStable,
     *     svars,nsvars,
     *     energy,
c          
     *     nnode,ndofel,
     *     props,nprops,
     *     jprops,njprops,
     *     coords,ncrd,
     *     u,du,v,a,
     *     jtype,jelem,
     *     time,period,dtimeCur,dtimePrev,kstep,kinc,lflags,
     *     dMassScaleFactor,
     *     predef,npredef,
     *     jdltyp,adlmag)
C     
      use vars_module
      include 'vaba_param.inc'

      parameter ( zero = 0.d0, half = 0.5d0, one = 1.d0, two=2.d0 )

c     operation code
      parameter ( jMassCalc            = 1,
     *            jIntForceAndDtStable = 2,
     *            jExternForce         = 3)

c     flags
      parameter (iProcedure = 1,
     *           iNlgeom    = 2,
     *           iOpCode    = 3,
     *           nFlags     = 3)

c     time
      parameter (iStepTime  = 1,
     *           iTotalTime = 2,
     *           nTime      = 2)

c     procedure flags
      parameter ( jDynExplicit = 17 )

c     energies 
      parameter ( iElPd = 1,
     *            iElCd = 2,
     *            iElIe = 3,
     *            iElTs = 4,
     *            iElDd = 5,
     *            iElBv = 6,
     *            iElDe = 7,
     *            iElHe = 8,
     *            iElKe = 9,
     *            iElTh = 10,
     *            iElDmd = 11,
     *            iElDc = 12,
     *            nElEnergy = 12)

c     predefined variables
      parameter ( iPredValueNew = 1,
     *            iPredValueOld = 2,
     *            nPred         = 2)    

c     indexing in a 3-long vector

      parameter (factorStable = 0.99d0)
      
      dimension rhs(nblock,ndofel), amass(nblock,ndofel,ndofel),
     *     dtimeStable(nblock),
     *     svars(nblock,nsvars), energy(nblock,nElEnergy),
     *     props(nprops), jprops(njprops),
     *     jelem(nblock), time(nTime), lflags(nFlags),
     *     coords(nblock,nnode,ncrd), u(nblock,ndofel),
     *     du(nblock,ndofel), v(nblock,ndofel), a(nblock, ndofel),
     *     dMassScaleFactor(nblock),
     *     predef(nblock, nnode, npredef, nPred), adlmag(nblock)
      
      dimension Rot(ndofel,ndofel), u_loc(ndofel), u_glo(ndofel)
      dimension B(nnode,ndofel),v_loc(ndofel), v_glo(ndofel)
      dimension Rot_tran(ndofel,ndofel)
      dimension force_loc(2),force_loc_elem(4), force_glo(4)
C      
      dimension UD(2),DUD(2)
c      
      dimension S(3),SPos(3),SNeg(3)
C
      integer i,j,k
C
      real*8  lc,lb, Gf, ea, ft 
      real*8  p, a1, a2, a3
      real*8  phase, dalpha, dphase
      real*8  omega,domega
      real*8  hist
      real*8  c0
*     ================================================================================================     *
*     parameter initilization
*
      eta = 1.0d-4
*          
      ea   =  props(1)      ! props(1) -- Young's modulus
      ft   =  props(2)      ! props(3) -- failure strength
      Gf   =  props(3)      ! props(4) -- fracture energy
      lb   =  props(4)      ! props(5) -- length scale  
*      
      c0   =  3.1415926535897932384626433832d0
*
      area0 = 1.0d0
      rho   = 1.0d-4
*
      eDampTra    = 0.1d0        
      amassFact0  = half*area0*rho
*
      if ( lflags(iOpCode).eq.jMassCalc ) then
*     ================================================================================================     *
          do kblock = 1, nblock
*
*             coords(nblock,nnode,ncrd)
              alenX0 = (coords(kblock,2,1) - coords(kblock,1,1))
              alenY0 = (coords(kblock,2,2) - coords(kblock,1,2))
              alen0  = sqrt(alenX0*alenX0 + alenY0*alenY0)
*               
              am0               = amassFact0*alen0
              amass(kblock,1,1) = am0
              amass(kblock,2,2) = am0
              amass(kblock,4,4) = am0
              amass(kblock,5,5) = am0
*              
              amass(kblock,3,3) = eta * alen0 / 2.0d0
              amass(kblock,6,6) = eta * alen0 / 2.0d0
*              
          end do
*     ================================================================================================     *
      else if ( lflags(iOpCode) .eq. jIntForceAndDtStable) then
*     ================================================================================================     *
          do kblock = 1, nblock
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
*         here begins the displacement field
*              
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             Parameter initialization 
              Rot = 0.0d0
              Rot_tran = 0.0d0
              force_loc = 0.0d0
              force_loc_elem = 0.0d0
              force_glo = 0.0d0
              u_glo = 0.0d0
              u_loc = 0.0d0
*              
              a1   =  4.d0/(c0*lb)*ea*Gf/(ft*ft)
              p    =  2.0d0
              a2   =  1.3868d0
              a3   =  0.6567d0  
              
              WF = Gf/2.0/lb      
              WO = zero    
*          
              area = allArea( jelem(kblock) - NumEle)
             
              alenX0 = (coords(kblock,2,1) - coords(kblock,1,1))
              alenY0 = (coords(kblock,2,2) - coords(kblock,1,2))
              alen0  = sqrt(alenX0*alenX0 + alenY0*alenY0)
*
              cos_theta = alenX0/alen0
              sin_theta = alenY0/alen0
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             Rotation matrix 
              Rot(1,1) =   cos_theta 
              Rot(1,2) =   sin_theta
              Rot(2,1) = - sin_theta
              Rot(2,2) =   cos_theta 
              Rot(3,3) =   cos_theta 
              Rot(3,4) =   sin_theta
              Rot(4,3) = - sin_theta
              Rot(4,4) =   cos_theta 
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             global displacement 
*
              u_glo(1) = u(kblock,1)
              u_glo(2) = u(kblock,2)  
              u_glo(3) = u(kblock,4)
              u_glo(4) = u(kblock,5)
*             global velocity  
*
              v_glo(1) = v(kblock,1)
              v_glo(2) = v(kblock,2)
              v_glo(3) = v(kblock,4)
              v_glo(4) = v(kblock,5)
*             local displacement & velocity 
*
              u_loc = matmul(Rot, u_glo)   
              v_loc = matmul(Rot, v_glo) 
*             local strain    
*
              strainLocX = (u_loc(3) - u_loc(1)) / alen0
              strainLocY = (u_loc(4) - u_loc(2)) / alen0
*             local velocity 
*
              aVelX = (v_loc(3) - v_loc(1))
              aVelY = (v_loc(4) - v_loc(2))
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             Local stress 
*
              phase = svars(kblock,15)
*              
              call energeticFunc(omega,domega,phase,a1,a2,a3,p)
*              
              sigLocX = ea * strainLocX /(1.0d0 - 0.2d0)
              sigLocY = ea * strainLocY *(1.0d0 - 0.6d0)/(1.0d0 + 0.2d0)
              
              S(1) = sigLocX
              S(2) = 0.0
              S(3) = sigLocY
              
              PsPos = S(1) / 2.0d0 + sqrt( (S(1) / 2.0d0 )**2.0d0 + 
     +                    S(3)**2.0d0 )
              PsPos = max(PsPos, 0.0d0)
              
              PsNeg = S(1) / 2.0d0 - sqrt( (S(1) / 2.0d0 )**2.0d0 + 
     +                    S(3)**2.0d0 )
              PsNeg = min(PsNeg, 0.0d0)
*     
              if( abs(S(1)) .ne. 0.0d0) then
                  stheta = atan(2.0*s(3) / s(1)) /2.0
              else
                  stheta = sign( 45.0/180.0*3.1415926 , S(3))
              end if
              
              SPos(1) = PsPos * cos (-stheta)**2.0d0 
              SPos(2) = PsPos * sin (-stheta)**2.0d0 
              SPos(3) = (-PsPos) * sin (-stheta)*cos (-stheta)
*
              SNeg(1) = S(1) - SPos(1)
              SNeg(2) = S(2) - SPos(2)
              SNeg(3) = S(3) - SPos(3)
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
              stressLocX = SPos(1) * omega + SNeg(1)
              stressLocY = SPos(3) * omega + SNeg(3)
*     
*             local force vector 
*
              force_loc(1)      = area * stressLocX + eDampTra * aVelX
              force_loc(2)      = area * stressLocY + eDampTra * aVelY
*              force_loc(2)      = 0.0d0
*               
              force_loc_elem(1) = - force_loc(1) 
              force_loc_elem(2) = - force_loc(2)
              force_loc_elem(3) =   force_loc(1)
              force_loc_elem(4) =   force_loc(2)
*               
              Rot_tran = transpose(Rot)
*     
*             global force vector 
*               
              force_glo = matmul(Rot_tran, force_loc_elem)
*
*     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *  
*             stress power 
*              
              stressPower = 0.5d0 * PsPos * PsPos / ea 
*                            
*              if (strainLocX .gt. 0.0) then
*                  stressPower = 0.5d0 * ea * strainLocX * strainLocX
*              else
*                  stressPower = 0.0d0
*              end if
*
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
*         here begins the phase field
*              
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
              hist =  svars(kblock,14)
              if(hist .lt. zero)  hist = zero 
*
              UD(1)  = U(kblock,3)
              UD(2)  = U(kblock,6)
*
              phase   = (UD(1) + UD(2)) / 2.0d0
              dphase  = (UD(2) - UD(1)) / alen0
*                           
              if (phase .le. 0.0) then
                 phase = 0.0d0
              end if
*             
              if (phase .ge. 1.0d0) then
                 phase = 1.0d0
                 dphase = 0.0d0
              end if
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             threshold of the power
*
              call geometricFunc(dalpha,0.0d0) 
              call energeticFunc(omega,domega,0.0d0,a1,a2,a3,p)  
              hist_0 = - Gf * dalpha / domega / c0 / lb
*            
              hist = max(hist, hist_0)             
*
              call geometricFunc(dalpha,phase)                         ! geometric function
              call energeticFunc(omega,domega,phase,a1,a2,a3,p)        ! energetic function
*
              phi_source  = domega *  hist + Gf/(c0*lb)*dalpha
*
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
*             right hand side vector 
*              
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
              rhs(kblock,1) = force_glo(1)
              rhs(kblock,2) = force_glo(2)
              rhs(kblock,4) = force_glo(3)
              rhs(kblock,5) = force_glo(4)
*
              rhs(kblock,3) =   2.d0*lb*Gf/c0*dphase/alen0 - phi_source
              rhs(kblock,6) = - 2.d0*lb*Gf/c0*dphase/alen0 - phi_source
*             
              rhs(kblock,3) =  - rhs(kblock,3) * alen0
              rhs(kblock,6) =  - rhs(kblock,6) * alen0
             
              phaseOld =  svars(kblock,15)
              if ( (hist .le. hist_0) )then
                  rhs(kblock,3) =  0.0d0
                  rhs(kblock,6) = 0.0d0
              end if
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
*             critical time  
*              
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             Critical time     
              amElem0 = two*amassFact0*alen0
*
              ak     = area0*ea/alen0
*             undamped stable time increment for translations
              dtTrialTransl = sqrt(amElem0/ak)

*             damped stable time increment; since eDampTra=0, the
*             stable time increment does not change because of damping
              critDampTransl = two*sqrt(amElem0*ak)
              csiTra = eDampTra/critDampTransl
              factDamp = sqrt(one+csiTra*csiTra) - csiTra
              dtTrialTransl = dtTrialTransl*factDamp*factorStable
              dtimeStable1 = dtTrialTransl
*
              dtimeStable2 = factorStable*eta/(500.0*(WF-WO))
*             
              dtimeStable(kblock) = min(dtimeStable1/10.0, dtimeStable2)
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*
*             Internal state variables  
*              
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
*             stress power
*
              svars(kblock,11) =  stressPower
              svars(kblock,14) =  max(svars(kblock,11),svars(kblock,14))
              allH(jelem(kblock) - NumEle) = svars(kblock,14)
*             phase field
*              
              svars(kblock,15) =  phase
              allD(jelem(kblock) - NumEle) = phase  
*             phase field
*
              allsigLocX(jelem(kblock) - NumEle)     =  SPos(1)
              allsigLocY(jelem(kblock) - NumEle)     =  SPos(3)
              allstressLocX(jelem(kblock) - NumEle)  =  stressLocX
              allstressLocY(jelem(kblock) - NumEle)  =  stressLocY
              allLE11(jelem(kblock) - NumEle)  =  strainLocX
              allLE22(jelem(kblock) - NumEle)  =  strainLocY
              
              allFt(jelem(kblock) - NumEle)  =  ft
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *          
          end do
*         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     *
      end if
*     ================================================================================================     *
      return
      end      
      



      
*     *************************************************************************************************
*
*     VUSDFLD for converting the state variables in VUMAT into global variables  
*     
*     *************************************************************************************************
*       
      subroutine vusdfld(
* Read only -
     *   nblock, nstatev, nfieldv, nprops, ndir, nshr, 
     *   jElemUid, kIntPt, kLayer, kSecPt, 
     *   stepTime, totalTime, dt, cmname, 
     *   coordMp, direct, T, charLength, props, 
     *   stateOld, 
* Write only -
     *   stateNew, field )
*
      include 'vaba_param.inc'
*
      dimension props(nprops),
     *          jElemUid(nblock), coordMp(nblock, *), 
     *          direct(nblock, 3, 3), T(nblock,3,3), 
     *          charLength(nblock),
     *          stateOld(nblock, nstatev), 
     *          stateNew(nblock, nstatev),
     *          field(nblock, nfieldv)
      character*80 cmname
*
      character*3 cData(maxblk)
      dimension jData(maxblk)
      dimension eqps(maxblk)
*
      parameter ( zero = 0.d0 )
*     ================================================================================================     *
      do k = 1, nblock
            stateOld(k,8)   =  jElemUid(k)
            stateNew(k,8)   =  jElemUid(k)
      end do
*     ================================================================================================     *
      return
      end subroutine
      
       
*     *************************************************************************************************
*
*     energetic degradation function  omega  
*     
*     *************************************************************************************************
*        
      subroutine energeticFunc(omega,domega,phi,a1,a2,a3,p)
*
      include 'vaba_param.inc'
      
      real*8:: omega, domega, phi
      real*8:: fac1, dfac1, ddfac1, fac2, dfac2, ddfac2
      real*8:: p, a1, a2, a3
*     ================================================================================================     *
*      
      fac1    =  (1.d0 - phi)**p
      dfac1   = - p*(1.d0 - phi)**(p - 1.d0); 
      ddfac1  =  p*(p - 1.d0)*(1.d0 - phi)**(p - 2.d0)
*        
      fac2   =  fac1   + a1*phi + a1*a2*phi**2.d0 + a1*a2*a3*phi**3.d0
      dfac2  =  dfac1  + a1 + 2.d0*a1*a2*phi + 3.d0*a1*a2*a3*phi**2.d0
      ddfac2  =  ddfac1 + 2.d0*a1*a2 + 6.d0*a1*a2*a3*phi
*        
      omega   =  fac1/fac2        
      domega  =  (dfac1*fac2  - fac1*dfac2)/(fac2**2.d0)   
*     ================================================================================================     *
      return
      end subroutine energeticFunc
      
***********************************************************************************************************
*
      subroutine geometricFunc(dalpha,phase)
*
***********************************************************************************************************
      include 'vaba_param.inc'
      real*8  dalpha, phase
*     ================================================================================================     *       
      dalpha  = 2.d0 - 2.d0*phase
*     ================================================================================================     *        
      return 
      end subroutine geometricFunc  


