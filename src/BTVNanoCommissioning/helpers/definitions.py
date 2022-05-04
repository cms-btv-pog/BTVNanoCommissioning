def definitions(defaults=True):
      input_names = ['Jet_eta',
      'Jet_pt',
      'Jet_DeepCSV_flightDistance2dSig',
      'Jet_DeepCSV_flightDistance2dVal',
      'Jet_DeepCSV_flightDistance3dSig',
      'Jet_DeepCSV_flightDistance3dVal',
      'Jet_DeepCSV_trackDecayLenVal_0',
      'Jet_DeepCSV_trackDecayLenVal_1',
      'Jet_DeepCSV_trackDecayLenVal_2',
      'Jet_DeepCSV_trackDecayLenVal_3',
      'Jet_DeepCSV_trackDecayLenVal_4',
      'Jet_DeepCSV_trackDecayLenVal_5',
      'Jet_DeepCSV_trackDeltaR_0',
      'Jet_DeepCSV_trackDeltaR_1',
      'Jet_DeepCSV_trackDeltaR_2',
      'Jet_DeepCSV_trackDeltaR_3',
      'Jet_DeepCSV_trackDeltaR_4',
      'Jet_DeepCSV_trackDeltaR_5',
      'Jet_DeepCSV_trackEtaRel_0',
      'Jet_DeepCSV_trackEtaRel_1',
      'Jet_DeepCSV_trackEtaRel_2',
      'Jet_DeepCSV_trackEtaRel_3',
      'Jet_DeepCSV_trackJetDistVal_0',
      'Jet_DeepCSV_trackJetDistVal_1',
      'Jet_DeepCSV_trackJetDistVal_2',
      'Jet_DeepCSV_trackJetDistVal_3',
      'Jet_DeepCSV_trackJetDistVal_4',
      'Jet_DeepCSV_trackJetDistVal_5',
      'Jet_DeepCSV_trackJetPt',
      'Jet_DeepCSV_trackPtRatio_0',
      'Jet_DeepCSV_trackPtRatio_1',
      'Jet_DeepCSV_trackPtRatio_2',
      'Jet_DeepCSV_trackPtRatio_3',
      'Jet_DeepCSV_trackPtRatio_4',
      'Jet_DeepCSV_trackPtRatio_5',
      'Jet_DeepCSV_trackPtRel_0',
      'Jet_DeepCSV_trackPtRel_1',
      'Jet_DeepCSV_trackPtRel_2',
      'Jet_DeepCSV_trackPtRel_3',
      'Jet_DeepCSV_trackPtRel_4',
      'Jet_DeepCSV_trackPtRel_5',
      'Jet_DeepCSV_trackSip2dSigAboveCharm',
      'Jet_DeepCSV_trackSip2dSig_0',
      'Jet_DeepCSV_trackSip2dSig_1',
      'Jet_DeepCSV_trackSip2dSig_2',
      'Jet_DeepCSV_trackSip2dSig_3',
      'Jet_DeepCSV_trackSip2dSig_4',
      'Jet_DeepCSV_trackSip2dSig_5',
      'Jet_DeepCSV_trackSip2dValAboveCharm',
      'Jet_DeepCSV_trackSip3dSigAboveCharm',
      'Jet_DeepCSV_trackSip3dSig_0',
      'Jet_DeepCSV_trackSip3dSig_1',
      'Jet_DeepCSV_trackSip3dSig_2',
      'Jet_DeepCSV_trackSip3dSig_3',
      'Jet_DeepCSV_trackSip3dSig_4',
      'Jet_DeepCSV_trackSip3dSig_5',
      'Jet_DeepCSV_trackSip3dValAboveCharm',
      'Jet_DeepCSV_trackSumJetDeltaR',
      'Jet_DeepCSV_trackSumJetEtRatio',
      'Jet_DeepCSV_vertexCategory',
      'Jet_DeepCSV_vertexEnergyRatio',
      'Jet_DeepCSV_vertexJetDeltaR',
      'Jet_DeepCSV_vertexMass',
      'Jet_DeepCSV_jetNSecondaryVertices',
      'Jet_DeepCSV_jetNSelectedTracks',
      'Jet_DeepCSV_jetNTracksEtaRel',
      'Jet_DeepCSV_vertexNTracks',]


      display_names = ['Jet $\eta$',
                  'Jet $p_T$',
                  'Flight Distance 2D Sig','Flight Distance 2D Val','Flight Distance 3D Sig', 'Flight Distance 3D Val',
                  'Track Decay Len Val [0]','Track Decay Len Val [1]','Track Decay Len Val [2]','Track Decay Len Val [3]','Track Decay Len Val [4]','Track Decay Len Val [5]',
                  'Track $\Delta R$ [0]','Track $\Delta R$ [1]','Track $\Delta R$ [2]','Track $\Delta R$ [3]','Track $\Delta R$ [4]','Track $\Delta R$ [5]',
                  'Track $\eta_{rel}$ [0]','Track $\eta_{rel}$ [1]','Track $\eta_{rel}$ [2]','Track $\eta_{rel}$ [3]',
                  'Track Jet Dist Val [0]','Track Jet Dist Val [1]','Track Jet Dist Val [2]','Track Jet Dist Val [3]','Track Jet Dist Val [4]','Track Jet Dist Val [5]',
                  'Track Jet $p_T$',
                  'Track $p_{T,rel}$ Ratio [0]','Track $p_{T,rel}$ Ratio [1]','Track $p_{T,rel}$ Ratio [2]','Track $p_{T,rel}$ Ratio [3]','Track $p_{T,rel}$ Ratio [4]','Track $p_{T,rel}$ Ratio [5]',
                  'Track $p_{T,rel}$ [0]','Track $p_{T,rel}$ [1]','Track $p_{T,rel}$ [2]','Track $p_{T,rel}$ [3]','Track $p_{T,rel}$ [4]','Track $p_{T,rel}$ [5]',
                  'Track SIP 2D Sig Above Charm',
                  'Track SIP 2D Sig [0]','Track SIP 2D Sig [1]','Track SIP 2D Sig [2]','Track SIP 2D Sig [3]','Track SIP 2D Sig [4]','Track SIP 2D Sig [5]',
                  'Track SIP 2D Val Above Charm',
                  'Track SIP 3D Sig Above Charm',
                  'Track SIP 3D Sig [0]','Track SIP 3D Sig [1]','Track SIP 3D Sig [2]','Track SIP 3D Sig [3]','Track SIP 3D Sig [4]','Track SIP 3D Sig [5]',
                  'Track SIP 3D Val Above Charm',
                  'Track Sum Jet $\Delta R$','Track Sum Jet $E_T$ Ratio',
                  'Vertex Category','Vertex Energy Ratio','Vertex Jet $\Delta R$','Vertex Mass',
                  'Jet N Secondary Vertices','Jet N Selected Tracks','Jet N Tracks $\eta_{rel}$','Vertex N Tracks',]

      jetINDEX = [0,1,28,41,48,49,56,57,58,59,63,64,65] 
      trackINDEX = [6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40,42,43,44,45,46,47,50,51,52,53,54,55,]
      svINDEX = [2,3,4,5,60,61,62,66]

      interger_variables = [59,63,64,65,66]

      manual_ranges = [[None,None],
                  [0,250],
                  [-0.1,100],
                  [-0.1,2.6],
                  [-0.1,100],
                  [-0.1,5],
                  [-0.1,1],
                  [-0.1,1],
                  [-0.1,1],
                  [-0.1,1],
                  [-0.1,1],
                  [-0.1,1],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.,9],
                  [-0.,9],
                  [-0.,9],
                  [-0.,9],
                  [-0.08,0.0025],
                  [-0.08,0.0025],
                  [-0.08,0.0025],
                  [-0.08,0.0025],
                  [-0.08,0.0025],
                  [-0.08,0.0025],
                  [0,250],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.001,0.301],
                  [-0.1,3.1],
                  [-0.1,3.1],
                  [-0.1,3.1],
                  [-0.1,3.1],
                  [-0.1,3.1],
                  [-0.1,3.1],
                  [-5.5,5.5],
                  [-4.5,16],
                  [-5,13],
                  [-5.5,10],
                  [-6,7],
                  [-6.5,4.5],
                  [-7,2],
                  [-0.06,0.06],
                  [-6.5,6.5],
                  [-25,50],
                  [-25,50],
                  [-25,50],
                  [-25,50],
                  [-25,50],
                  [-25,50],
                  [-0.06,0.06],
                  [-0.001,0.301],
                  [None,1.4],
                  [-0.6,2.6],
                  [0,2.5],
                  [-0.001,0.301],
                  [0,20],
                  [-0.5,None],
                  [-0.5,None],
                  [-0.5,None],
                  [-0.5,None],]

      bins=[25,#Jet_eta
            50,#Jet_pt
            101,#Jet_DeepCSV_flightDistance2dSig
            27,#Jet_DeepCSV_flightDistance2dVal
      100,#Jet_DeepCSV_flightDistance3dSig
      51,#Jet_DeepCSV_flightDistance3dVal
      22,#Jet_DeepCSV_trackDecayLenVal_0
      22,#Jet_DeepCSV_trackDecayLenVal_1
      22,#Jet_DeepCSV_trackDecayLenVal_2
      22,#Jet_DeepCSV_trackDecayLenVal_3
      22,#Jet_DeepCSV_trackDecayLenVal_4
      22,#Jet_DeepCSV_trackDecayLenVal_5
      30,#Jet_DeepCSV_trackDeltaR_0
      30,#Jet_DeepCSV_trackDeltaR_1
      30,#Jet_DeepCSV_trackDeltaR_2
      30,#Jet_DeepCSV_trackDeltaR_3
      30,#Jet_DeepCSV_trackDeltaR_4
      30,#Jet_DeepCSV_trackDeltaR_5
      30,#Jet_DeepCSV_trackEtaRel_0
      30,#Jet_DeepCSV_trackEtaRel_1
      30,#Jet_DeepCSV_trackEtaRel_2
      30,#Jet_DeepCSV_trackEtaRel_3
      35,#Jet_DeepCSV_trackJetDistVal_0
      35,#Jet_DeepCSV_trackJetDistVal_1
      35,#Jet_DeepCSV_trackJetDistVal_2
      35,#Jet_DeepCSV_trackJetDistVal_3
      35,#Jet_DeepCSV_trackJetDistVal_4
      35,#Jet_DeepCSV_trackJetDistVal_5
      50,#Jet_DeepCSV_trackJetPt
      30,#Jet_DeepCSV_trackPtRatio_0
      30,#Jet_DeepCSV_trackPtRatio_1
      30,#Jet_DeepCSV_trackPtRatio_2
      30,#Jet_DeepCSV_trackPtRatio_3
      30,#Jet_DeepCSV_trackPtRatio_4
      30,#Jet_DeepCSV_trackPtRatio_5
      32,#Jet_DeepCSV_trackPtRel_0
      32,#Jet_DeepCSV_trackPtRel_1
      32,#Jet_DeepCSV_trackPtRel_2
      32,#Jet_DeepCSV_trackPtRel_3
      32,#Jet_DeepCSV_trackPtRel_4
      32,#Jet_DeepCSV_trackPtRel_5
      22,#Jet_DeepCSV_trackSip2dSigAboveCharm
      21,#Jet_DeepCSV_trackSip2dSig_0
      18,#Jet_DeepCSV_trackSip2dSig_1
      16,#Jet_DeepCSV_trackSip2dSig_2
      26,#Jet_DeepCSV_trackSip2dSig_3
      22,#Jet_DeepCSV_trackSip2dSig_4
      18,#Jet_DeepCSV_trackSip2dSig_5
      24,#Jet_DeepCSV_trackSip2dValAboveCharm
      26,#Jet_DeepCSV_trackSip3dSigAboveCharm
      25,#Jet_DeepCSV_trackSip3dSig_0
      25,#Jet_DeepCSV_trackSip3dSig_1
      25,#Jet_DeepCSV_trackSip3dSig_2
      25,#Jet_DeepCSV_trackSip3dSig_3
      25,#Jet_DeepCSV_trackSip3dSig_4
      25,#Jet_DeepCSV_trackSip3dSig_5
      24,#Jet_DeepCSV_trackSip3dValAboveCharm
      30,#Jet_DeepCSV_trackSumJetDeltaR
      28,#Jet_DeepCSV_trackSumJetEtRatio
      26,#Jet_DeepCSV_vertexCategory
      25,#Jet_DeepCSV_vertexEnergyRatio
      30,#Jet_DeepCSV_vertexJetDeltaR
      20,#Jet_DeepCSV_vertexMass
      25,#Jet_DeepCSV_jetNSecondaryVertices
      25,#Jet_DeepCSV_jetNSelectedTracks
      25,#Jet_DeepCSV_jetNTracksEtaRel
      25,#Jet_DeepCSV_vertexNTracks
      ]

      inputVar_units = ['',
      'GeV',
      '',
      'cm',
      '',
      'cm',
      'cm',
      'cm',
      'cm',
      'cm',
      'cm',
      'cm',
      '',
      '',
      '',
      '',
      '',
      '',
      '',
      '',
      '',
      '',
      'cm',
      'cm',
      'cm',
      'cm',
      'cm',
      'cm',
      'GeV',
      '',
      '',
      '',
      '',
      '',
      '',
      'GeV',
      'GeV',
      'GeV',
      'GeV',
      'GeV',
      'GeV',
      '',
      '',
      '',
      '',
      '',
      '',
      '',
      'cm',
      '',
      '',
      '',
      '',
      '',
      '',
      '',
      'cm',
      '',
      '',
      '',
      '',
      '',
      'GeV',
      '',
      '',
      '',
      '']

      format_unit = ['2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f',
      '2f']

      format_unit_digits = [2, #
      0,  #
      1,  # 'Jet_DeepCSV_flightDist
      3,  # 'Jet_DeepCSV_flightDist
      1,  # 'Jet_DeepCSV_flightDist
      2,  # 'Jet_DeepCSV_flightDist
      3,  # 'Jet_DeepCSV_trackDecay
      3,  # 'Jet_DeepCSV_trackDecay
      3,  # 'Jet_DeepCSV_trackDecay
      3,  # 'Jet_DeepCSV_trackDecay
      3,  # 'Jet_DeepCSV_trackDecay
      3,  # 'Jet_DeepCSV_trackDecay
      3,  # 'Jet_DeepCSV_trackDelta
      3,  # 'Jet_DeepCSV_trackDelta
      3,  # 'Jet_DeepCSV_trackDelta
      3,  # 'Jet_DeepCSV_trackDelta
      3,  # 'Jet_DeepCSV_trackDelta
      3,  # 'Jet_DeepCSV_trackDelta
      2,  # 'Jet_DeepCSV_trackEtaRe
      2,  # 'Jet_DeepCSV_trackEtaRe
      2,  # 'Jet_DeepCSV_trackEtaRe
      2,  # 'Jet_DeepCSV_trackEtaRe
      4,  # 'Jet_DeepCSV_trackJetDi
      4,  # 'Jet_DeepCSV_trackJetDi
      4,  # 'Jet_DeepCSV_trackJetDi
      4,  # 'Jet_DeepCSV_trackJetDi
      4,  # 'Jet_DeepCSV_trackJetDi
      4,  # 'Jet_DeepCSV_trackJetDi
      1,  # 'Jet_DeepCSV_trackJetPt
      3,  # 'Jet_DeepCSV_trackPtRat
      3,  # 'Jet_DeepCSV_trackPtRat
      3,  # 'Jet_DeepCSV_trackPtRat
      3,  # 'Jet_DeepCSV_trackPtRat
      3,  # 'Jet_DeepCSV_trackPtRat
      3,  # 'Jet_DeepCSV_trackPtRat
      2,  # 'Jet_DeepCSV_trackPtRel
      2,  # 'Jet_DeepCSV_trackPtRel
      2,  # 'Jet_DeepCSV_trackPtRel
      2,  # 'Jet_DeepCSV_trackPtRel
      2,  # 'Jet_DeepCSV_trackPtRel
      2,  # 'Jet_DeepCSV_trackPtRel
      2,  # 'Jet_DeepCSV_trackSip2d
      2,  # 'Jet_DeepCSV_trackSip2d
      2,  # 'Jet_DeepCSV_trackSip2d
      2,  # 'Jet_DeepCSV_trackSip2d
      2,  # 'Jet_DeepCSV_trackSip2d
      2,  # 'Jet_DeepCSV_trackSip2d
      2,  # 'Jet_DeepCSV_trackSip2d
      3,  # 'Jet_DeepCSV_trackSip2d
      2,  # 'Jet_DeepCSV_trackSip3d
      2,  # 'Jet_DeepCSV_trackSip3d
      2,  # 'Jet_DeepCSV_trackSip3d
      2,  # 'Jet_DeepCSV_trackSip3d
      2,  # 'Jet_DeepCSV_trackSip3d
      2,  # 'Jet_DeepCSV_trackSip3d
      2,  # 'Jet_DeepCSV_trackSip3d
      3,  # 'Jet_DeepCSV_trackSip3d
      3,  # 'Jet_DeepCSV_trackSumJe
      2,  # 'Jet_DeepCSV_trackSumJe
      2,  # 'Jet_DeepCSV_vertexCate
      2,  # 'Jet_DeepCSV_vertexEner
      3,  # 'Jet_DeepCSV_vertexJetD
      2,  # 'Jet_DeepCSV_vertexMass
      2,  # 'Jet_DeepCSV_jetNSecond
      2,  # 'Jet_DeepCSV_jetNSelect
      2,  # 'Jet_DeepCSV_jetNTracks
      2]  # 'Jet_DeepCSV_vertexNTra

      ylabel_text = ['Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Jets',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Tracks',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets',
      'Jets']
      if defaults:return input_names,manual_ranges,bins
