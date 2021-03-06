

class hist_class {
 public:
  hist_class(TString the_desc, bool is_it_data);
  void Delete();
  void NormalizeHists();
  void Write(int mc_type_i);
  void AddHists(hist_class *more_hists, float wt);
  bool is_data;
  TString desc;
  int n_evt_raw ;

  TH1F* NEvents;
  TH1F* NEvents_test;
  TH1F* NEvents_after_noise;
  TH1F* NEvents_dijets;
  TH1F* Vz; 
  TH1F* Centrality;
  TH1F* Vz_new;
  TH1F* Centrality_new;
  TH1F* dPhi_hist[nCBins];
  TH1F* Aj[nCBins];

  
  TH1F* all_jets_corrpT[nCBins][nPtBins];
  TH1F* all_jets_phi[nCBins][nPtBins];
  TH1F* all_jets_eta[nCBins][nPtBins];
  TH1F* only_leadingjets_corrpT[nCBins][nPtBins];
  TH1F* only_leadingjets_phi[nCBins][nPtBins];
  TH1F* only_leadingjets_eta[nCBins][nPtBins];

  TH1F* only_subleadingjets_corrpT[nCBins][nPtBins];
  TH1F* only_subleadingjets_phi[nCBins][nPtBins];
  TH1F* only_subleadingjets_eta[nCBins][nPtBins];
  TH1F* only_nonleadingjets_corrpT[nCBins][nPtBins];
  TH1F* only_nonleadingjets_phi[nCBins][nPtBins];
  TH1F* only_nonleadingjets_eta[nCBins][nPtBins];

 
  TH1F* TrkPhi[nCBins][nPtBins][nTrkPtBins];
  TH1F* TrkPt[nCBins][nPtBins][nTrkPtBins];
  TH1F* TrkEta[nCBins][nPtBins][nTrkPtBins];

  TH1F* TrkPhi_weighted[nCBins][nPtBins][nTrkPtBins];
  TH1F* TrkPt_weighted[nCBins][nPtBins][nTrkPtBins];
  TH1F* TrkEta_weighted[nCBins][nPtBins][nTrkPtBins];
 

  TH1F* ME_TrkPhi[nCBins][nPtBins][nTrkPtBins];
  TH1F* ME_TrkPt[nCBins][nPtBins][nTrkPtBins];
  TH1F* ME_TrkEta[nCBins][nPtBins][nTrkPtBins];

  TH1F* ME_TrkPhi_weighted[nCBins][nPtBins][nTrkPtBins];
  TH1F* ME_TrkPt_weighted[nCBins][nPtBins][nTrkPtBins];
  TH1F* ME_TrkEta_weighted[nCBins][nPtBins][nTrkPtBins];
 


  TH2D* hJetTrackSignalBackground[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackground_notrkcorr[nCBins][nPtBins][nTrkPtBins];
 
  TH2D* hJetTrackSignalBackgroundLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundLeading_notrkcorr[nCBins][nPtBins][nTrkPtBins];
    
  TH2D* hJetTrackSignalBackgroundSubLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundSubLeading_notrkcorr[nCBins][nPtBins][nTrkPtBins];
   
   
  TH2D* hJetTrackSignalBackgroundNonLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackSignalBackgroundNonLeading_notrkcorr[nCBins][nPtBins][nTrkPtBins];
   


  TH2D* hJetTrackME[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackME_notrkcorr[nCBins][nPtBins][nTrkPtBins];
 
  TH2D* hJetTrackMELeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackMELeading_notrkcorr[nCBins][nPtBins][nTrkPtBins];
    
  TH2D* hJetTrackMESubLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackMESubLeading_notrkcorr[nCBins][nPtBins][nTrkPtBins];
   
  
  TH2D* hJetTrackMENonLeading[nCBins][nPtBins][nTrkPtBins];
  TH2D* hJetTrackMENonLeading_notrkcorr[nCBins][nPtBins][nTrkPtBins];
  
};


hist_class::hist_class(TString the_desc, bool is_it_data) {
  n_evt_raw = 0;
  desc = the_desc;
  is_data = is_it_data;

  NEvents = new TH1F((TString) (desc + "_Nevents"), "", 100, 0., 100.);     NEvents->Sumw2(); 
  NEvents_test = new TH1F((TString) (desc + "_Nevents_test"), "", 100, 0., 100.);     NEvents_test->Sumw2();
  NEvents_after_noise = new TH1F((TString) (desc + "_Nevents_after_noise"), "", 100, 0., 100.);     NEvents_after_noise->Sumw2();
  NEvents_dijets = new TH1F((TString) (desc + "_Nevents_dijets"), "", 100, 0., 100.);     NEvents_dijets->Sumw2();
  Centrality = new TH1F((TString) (desc + "_Centrality"), "", 40,0.,200);     Centrality->Sumw2();
  Vz = new TH1F((TString) (desc + "_Vz"), "", 80, -20., 20.); Vz->Sumw2();
  Centrality_new = new TH1F((TString) (desc + "_Centrality_Reweighted"), "", 40,0.,200);     Centrality_new->Sumw2();
  Vz_new = new TH1F((TString) (desc + "_Vz_Reweighted"), "", 80, -20., 20.); Vz_new->Sumw2();

  
  for (int ibin=0;ibin<nCBins;ibin++){
 
     dPhi_hist[ibin] = new TH1F((TString) (desc + "_dPhi_hist"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 30, 0.,TMath::Pi());     dPhi_hist[ibin]->Sumw2();
    Aj[ibin] = new TH1F((TString) (desc + "_Aj"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1]), "", 25,0.,1.25);     Aj[ibin]->Sumw2();

    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
      
      all_jets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 50, 0., 500.);     all_jets_corrpT[ibin][ibin2]->Sumw2();
      all_jets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 72, -TMath::Pi(), TMath::Pi());     all_jets_phi[ibin][ibin2]->Sumw2();

      all_jets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_all_jets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, -5., 5.);     all_jets_eta[ibin][ibin2]->Sumw2();

      //// leading jet histograms
      only_leadingjets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     only_leadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_leadingjets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, -TMath::Pi(),TMath::Pi());     only_leadingjets_phi[ibin][ibin2]->Sumw2();
      only_leadingjets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_only_leadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, -5., 5.);     only_leadingjets_eta[ibin][ibin2]->Sumw2();


      //// subleading jet histograms
      only_subleadingjets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     only_subleadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_subleadingjets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, -TMath::Pi(),TMath::Pi());     only_subleadingjets_phi[ibin][ibin2]->Sumw2();
      only_subleadingjets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_only_subleadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, -5., 5.);     only_subleadingjets_eta[ibin][ibin2]->Sumw2();



 only_nonleadingjets_corrpT[ibin][ibin2] = new TH1F((TString) (desc + "_only_nonleadingjets_corrpT"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, 0., 500.);     only_nonleadingjets_corrpT[ibin][ibin2]->Sumw2();
      only_nonleadingjets_phi[ibin][ibin2] = new TH1F((TString) (desc + "_only_nonleadingjets_phi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 36, -TMath::Pi(),TMath::Pi());     only_nonleadingjets_phi[ibin][ibin2]->Sumw2();
      only_nonleadingjets_eta[ibin][ibin2] = new TH1F((TString) (desc + "_only_nonleadingjets_eta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]), "", 100, -5., 5.);     only_nonleadingjets_eta[ibin][ibin2]->Sumw2();






      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){
     
	hJetTrackSignalBackground[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackground"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-TMath::Pi()/2,3*TMath::Pi()/2);     hJetTrackSignalBackground[ibin][ibin2][ibin3]->Sumw2();

    
	hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackground_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Sumw2();


	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackgroundLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Sumw2();

    	hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackgroundLeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Sumw2();


	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackgroundSubLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Sumw2();

    
	hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackgroundSubLeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Sumw2();



	hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackgroundNonLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3]->Sumw2();

    
	hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackSignalBackgroundNonLeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Sumw2();






	TrkPt[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 500, 0., 20.);     TrkPt[ibin][ibin2][ibin3]->Sumw2();

	TrkEta[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkEta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, -5., 5. );     TrkEta[ibin][ibin2][ibin3]->Sumw2();

	TrkPhi[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPhi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi());     TrkPhi[ibin][ibin2][ibin3]->Sumw2();



	TrkPt_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPt_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 500, 0., 20.);     TrkPt_weighted[ibin][ibin2][ibin3]->Sumw2();

	TrkEta_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkEta_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, -5., 5.);     TrkEta_weighted[ibin][ibin2][ibin3]->Sumw2();

	TrkPhi_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_TrkPhi_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "",100, -0.5*TMath::Pi(), 1.5*TMath::Pi());     TrkPhi_weighted[ibin][ibin2][ibin3]->Sumw2();





	ME_TrkPt[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_ME_TrkPt"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, 0., 100.);     ME_TrkPt[ibin][ibin2][ibin3]->Sumw2();

	ME_TrkEta[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_ME_TrkEta"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, -5., 5. );     ME_TrkEta[ibin][ibin2][ibin3]->Sumw2();

	ME_TrkPhi[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_ME_TrkPhi"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, -0.5*TMath::Pi(), 1.5*TMath::Pi());     ME_TrkPhi[ibin][ibin2][ibin3]->Sumw2();



	ME_TrkPt_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_ME_TrkPt_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, 0., 100.);     ME_TrkPt_weighted[ibin][ibin2][ibin3]->Sumw2();

	ME_TrkEta_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_ME_TrkEta_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100, -5., 5.);     ME_TrkEta_weighted[ibin][ibin2][ibin3]->Sumw2();

	ME_TrkPhi_weighted[ibin][ibin2][ibin3] = new TH1F((TString) (desc + "_ME_TrkPhi_weighted"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "",100, -0.5*TMath::Pi(), 1.5*TMath::Pi());     ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Sumw2();






	hJetTrackME[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackME"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-TMath::Pi()/2,3*TMath::Pi()/2);     hJetTrackME[ibin][ibin2][ibin3]->Sumw2();

    
	hJetTrackME_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackME_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Sumw2();


	hJetTrackMELeading[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackMELeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackMELeading[ibin][ibin2][ibin3]->Sumw2();

    	hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackMELeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Sumw2();


	hJetTrackMESubLeading[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackMESubLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackMESubLeading[ibin][ibin2][ibin3]->Sumw2();

    
	hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackMESubLeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Sumw2();



	hJetTrackMENonLeading[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackMENonLeading"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackMENonLeading[ibin][ibin2][ibin3]->Sumw2();

    
	hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3] = new TH2D((TString) (desc + "_hJetTrackMENonLeading_notrkcorr"+ CBin_strs[ibin] + "_" + CBin_strs[ibin+1] + "_" + PtBin_strs[ibin2] + "_" + PtBin_strs[ibin2+1]+ "_" + TrkPtBin_strs[ibin3] + "_" + TrkPtBin_strs[ibin3+1]), "", 100,-5,5,100,-0.5*TMath::Pi(),1.5*TMath::Pi());     hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Sumw2();

     
     
      } /// ibin3
    } // pt bin loop
  } // centrality bin loop
} // hist class loop




void hist_class::AddHists(hist_class *more_hists, float wt)
{


  
  NEvents->Add(more_hists->NEvents, wt);
  NEvents_test->Add(more_hists->NEvents_test, wt);
  NEvents_after_noise->Add(more_hists->NEvents_after_noise, wt);
  NEvents_dijets->Add(more_hists->NEvents_dijets, wt);


  Centrality->Add(more_hists->Centrality, wt);
  Vz->Add(more_hists->Vz, wt);

  Centrality_new->Add(more_hists->Centrality_new, wt);
  Vz_new->Add(more_hists->Vz_new, wt);


  for (int ibin=0;ibin<nCBins;ibin++){
    
    dPhi_hist[ibin]->Add(more_hists->dPhi_hist[ibin],wt);
    Aj[ibin]->Add(more_hists->Aj[ibin],wt);


    for (int ibin2=0;ibin2<nPtBins;ibin2++){ 
    
      all_jets_corrpT[ibin][ibin2]->Add(more_hists->all_jets_corrpT[ibin][ibin2], wt);
      all_jets_phi[ibin][ibin2]->Add(more_hists->all_jets_phi[ibin][ibin2], wt);
      all_jets_eta[ibin][ibin2]->Add(more_hists->all_jets_eta[ibin][ibin2], wt);

      //// leading jet histograms
      only_leadingjets_corrpT[ibin][ibin2]->Add(more_hists->only_leadingjets_corrpT[ibin][ibin2], wt);
      only_leadingjets_phi[ibin][ibin2]->Add(more_hists->only_leadingjets_phi[ibin][ibin2], wt);
      only_leadingjets_eta[ibin][ibin2]->Add(more_hists->only_leadingjets_eta[ibin][ibin2], wt);
     
      //// subleading jet histograms
      only_subleadingjets_corrpT[ibin][ibin2]->Add(more_hists->only_subleadingjets_corrpT[ibin][ibin2], wt);
      only_subleadingjets_phi[ibin][ibin2]->Add(more_hists->only_subleadingjets_phi[ibin][ibin2], wt);
      only_subleadingjets_eta[ibin][ibin2]->Add(more_hists->only_subleadingjets_eta[ibin][ibin2], wt);
     

      //// nonleading jet histograms
      only_nonleadingjets_corrpT[ibin][ibin2]->Add(more_hists->only_nonleadingjets_corrpT[ibin][ibin2], wt);
      only_nonleadingjets_phi[ibin][ibin2]->Add(more_hists->only_nonleadingjets_phi[ibin][ibin2], wt);
      only_nonleadingjets_eta[ibin][ibin2]->Add(more_hists->only_nonleadingjets_eta[ibin][ibin2], wt);
     


 
      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){


	TrkPt[ibin][ibin2][ibin3]->Add(more_hists->TrkPt[ibin][ibin2][ibin3], wt);
	TrkEta[ibin][ibin2][ibin3]->Add(more_hists->TrkEta[ibin][ibin2][ibin3], wt);
	TrkPhi[ibin][ibin2][ibin3]->Add(more_hists->TrkPhi[ibin][ibin2][ibin3], wt);

	TrkPt_weighted[ibin][ibin2][ibin3]->Add(more_hists->TrkPt_weighted[ibin][ibin2][ibin3], wt);
	TrkEta_weighted[ibin][ibin2][ibin3]->Add(more_hists->TrkEta_weighted[ibin][ibin2][ibin3], wt);
	TrkPhi_weighted[ibin][ibin2][ibin3]->Add(more_hists->TrkPhi_weighted[ibin][ibin2][ibin3], wt);


	ME_TrkPt[ibin][ibin2][ibin3]->Add(more_hists->ME_TrkPt[ibin][ibin2][ibin3], wt);
	ME_TrkEta[ibin][ibin2][ibin3]->Add(more_hists->ME_TrkEta[ibin][ibin2][ibin3], wt);
	ME_TrkPhi[ibin][ibin2][ibin3]->Add(more_hists->ME_TrkPhi[ibin][ibin2][ibin3], wt);

	ME_TrkPt_weighted[ibin][ibin2][ibin3]->Add(more_hists->ME_TrkPt_weighted[ibin][ibin2][ibin3], wt);
	ME_TrkEta_weighted[ibin][ibin2][ibin3]->Add(more_hists->ME_TrkEta_weighted[ibin][ibin2][ibin3], wt);
	ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Add(more_hists->ME_TrkPhi_weighted[ibin][ibin2][ibin3], wt);




	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackground[ibin][ibin2][ibin3], wt);
	hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3], wt);

	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3], wt);
	hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3], wt);

	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3], wt);
	hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3], wt);

	hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3], wt);
	hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3], wt);




	hJetTrackME[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackME[ibin][ibin2][ibin3], wt);
	hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackME_notrkcorr[ibin][ibin2][ibin3], wt);

	hJetTrackMELeading[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackMELeading[ibin][ibin2][ibin3], wt);
	hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3], wt);

	hJetTrackMESubLeading[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackMESubLeading[ibin][ibin2][ibin3], wt);
	hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3], wt);


	hJetTrackMENonLeading[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackMENonLeading[ibin][ibin2][ibin3], wt);
	hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Add(more_hists->hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3], wt);


      } /// ibin3
    }
  }
}


void hist_class::Delete()
{
  delete NEvents;
  delete NEvents_test;
  delete NEvents_after_noise;
  delete NEvents_dijets;
  delete Centrality;
  delete Vz;
  delete Centrality_new;
  delete Vz_new;

  for (int ibin=0;ibin<nCBins;ibin++){

    delete dPhi_hist[ibin];
    delete Aj[ibin];

    for (int ibin2=0;ibin2<nPtBins;ibin2++){

      delete all_jets_corrpT[ibin][ibin2];
      delete all_jets_phi[ibin][ibin2];
      delete all_jets_eta[ibin][ibin2];

      delete only_leadingjets_corrpT[ibin][ibin2];
      delete only_leadingjets_phi[ibin][ibin2];
      delete only_leadingjets_eta[ibin][ibin2];
   
      delete only_subleadingjets_corrpT[ibin][ibin2];
      delete only_subleadingjets_phi[ibin][ibin2];
      delete only_subleadingjets_eta[ibin][ibin2];


      delete only_nonleadingjets_corrpT[ibin][ibin2];
      delete only_nonleadingjets_phi[ibin][ibin2];
      delete only_nonleadingjets_eta[ibin][ibin2];


      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	delete hJetTrackSignalBackground[ibin][ibin2][ibin3];
	delete hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3];

	delete hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3];
	delete hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3];

	delete hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3];
	delete hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3];


	delete hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3];
	delete hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3];


	delete hJetTrackME[ibin][ibin2][ibin3];
	delete hJetTrackME_notrkcorr[ibin][ibin2][ibin3];

	delete hJetTrackMELeading[ibin][ibin2][ibin3];
	delete hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3];

	delete hJetTrackMESubLeading[ibin][ibin2][ibin3];
	delete hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3];

	delete hJetTrackMENonLeading[ibin][ibin2][ibin3];
	delete hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3];



	delete TrkPt[ibin][ibin2][ibin3];
	delete TrkEta[ibin][ibin2][ibin3];
	delete TrkPhi[ibin][ibin2][ibin3];

	delete TrkPt_weighted[ibin][ibin2][ibin3];
	delete TrkEta_weighted[ibin][ibin2][ibin3];
	delete TrkPhi_weighted[ibin][ibin2][ibin3];




	delete ME_TrkPt[ibin][ibin2][ibin3];
	delete ME_TrkEta[ibin][ibin2][ibin3];
	delete ME_TrkPhi[ibin][ibin2][ibin3];

	delete ME_TrkPt_weighted[ibin][ibin2][ibin3];
	delete ME_TrkEta_weighted[ibin][ibin2][ibin3];
	delete ME_TrkPhi_weighted[ibin][ibin2][ibin3];

      } /// ibin3
    } // ibin2
  } // ibin
}









void hist_class::Write(int mc_type_i)
{

  TString parti_str = "";
  if( parti >= 0 ) {
    parti_str += "part";
    parti_str +=  parti;
  }

  TString pT_str = "";
  if( trkPtCut >= 0.49 && trkPtCut < 1.5 ) pT_str = "trkPtCut1";
  else if( trkPtCut >= 1.5 && trkPtCut < 2.5 ) pT_str = "trkPtCut2";
  else if( trkPtCut >= 2.5 && trkPtCut < 3.5 ) pT_str = "trkPtCut3";
  else if( trkPtCut >= 3.5 && trkPtCut < 4.5 ) pT_str = "trkPtCut4";
  else assert(0);  

  TString out_name = (TString) ("root_output/" + dataset_type_strs[dataset_type_code] + "_" + data_mc_type_strs[mc_type_i]+"_"+ parti_str + ".root");
  TFile *out_file = new TFile(out_name, "RECREATE");

  NEvents->Write();
  NEvents_test->Write();
  NEvents_after_noise->Write();
  NEvents_dijets->Write();

  Vz->Write();
  Centrality->Write();
  Vz_new->Write();
  Centrality_new->Write();

  for (int ibin=0;ibin<nCBins;ibin++){
    
    dPhi_hist[ibin]->Write();
    Aj[ibin]->Write();

    for (int ibin2=0;ibin2<nPtBins;ibin2++){

      all_jets_corrpT[ibin][ibin2]->Write();
      all_jets_phi[ibin][ibin2]->Write();
      all_jets_eta[ibin][ibin2]->Write();
      only_leadingjets_corrpT[ibin][ibin2]->Write();
      only_leadingjets_phi[ibin][ibin2]->Write();
      only_leadingjets_eta[ibin][ibin2]->Write();
      only_subleadingjets_corrpT[ibin][ibin2]->Write();
      only_subleadingjets_phi[ibin][ibin2]->Write();
      only_subleadingjets_eta[ibin][ibin2]->Write();

      only_nonleadingjets_corrpT[ibin][ibin2]->Write();
      only_nonleadingjets_phi[ibin][ibin2]->Write();
      only_nonleadingjets_eta[ibin][ibin2]->Write();

      for (int ibin3=0;ibin3<nTrkPtBins;ibin3++){

	hJetTrackSignalBackground[ibin][ibin2][ibin3]->Write();
	hJetTrackSignalBackground_notrkcorr[ibin][ibin2][ibin3]->Write();


	hJetTrackSignalBackgroundLeading[ibin][ibin2][ibin3]->Write();
	hJetTrackSignalBackgroundLeading_notrkcorr[ibin][ibin2][ibin3]->Write();


	hJetTrackSignalBackgroundSubLeading[ibin][ibin2][ibin3]->Write();
	hJetTrackSignalBackgroundSubLeading_notrkcorr[ibin][ibin2][ibin3]->Write();


	hJetTrackSignalBackgroundNonLeading[ibin][ibin2][ibin3]->Write();
	hJetTrackSignalBackgroundNonLeading_notrkcorr[ibin][ibin2][ibin3]->Write();


	hJetTrackME[ibin][ibin2][ibin3]->Write();
	hJetTrackME_notrkcorr[ibin][ibin2][ibin3]->Write();


	hJetTrackMELeading[ibin][ibin2][ibin3]->Write();
	hJetTrackMELeading_notrkcorr[ibin][ibin2][ibin3]->Write();


	hJetTrackMESubLeading[ibin][ibin2][ibin3]->Write();
	hJetTrackMESubLeading_notrkcorr[ibin][ibin2][ibin3]->Write();



	hJetTrackMENonLeading[ibin][ibin2][ibin3]->Write();
	hJetTrackMENonLeading_notrkcorr[ibin][ibin2][ibin3]->Write();




	TrkPt[ibin][ibin2][ibin3]->Write();
	TrkEta[ibin][ibin2][ibin3]->Write();
	TrkPhi[ibin][ibin2][ibin3]->Write();

	TrkPt_weighted[ibin][ibin2][ibin3]->Write();
	TrkEta_weighted[ibin][ibin2][ibin3]->Write();
	TrkPhi_weighted[ibin][ibin2][ibin3]->Write();


	ME_TrkPt[ibin][ibin2][ibin3]->Write();
	ME_TrkEta[ibin][ibin2][ibin3]->Write();
	ME_TrkPhi[ibin][ibin2][ibin3]->Write();

	ME_TrkPt_weighted[ibin][ibin2][ibin3]->Write();
	ME_TrkEta_weighted[ibin][ibin2][ibin3]->Write();
	ME_TrkPhi_weighted[ibin][ibin2][ibin3]->Write();


      } /// ibin3
    } /// ptbin
  }  //centralitybin
  out_file->Close();
} 


