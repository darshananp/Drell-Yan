#ifndef SAMODULESIMPLEDIMUON_H_
#define SAMODULESIMPLEDIMUON_H_

#include <vector>
#include <map>
#include <string>

#include <SubsysReco.h>
#include <DiMuonContainer.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <DiMuon.h>
#include <FvtxPrimVtxContainer.h>
#include <SingleMuonContainer.h>
#include <SingleMuon.h>
#include <DSTReader_v1.h>
#include <PHPythiaHeader.h>
#include <PHPythiaContainerV2.h>
#include <PHPythia.h>
#include "saModuleBase.h"


class saModuleSimpleDimuon : public saModuleBase{

public:
  saModuleSimpleDimuon(const std::string &name = "SimpleDiMuon");
  virtual
  ~saModuleSimpleDimuon();

	bool _use_2_hit_tracklet;
	bool _use_cut_tracklet_chi2;
	float _cut_tracklet_chi2;
	float _cut_tracklet_dcar;

  bool passesCuts() const;

  int GetNumberofTracklets(DSTReader *fvtx_trk_map, const DiMuon *dimuon);

#ifndef __CINT__

  //! global initialization
  virtual int
  init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! Run initialization
  virtual int
  init_run(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! event method
  virtual int
  event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

  //! global termination
  virtual int
  end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm);

#endif


  saHist * _h_mass;

private:

  DSTReader *fvtx_trk_map;
  DiMuonContainer *dimuons;
  SingleMuonContainer *singlemuoncontainer;
  DiMuon* dimuon;
  int nTracklets;


};

#endif /* SAMODULESIMPLEDIMUON_H_ */
