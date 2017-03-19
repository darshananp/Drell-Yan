#ifndef __DNPFILTER_H__
#define __DNPFILTER_H__

#include <vector>
#include <map>
#include <string>

#include <SubsysReco.h>

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
class TClonesArray;

class VtxOut;
class TFvtxCompactTrk;
class TFvtxCompactTrkMap;
class FvtxPrimVtxContainer;
class FvtxPrimVtxContainer_v2;
class DiMuonContainer;
class DiMuonContainer_v9;



#include <iostream>
#include <fstream>

using namespace std;

class dnpfilter : public SubsysReco
{
	public:
		//! default constructor
		dnpfilter(const std::string &out_name = "femtoDST");

		//! default destructor
		virtual ~dnpfilter();

		int Init(PHCompositeNode *topNode);

		//! event method
		int process_event(PHCompositeNode *topNode);

  		int End(PHCompositeNode *topNode);

		void init_tree();
		void init_histos();
		bool set_dimuon_vars(const DiMuon *dimuon);
		bool pass_dimuon_cut(const DiMuon *dimuon) const;

		bool  momentum_matching(const DiMuon *dimuon);
		float set_singlemuon_vars(unsigned int i) const;
		int set_singlemuon_charge(unsigned int index) const;

		bool fill_dimuon_hist(const DiMuon *dimuon);
		float set_singlemuon_x0(unsigned int index) const;
		float set_singlemuon_y0(unsigned int index) const;
		float set_singlemuon_fvtxmutrz(unsigned int index) const;
		float set_singlemuon_fvtxmutrx(unsigned int index) const;
		float set_singlemuon_fvtxmutry(unsigned int index) const;
		int set_processtype(PHPythiaContainer *phpythia);
		int set_processtypetwo(PHPythiaContainer *phpythia);
		int GetNumberofTracklets(DSTReader *fvtx_trk_map, const DiMuon *dimuon);
		int GetNumberofTracklets_all(DSTReader *fvtx_trk_map, const DiMuon *dimuon);
		int GetNumberofTracklets_samearm(DSTReader *fvtx_trk_map, const DiMuon *dimuon);
		int GetNumberofTracklets_samearm_all(DSTReader *fvtx_trk_map, const DiMuon *dimuon);
		void set_is_sim(const bool a) {_is_sim = a;}
		void set_is_mb(const bool b) {_is_mb = b;}
		void set_out_name(const std::string &name)
		{
			_out_name = name;
		}
		//!
		std::string get_out_name() const
		{
			return _out_name;
		}
		enum DIMUON_CHARGE_TYPE {Unlike_Sign, Like_Sign};
		bool _use_2_hit_tracklet;
		bool _use_cut_tracklet_chi2;
		float _cut_tracklet_chi2;
		float _cut_tracklet_dcar;

	private:
		bool _is_sim;
		bool _is_mb;
		int nevents;	
		TTree *T;   
		float Evt_vtxchi2;     
		float Evt_vtxoor;
		float mass;
		float mass_fvtx;
		float mass_fvtxmutr;
		float pT;
		float pz;
		float rapidity;
		int charge;
		unsigned int lvl1_trigscaled;
		unsigned int lvl1scaledown;

		bool same_event;
		float Tr0_fvtx_dca_z;
		float Tr1_fvtx_dca_z;
		short Tr0_idhits;
		short Tr1_idhits;
		short Tr0_nidhits;
		short Tr1_nidhits;
		short Tr0_ntrhits;
		short Tr1_ntrhits;
		short Tr0_lastgap;
		short Tr1_lastgap;
		float Tr0_DG0;
		float Tr1_DG0;
		float Tr0_DDG0;
		float Tr1_DDG0;
		float Tr0_trchi2;
		float Tr1_trchi2;

		float dca_r;
		float dca_z;
		float Tr0_chi2_fvtx;
		float Tr1_chi2_fvtx;
		float Tr0_chi2_fvtxmutr;
		float Tr1_chi2_fvtxmutr;

		float Tr0_x0;
		float Tr1_x0;
		float Tr0_y0;
		float Tr1_y0;

		float Tr0_px;
		float Tr1_px;
		float Tr0_py;
		float Tr1_py;
		float Tr0_pz;
		float Tr1_pz;

		float Tr0_x0_fvtxmutr ;
		float Tr0_y0_fvtxmutr ;
		float Tr0_z0_fvtxmutr ;
		float Tr0_px_fvtxmutr ;
		float Tr0_py_fvtxmutr ;
		float Tr0_pz_fvtxmutr ;
		float Tr0_x0_fvtx;
		float Tr0_y0_fvtx;
		float Tr0_z0_fvtx;
		float Tr0_px_fvtx;
		float Tr0_py_fvtx;
		float Tr0_pz_fvtx;

		float Tr1_x0_fvtxmutr ;
		float Tr1_y0_fvtxmutr ;
		float Tr1_z0_fvtxmutr ;
		float Tr1_px_fvtxmutr ;
		float Tr1_py_fvtxmutr ;
		float Tr1_pz_fvtxmutr ;
		float Tr1_x0_fvtx;
		float Tr1_y0_fvtx;
		float Tr1_z0_fvtx;
		float Tr1_px_fvtx;
		float Tr1_py_fvtx;
		float Tr1_pz_fvtx;

		float Tr0_xst1;
		float Tr1_xst1;
		float Tr0_yst1;
		float Tr1_yst1;
		float Tr0_dca_r;
		float Tr1_dca_r;
		float Tr0_dca_z;
		float Tr1_dca_z;
		float Tr0_dr_fvtx;
		float Tr1_dr_fvtx;
		float Tr0_dphi_fvtx;
		float Tr1_dphi_fvtx;
		float Tr0_dtheta_fvtx;
		float Tr1_dtheta_fvtx;
		float Evt_fvtxX;
		float Evt_fvtxY;
		float Evt_fvtxZ;
		float Evt_bbcZ;
	        bool from_same_vertex;
		float Mutr_tr0_closeZ0;
		float Mutr_tr1_closeZ0;
		int Mutr_tr0_charge;
		int Mutr_tr1_charge;
		unsigned int Tr0_sngl_index;
		unsigned int Tr1_sngl_index;
		int nTracklets;
		int nTracklets_all;
		int nTracklets_samearm_all;
		int nTracklets_samearm;
		int bbindex;
		int ccindex;
		int flavor_index;

		std::string _out_name;
		DiMuonContainer *dimu_container;
		SingleMuonContainer *singlemuoncontainer;
		FvtxSnglPrimVtx * fvtx;
		DSTReader *fvtx_trk_map;
		VtxOut *vtxout;
		PHPythiaContainer *phpythia;
		TFile *_out_file;
 		TH1 *Pz_dist;
  		char namest[64];




};

#endif
