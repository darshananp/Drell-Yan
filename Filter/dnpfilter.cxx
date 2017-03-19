#include <PHGeometry.h>
#include <PHGlobal.h>
#include <Fun4AllReturnCodes.h>
#include <Fun4AllServer.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <MUTOO.h>
#include <RunHeader.h>
#include <TriggerHelper.h>
#include <TrigRunLvl1.h>
#include <TrigLvl1.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <PHPythiaHeader.h>
#include <PHPythiaContainerV2.h>
#include <VtxOut.h>
#include <TFvtxCompactTrk.h>
#include <TFvtxCompactTrkMap.h>
#include <FvtxPrimVtxContainer.h>
#include <FvtxPrimVtxContainer_v2.h>
#include "FvtxPrimVtxContainer_v1.h"
#include "FvtxSnglPrimVtx.h"
#include "FvtxSnglPrimVtx_v1.h"
#include "FvtxSnglPrimVtx_v2.h"
#include <SingleMuonContainer.h>
#include <DiMuonContainer.h>
//#include <DiMuonContainer_v9.h>
#include <SingleMuon.h>
#include <DiMuon.h>
//#include <DiMuon_v7.h>
//#include <DiMuon_v5.h>
#include "dnpfilter.h"

using namespace std;
//typedef PHIODataNode<PHObject> PHObjectNode_t;
//! global initialization - 1st step in Fun 4 all, called once at startup,Introduce histos and trees
int dnpfilter::Init(PHCompositeNode *topNode)
{
	init_tree();
	init_histos();

	nevents=0;

	return 0;
}

int dnpfilter::process_event(PHCompositeNode *topNode)
{
	if(nevents%100000==0){
		cout << "dnpfilter::process_event Event: " << nevents << endl;
	}
	nevents++;
	int runnumber = 0;
	RunHeader* runheader = findNode::getClass<RunHeader>(topNode, "RunHeader" );
	if(!runheader && !_is_sim){
		cout<<"WARN: "<<PHWHERE<<" RunHeader is not found"<<endl;
	}
	if(runheader) runnumber = runheader->get_RunNumber();

//---------------------------------------------------------------------------------------------------------

	singlemuoncontainer = findNode::getClass<SingleMuonContainer>(topNode,"SingleMuonContainer");
	if (!singlemuoncontainer){
		cout << "dnpfilter:: SingleMuonContainer not in Node Tree" << endl;
		return ABORTRUN;
	}
//---------------------------------------------------------------------------------------------------------
	
	dimu_container = findNode::getClass<DiMuonContainer>(topNode,"DiMuonContainer");
	if (!dimu_container){
		cout << "dnpfilter:: DiMuonContainer not in Node Tree" << endl;
		return ABORTRUN;
	}

//--------------------------------------------------------------------------------------------------------
	fvtx_trk_map = findNode::getClass<DSTReader>(topNode,"DSTReader");
	if(!fvtx_trk_map){
		cout << PHWHERE << "dnpfilter:: DSTReader not in Node Tree" << endl;
		return ABORTRUN;
	}
//--------------------------------------------------------------------------------------------------------
	vtxout = findNode::getClass<VtxOut>(topNode, "VtxOut" );
	if (!vtxout)
	{
		cout << PHWHERE << "dnpfilter:: VtxOut not in Node Tree" << endl;
		return ABORTRUN;
	}
//--------------------------------------------------------------------------------------------------------
	if (_is_sim && _is_mb){
        phpythia = findNode::getClass<PHPythiaContainer>(topNode,"PHPythia");

	if (!phpythia){
		cout << PHWHERE << "dnpfilter:: phpythia not in Node Tree" << endl;
    		return ABORTEVENT;
	}
	}
//--------------------------------------------------------------------------------------------------------
	TriggerHelper *trhelper = NULL; 
	TrigRunLvl1 *trigrunlvl1 = NULL;
	TrigLvl1 *triglvl1 = NULL;

	if(!_is_sim)
	{
		trhelper = new TriggerHelper(topNode);

		trigrunlvl1 = trhelper->get_trigRunLvl1();
		if(!trigrunlvl1){
			cout << PHWHERE << "WARN: Can not find TrigRunLvl1"<<endl;
			return ABORTRUN;
		}

		triglvl1 = trhelper->get_trigLvl1();
		if(!triglvl1){
			cout << PHWHERE << "WARN: Can not find TrigLvl1"<<endl;
			return ABORTRUN;
		}
		lvl1_trigscaled = triglvl1->get_lvl1_trigscaled();
	}

	lvl1scaledown = 0;
	bool is_good_event = true;
	if(!_is_sim){
		std::string trigger_name = "((MUIDLL1_N2D||S2D)||(N1D&S1D))&BBCLL1(noVtx)";
		//! check whether this trigger name is enabled
		if(!trhelper->getLevel1BitMask(trigger_name)) is_good_event = false;
		//! check whether this trigger bit is fired
		unsigned int bit_mask_position = trhelper->getLevel1BitNumber(trigger_name);
		if(bit_mask_position>31 || bit_mask_position<0) is_good_event = false;
		if(!triglvl1->get_lvl1_trigscaled_bit(bit_mask_position)) is_good_event = false;
		lvl1scaledown = trhelper->getLevel1Scaledown(trigger_name.data());

	}

	if(!is_good_event){
		delete trhelper;
		trhelper = NULL;
		trigrunlvl1 = NULL;
		return DISCARDEVENT;
	}


//--------------------------------------------------------------------------------------------------------
	for(unsigned int idimu=0;idimu<dimu_container->get_nDiMuons();idimu++){
		DiMuon *dimuon = dynamic_cast<DiMuon*> (dimu_container->get_DiMuon(idimu));
		if(verbosity>1) dimuon->identify(cout);
		if(!dimuon) continue;


		if(!set_dimuon_vars(dimuon)) continue;

	        nTracklets = GetNumberofTracklets(fvtx_trk_map,dimuon);
	        nTracklets_all = GetNumberofTracklets_all(fvtx_trk_map,dimuon);
		nTracklets_samearm = GetNumberofTracklets_samearm(fvtx_trk_map,dimuon);
		nTracklets_samearm_all = GetNumberofTracklets_samearm_all(fvtx_trk_map,dimuon);
//----------------------------------------------start of cut---------------------------------------------

                if(!pass_dimuon_cut(dimuon)) continue;

		short  Tr0_vtx_index = dimuon->get_Tr0_vtx_index();
		short  Tr1_vtx_index = dimuon->get_Tr1_vtx_index();
		int index_diff = Tr0_vtx_index - Tr1_vtx_index;
		if(index_diff==0) {from_same_vertex = true;}
		if(index_diff!=0) {from_same_vertex = false;}

//----------------------------------------------Single Muon Calculations--------------------------------
                if(momentum_matching(dimuon)) {
		Mutr_tr0_closeZ0 = set_singlemuon_vars(Tr0_sngl_index);
		Mutr_tr1_closeZ0 = set_singlemuon_vars(Tr1_sngl_index);
		Mutr_tr0_charge  = set_singlemuon_charge(Tr0_sngl_index);
		Mutr_tr1_charge  = set_singlemuon_charge(Tr1_sngl_index);
		Tr0_x0 		 = set_singlemuon_x0(Tr0_sngl_index);
		Tr1_x0 		 = set_singlemuon_x0(Tr1_sngl_index);
		Tr0_y0 		 = set_singlemuon_y0(Tr0_sngl_index);
		Tr1_y0 		 = set_singlemuon_y0(Tr1_sngl_index);
		Tr0_pz_fvtxmutr  = set_singlemuon_fvtxmutrz(Tr0_sngl_index);
		Tr1_pz_fvtxmutr  = set_singlemuon_fvtxmutrz(Tr1_sngl_index);
		Tr0_px_fvtxmutr  = set_singlemuon_fvtxmutrx(Tr0_sngl_index);
		Tr1_px_fvtxmutr  = set_singlemuon_fvtxmutrx(Tr1_sngl_index);
		Tr0_py_fvtxmutr  = set_singlemuon_fvtxmutry(Tr0_sngl_index);
		Tr1_py_fvtxmutr  = set_singlemuon_fvtxmutry(Tr1_sngl_index);
		}
                if(!momentum_matching(dimuon)) {
		Mutr_tr0_closeZ0 = -9999;
		Mutr_tr1_closeZ0 = -9999;
		Mutr_tr0_charge  = -9999;
		Mutr_tr1_charge  = -9999;
		Tr0_pz_fvtxmutr  = -9999;
		Tr1_pz_fvtxmutr  = -9999;
		}

		if (_is_sim && _is_mb){
		bbindex = set_processtype(phpythia);
		ccindex = set_processtypetwo(phpythia);
		if (bbindex>0){flavor_index = 1;}
		if (bbindex<1 && ccindex>0){flavor_index = 2;}
		if (bbindex<1 && ccindex<1){flavor_index = 3;}
		}
		//if(!fill_dimuon_hist(dimuon)) continue;	        
                T->Fill();


	}		
		return EVENT_OK;
}
//-----------------------------------------Set values for Dimuon Variables------------------------------
bool dnpfilter::set_dimuon_vars(const DiMuon *dimuon)
{
	Evt_vtxchi2 = dimuon->get_Evt_vtxchi2();
	Evt_vtxoor = dimuon->get_Evt_vtxoor();
	mass = dimuon->get_mass();
	mass_fvtx = dimuon->get_mass_fvtx();
	mass_fvtxmutr = dimuon->get_mass_fvtxmutr();
	pT = dimuon->get_pT();
	pz = dimuon->get_pz();
	rapidity = dimuon->get_rapidity();
	charge = dimuon->get_charge();

	same_event = dimuon->get_same_event();
	dca_r = dimuon->get_dca_r();
	dca_z = dimuon->get_dca_z();

	Tr0_idhits = dimuon->get_Tr0_idhits();
	Tr1_idhits = dimuon->get_Tr1_idhits();
	Tr0_nidhits = dimuon->get_Tr0_nidhits();
	Tr1_nidhits = dimuon->get_Tr1_nidhits();
	Tr0_ntrhits = dimuon->get_Tr0_ntrhits();
	Tr1_ntrhits = dimuon->get_Tr1_ntrhits();
	Tr0_lastgap = dimuon->get_Tr0_lastgap();
	Tr1_lastgap = dimuon->get_Tr1_lastgap();
	Tr0_DG0 = dimuon->get_Tr0_DG0();
	Tr1_DG0 = dimuon->get_Tr1_DG0();
	Tr0_DDG0 = dimuon->get_Tr0_DDG0();
	Tr1_DDG0 = dimuon->get_Tr1_DDG0();
	Tr0_trchi2 = dimuon->get_Tr0_trchi2();
	Tr1_trchi2 = dimuon->get_Tr1_trchi2();

	Tr0_px = dimuon->get_Tr0_px();
	Tr1_py = dimuon->get_Tr1_py();
	Tr0_pz = dimuon->get_Tr0_pz();
	Tr1_px = dimuon->get_Tr1_px();
	Tr1_py = dimuon->get_Tr1_py();
	Tr1_pz = dimuon->get_Tr1_pz();


//-----------------------------------------------fvtx px,py,pz,x,y,z------------------------------------
	Tr0_x0_fvtxmutr   = dimuon->get_Tr0_x0_fvtxmutr();
	Tr0_y0_fvtxmutr   = dimuon->get_Tr0_y0_fvtxmutr();
	Tr0_z0_fvtxmutr   = dimuon->get_Tr0_z0_fvtxmutr();

	Tr0_x0_fvtx	  = dimuon->get_Tr0_x0_fvtx();
	Tr0_y0_fvtx       = dimuon->get_Tr0_y0_fvtx();
	Tr0_z0_fvtx       = dimuon->get_Tr0_z0_fvtx();
	Tr0_px_fvtx       = dimuon->get_Tr0_px_fvtx();
	Tr0_py_fvtx       = dimuon->get_Tr0_py_fvtx();
	Tr0_pz_fvtx       = dimuon->get_Tr0_pz_fvtx();

	Tr1_x0_fvtxmutr   = dimuon->get_Tr1_x0_fvtxmutr();
	Tr1_y0_fvtxmutr   = dimuon->get_Tr1_y0_fvtxmutr();
	Tr1_z0_fvtxmutr   = dimuon->get_Tr1_z0_fvtxmutr();

	Tr1_x0_fvtx	  = dimuon->get_Tr1_x0_fvtx();
	Tr1_y0_fvtx       = dimuon->get_Tr1_y0_fvtx();
	Tr1_z0_fvtx       = dimuon->get_Tr1_z0_fvtx();
	Tr1_px_fvtx       = dimuon->get_Tr1_px_fvtx();
	Tr1_py_fvtx       = dimuon->get_Tr1_py_fvtx();
	Tr1_pz_fvtx       = dimuon->get_Tr1_pz_fvtx();
//--------------------------------------------Matching Variables----------------------------------------
	Tr0_chi2_fvtx = dimuon->get_Tr0_chi2_fvtx();
	Tr1_chi2_fvtx = dimuon->get_Tr1_chi2_fvtx();
	Tr0_chi2_fvtxmutr = dimuon->get_Tr0_chi2_fvtxmutr();
	Tr1_chi2_fvtxmutr = dimuon->get_Tr1_chi2_fvtxmutr();
	Tr0_xst1 = dimuon->get_Tr0_xst1();
	Tr1_xst1 = dimuon->get_Tr1_xst1();
	Tr0_yst1 = dimuon->get_Tr0_yst1();
	Tr1_yst1 = dimuon->get_Tr1_yst1();
	Tr0_dca_r = dimuon->get_Tr0_dca_r();
	Tr1_dca_r = dimuon->get_Tr1_dca_r();
	Tr0_dca_z = dimuon->get_Tr0_dca_z();
	Tr1_dca_z = dimuon->get_Tr1_dca_z();
	Tr0_dr_fvtx = dimuon->get_Tr0_dr_fvtx();
	Tr1_dr_fvtx = dimuon->get_Tr1_dr_fvtx();
	Tr0_dphi_fvtx = dimuon->get_Tr0_dphi_fvtx();
	Tr1_dphi_fvtx = dimuon->get_Tr1_dphi_fvtx();
	Tr0_dtheta_fvtx = dimuon->get_Tr0_dtheta_fvtx();
	Tr1_dtheta_fvtx = dimuon->get_Tr1_dtheta_fvtx();

//------------------------------------------- vertex information ---------------------------------------
	Evt_fvtxX = dimu_container->get_Evt_fvtxX();
	Evt_fvtxY = dimu_container->get_Evt_fvtxY();
	Evt_fvtxZ = dimu_container->get_Evt_fvtxZ();
	Evt_bbcZ  = dimu_container->get_Evt_bbcZ();

	return true;
}

bool dnpfilter::momentum_matching(const DiMuon *dimuon)
{	Tr0_sngl_index = 9999;
	Tr1_sngl_index = 9999;
	bool is_matching = false;
        float momentum_matching_cut = 0.0001;
	unsigned int Tr0_sngl_index_best = 9999;
	float Tr0_diff_this = 9999;
	float Tr0_diff_best = 9999;
	unsigned int Tr1_sngl_index_best = 9999;
	float Tr1_diff_this = 9999;
	float Tr1_diff_best = 9999;

	for(unsigned int isngl=0;isngl<singlemuoncontainer->get_nSingleMuons();isngl++)
	{
		SingleMuon *singlemuon = singlemuoncontainer->get_SingleMuon(isngl);
		if(!singlemuon){
			continue;
		}

		Tr0_diff_this = 
			sqrt(
					pow((dimuon->get_Tr0_px() - singlemuon->get_px()),2) +
					pow((dimuon->get_Tr0_py() - singlemuon->get_py()),2) +
					pow((dimuon->get_Tr0_pz() - singlemuon->get_pz()),2)
					);

		Tr1_diff_this = 
			sqrt(
					pow((dimuon->get_Tr1_px() - singlemuon->get_px()),2) +
					pow((dimuon->get_Tr1_py() - singlemuon->get_py()),2) +
					pow((dimuon->get_Tr1_pz() - singlemuon->get_pz()),2)
					);

		if(Tr0_diff_this<Tr0_diff_best && Tr0_diff_this<Tr1_diff_this){
			Tr0_diff_best = Tr0_diff_this;
			Tr0_sngl_index_best = isngl;
		}

		if(Tr1_diff_this<Tr1_diff_best && Tr1_diff_this<Tr0_diff_this){
			Tr1_diff_best = Tr1_diff_this;
			Tr1_sngl_index_best = isngl;
		}
	}

	if(Tr0_diff_best < momentum_matching_cut &&
			Tr1_diff_best < momentum_matching_cut){
		Tr0_sngl_index = Tr0_sngl_index_best;
		Tr1_sngl_index = Tr1_sngl_index_best;
		is_matching = true;
	}


	return is_matching;
}


float dnpfilter::set_singlemuon_vars(unsigned int index) const
{
	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(index);
	if(!muon){
		cout<<"WARN: "<<PHWHERE<<"!muon"<<endl;
		return -9999;
	}

	float dca_vertex = (muon->get_z0()-(((muon->get_x0()-0.1484)*muon->get_px()+(muon->get_y0()-0.04226)*muon->get_py())*muon->get_pz())/(muon->get_px()*muon->get_px()+muon->get_py()*muon->get_py()))-dimu_container->get_Evt_fvtxZ();

	//float dca_vertex = ((((muon->get_x0()-0.1484)*muon->get_px()+(muon->get_y0()-0.04226)*muon->get_py())*muon->get_pz())/(muon->get_px()*muon->get_px()+muon->get_py()*muon->get_py()))-dimu_container->get_Evt_fvtxZ();

	return dca_vertex;
}

int dnpfilter::set_singlemuon_charge(unsigned int index) const
{
	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(index);
	if(!muon){
		cout<<"WARN: "<<PHWHERE<<"!muon"<<endl;
		return -9999;
	}
	int sigl_charge = muon->get_charge();
	return sigl_charge;
}

//***************
float dnpfilter::set_singlemuon_x0(unsigned int index) const
{
	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(index);
	if(!muon){
		cout<<"WARN: "<<PHWHERE<<"!muon"<<endl;
		return -9999;
	}
	double mutr_x = muon->get_x0();
	return mutr_x;
}

float dnpfilter::set_singlemuon_y0(unsigned int index) const
{
	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(index);
	if(!muon){
		cout<<"WARN: "<<PHWHERE<<"!muon"<<endl;
		return -9999;
	}
	double mutr_y = muon->get_y0();
	return mutr_y;
}

//*********

float dnpfilter::set_singlemuon_fvtxmutrz(unsigned int index) const
{
	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(index);
	if(!muon){
		cout<<"WARN: "<<PHWHERE<<"!muon"<<endl;
		return -9999;
	}
	double fvtxmutr_z = muon->get_pz_fvtxmutr();
	return fvtxmutr_z;
}

float dnpfilter::set_singlemuon_fvtxmutrx(unsigned int index) const
{
	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(index);
	if(!muon){
		cout<<"WARN: "<<PHWHERE<<"!muon"<<endl;
		return -9999;
	}
	double fvtxmutr_x = muon->get_px_fvtxmutr();
	return fvtxmutr_x;
}

float dnpfilter::set_singlemuon_fvtxmutry(unsigned int index) const
{
	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(index);
	if(!muon){
		cout<<"WARN: "<<PHWHERE<<"!muon"<<endl;
		return -9999;
	}
	double fvtxmutr_y = muon->get_py_fvtxmutr();
	return fvtxmutr_y;
}


int dnpfilter::set_processtype(PHPythiaContainer *phpythia)
{
	int bbreturn=0;
        
	int bbarray[53]={511,521,10511,10521,513,523,10513,10523,20513,20523,515,525,531,10531,533,10533,
	20533,535,541,10541,543,10543,20543,545,551,10551,100551,110551,200551,210551 ,553,10553,20553,
	30553,100553,110553,120553,130553,200553,210553,220553,300553,9000553,9010553,555,10555,20555,
	100555,110555,120555,200555,557,100557};
	int nbbevents = 0;
	int npart = phpythia->size();
        for (int i = 0; i < npart; i++) {
        TMCParticle *part = phpythia->getParticle(i);

	for (int j=0;j<53;j++){
        if ( abs(part->GetKF())!=bbarray[j] ) continue;
	cout<<"this is bb event = "<<part->GetKF()<< "  in parton ="<< i <<endl;
	nbbevents++;
	}

	}
	if (nbbevents>0){bbreturn = 1;}
  	return bbreturn;
}


int dnpfilter::set_processtypetwo(PHPythiaContainer *phpythia)
{
	int ccreturn=0;
        
	int ccarray[31]={411,421,10411,10421,413,423,10413,10423,20413,20423,415,425,431,10431,433,
	10433,20433,435,441,10441,100441,443,10443,20443,100443,30443,9000443,9010443,9020443,445,9000445};
	int nccevents = 0;
	int npart = phpythia->size();
        for (int i = 0; i < npart; i++) {
        TMCParticle *part = phpythia->getParticle(i);

	for (int j=0;j<31;j++){
        if ( abs(part->GetKF())!=ccarray[j] ) continue;
	cout<<"this is cc event = "<<part->GetKF()<< "  in parton ="<< i <<endl;
	nccevents++;
	}

	}
	if (nccevents>0){ccreturn = 1;}
  	return ccreturn;
}



bool dnpfilter::fill_dimuon_hist(const DiMuon *dimuon)
{
  if(fabs(pz) <= 150.0){
	Pz_dist->Fill(pz);
    }
	return true;
}



dnpfilter::dnpfilter(const std::string &out_name):_out_name(out_name)

{  	_use_cut_tracklet_chi2=false;
	_cut_tracklet_chi2=4.;
	_cut_tracklet_dcar=1.5;
	_use_2_hit_tracklet=false;
	_is_sim = false;
	_is_mb = false;
	_out_file = new TFile((_out_name+"_"+".root").data(),"RECREATE");
}

//------------------------------------------------Need Work---------------------------------------------

int dnpfilter::GetNumberofTracklets(DSTReader *fvtx_trk_map, const DiMuon *dimuon)
{
	if(!fvtx_trk_map){
	cout<<"EXCEPTION: "<<PHWHERE<<endl;
	return NULL;
	}
	
	int ntrklets = 0;
	TClonesArray *array = fvtx_trk_map->get_FvtxCompactTrk();
	for (int i = 0; i < array->GetSize(); i++) {
	TFvtxCompactTrk *tracklet = dynamic_cast<TFvtxCompactTrk*> (array->At(i));
	
	if(!tracklet){
	cout<<"No tracklet"<<__LINE__<<"size: "<<array->GetSize()<<endl;
	break;
	}
	
	if(_use_cut_tracklet_chi2 && (tracklet->get_chi2_ndf() > _cut_tracklet_chi2)) continue;
	if(!_use_2_hit_tracklet && tracklet->get_nhits() <= 2) continue;
	/*float phi1 = tracklet->get_fvtx_phi();
	float theta1 = tracklet->get_fvtx_theta();
	PHVector vector1(cos(phi1)*sin(theta1),sin(phi1)*sin(theta1),cos(theta1));
	PHPoint point1 = tracklet->get_fvtx_vtx();
	PHLine track = PHLine(point1, vector1);
	float px = track.getDirection().getX();
	float py = track.getDirection().getY();
	float pt = sqrt(px*px+py*py);
*/
	SingleMuon *muon0 = singlemuoncontainer->get_SingleMuon(0);
	SingleMuon *muon1 = singlemuoncontainer->get_SingleMuon(1);

	//if (dimuon->get_mass_fvtxmutr()!=0){
	//cout << muon0->get_pz_fvtxmutr()<<"    "<<muon1->get_pz_fvtxmutr()<<endl;}

	float xx0 = tracklet->get_fvtx_vtx().getX()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*cos(tracklet->get_fvtx_phi());
	float yy0 = tracklet->get_fvtx_vtx().getY()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*sin(tracklet->get_fvtx_phi());
	float x0y0=sqrt((xx0 - dimu_container->get_Evt_fvtxX())*(xx0 - dimu_container->get_Evt_fvtxX()) +
	(yy0 - dimu_container->get_Evt_fvtxY())*(yy0 - dimu_container->get_Evt_fvtxY()));

	if (abs(TMath::ATan2(sqrt(muon0->get_px_fvtxmutr()*muon0->get_px_fvtxmutr()+
	muon0->get_py_fvtxmutr()*muon0->get_py_fvtxmutr()),muon0->get_pz_fvtxmutr())-
	tracklet->get_fvtx_theta()) >0.001 && abs(TMath::ATan2(sqrt(muon1->get_px_fvtxmutr()*
	muon1->get_px_fvtxmutr()+muon1->get_py_fvtxmutr()*muon1->get_py_fvtxmutr()),
	muon1->get_pz_fvtxmutr())-tracklet->get_fvtx_theta()) >0.001 && abs(TMath::ATan2(muon0->get_py_fvtxmutr(),
	muon0->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && abs(TMath::ATan2(muon1->get_py_fvtxmutr(),
	muon1->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && tracklet->get_fvtx_theta()+0!=0 && 
	x0y0 < _cut_tracklet_dcar){
	ntrklets++;
	}

	}

	return ntrklets;
        
}


int dnpfilter::GetNumberofTracklets_all(DSTReader *fvtx_trk_map, const DiMuon *dimuon)
{
	if(!fvtx_trk_map){
	cout<<"EXCEPTION: "<<PHWHERE<<endl;
	return NULL;
	}
	
	int ntrklets = 0;
	TClonesArray *array = fvtx_trk_map->get_FvtxCompactTrk();
	for (int i = 0; i < array->GetSize(); i++) {
	TFvtxCompactTrk *tracklet = dynamic_cast<TFvtxCompactTrk*> (array->At(i));
	
	if(!tracklet){
	cout<<"No tracklet"<<__LINE__<<"size: "<<array->GetSize()<<endl;
	break;
	}
	
	if(_use_cut_tracklet_chi2 && (tracklet->get_chi2_ndf() > _cut_tracklet_chi2)) continue;
	if(!_use_2_hit_tracklet && tracklet->get_nhits() <= 2) continue;

	//cout<<" x=   "<<tracklet->get_fvtx_vtx().getX()<<" y=   "<<tracklet->get_fvtx_vtx().getZ()<<"  z=  "
	//<<tracklet->get_fvtx_vtx().getY()<<endl;

	//cout<<" theta=   "<<tracklet->get_fvtx_theta()<<" phi=   "<<tracklet->get_fvtx_phi()<<"  zevt=  "
	//<<dimu_container->get_Evt_fvtxZ()<<endl;

	float xx0 = tracklet->get_fvtx_vtx().getX()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*cos(tracklet->get_fvtx_phi());
	float yy0 = tracklet->get_fvtx_vtx().getY()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*sin(tracklet->get_fvtx_phi());
	float x0y0=sqrt((xx0 - dimu_container->get_Evt_fvtxX())*(xx0 - dimu_container->get_Evt_fvtxX()) +
	(yy0 - dimu_container->get_Evt_fvtxY())*(yy0 - dimu_container->get_Evt_fvtxY()));

	PHPoint vertex;
	vertex = vtxout->get_Vertex("FVTX");
	//cout<<	vertex.getX()<<"  z=  "<<	vertex.getZ()<<"  zvtz=  "<<dimu_container->get_Evt_fvtxZ()<<endl;
	//cout<<"  xx0=  "<<xx0<<"  yy0=  "<<yy0<<"  x0y0=  "<<x0y0<<endl;
        //cout<<" before z ="<<dimuon->get_Tr1_pz_fvtxmutr()<<endl;
	if (tracklet->get_fvtx_theta()+0!=0 && x0y0 < _cut_tracklet_dcar){
	ntrklets++;
	}

	}

	return ntrklets;
        
}


int dnpfilter::GetNumberofTracklets_samearm_all(DSTReader *fvtx_trk_map, const DiMuon *dimuon)
{
	if(!fvtx_trk_map){
	cout<<"EXCEPTION: "<<PHWHERE<<endl;
	return NULL;
	}

	SingleMuon *muon = singlemuoncontainer->get_SingleMuon(0);

	int ntrklets = 0;
	TClonesArray *array = fvtx_trk_map->get_FvtxCompactTrk();
	for (int i = 0; i < array->GetSize(); i++) {
	TFvtxCompactTrk *tracklet = dynamic_cast<TFvtxCompactTrk*> (array->At(i));
	
	if(!tracklet){
	cout<<"No tracklet"<<__LINE__<<"size: "<<array->GetSize()<<endl;
	break;
	}
	
	if(_use_cut_tracklet_chi2 && (tracklet->get_chi2_ndf() > _cut_tracklet_chi2)) continue;
	if(!_use_2_hit_tracklet && tracklet->get_nhits() <= 2) continue;



	double track_eta= TMath::ATanH(muon->get_pz_fvtxmutr()/(sqrt(muon->get_px_fvtxmutr()*muon->get_px_fvtxmutr()+
	muon->get_py_fvtxmutr()*muon->get_py_fvtxmutr()+muon->get_pz_fvtxmutr()*muon->get_pz_fvtxmutr())));

	double tracklet_eta=0-log(tan((tracklet->get_fvtx_theta())/2));

	if(track_eta*tracklet_eta<0) continue;

	float xx0 = tracklet->get_fvtx_vtx().getX()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*cos(tracklet->get_fvtx_phi());
	float yy0 = tracklet->get_fvtx_vtx().getY()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*sin(tracklet->get_fvtx_phi());
	float x0y0=sqrt((xx0 - dimu_container->get_Evt_fvtxX())*(xx0 - dimu_container->get_Evt_fvtxX()) +
	(yy0 - dimu_container->get_Evt_fvtxY())*(yy0 - dimu_container->get_Evt_fvtxY()));

	//PHPoint vertex;
	//vertex = vtxout->get_Vertex("FVTX");

	if (tracklet->get_fvtx_theta()+0!=0 && x0y0 < _cut_tracklet_dcar){
	ntrklets++;
	}

	}

	return ntrklets;
        
}

int dnpfilter::GetNumberofTracklets_samearm(DSTReader *fvtx_trk_map, const DiMuon *dimuon)
{
	if(!fvtx_trk_map){
	cout<<"EXCEPTION: "<<PHWHERE<<endl;
	return NULL;
	}

	SingleMuon *muon0 = singlemuoncontainer->get_SingleMuon(0);
	SingleMuon *muon1 = singlemuoncontainer->get_SingleMuon(1);

	int ntrklets = 0;
	TClonesArray *array = fvtx_trk_map->get_FvtxCompactTrk();
	for (int i = 0; i < array->GetSize(); i++) {
	TFvtxCompactTrk *tracklet = dynamic_cast<TFvtxCompactTrk*> (array->At(i));
	
	if(!tracklet){
	cout<<"No tracklet"<<__LINE__<<"size: "<<array->GetSize()<<endl;
	break;
	}
	
	if(_use_cut_tracklet_chi2 && (tracklet->get_chi2_ndf() > _cut_tracklet_chi2)) continue;
	if(!_use_2_hit_tracklet && tracklet->get_nhits() <= 2) continue;



	double track_eta= TMath::ATanH(muon0->get_pz_fvtxmutr()/(sqrt(muon0->get_px_fvtxmutr()*muon0->get_px_fvtxmutr()+
	muon0->get_py_fvtxmutr()*muon0->get_py_fvtxmutr()+muon0->get_pz_fvtxmutr()*muon0->get_pz_fvtxmutr())));

	double tracklet_eta=0-log(tan((tracklet->get_fvtx_theta())/2));

	if(track_eta*tracklet_eta<0) continue;

	float xx0 = tracklet->get_fvtx_vtx().getX()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*cos(tracklet->get_fvtx_phi());
	float yy0 = tracklet->get_fvtx_vtx().getY()-(tracklet->get_fvtx_vtx().getZ()- dimu_container->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*sin(tracklet->get_fvtx_phi());
	float x0y0=sqrt((xx0 - dimu_container->get_Evt_fvtxX())*(xx0 - dimu_container->get_Evt_fvtxX()) +
	(yy0 - dimu_container->get_Evt_fvtxY())*(yy0 - dimu_container->get_Evt_fvtxY()));

	//PHPoint vertex;
	//vertex = vtxout->get_Vertex("FVTX");

	if (abs(TMath::ATan2(sqrt(muon0->get_px_fvtxmutr()*muon0->get_px_fvtxmutr()+
	muon0->get_py_fvtxmutr()*muon0->get_py_fvtxmutr()),muon0->get_pz_fvtxmutr())-
	tracklet->get_fvtx_theta()) >0.001 && abs(TMath::ATan2(sqrt(muon1->get_px_fvtxmutr()*
	muon1->get_px_fvtxmutr()+muon1->get_py_fvtxmutr()*muon1->get_py_fvtxmutr()),
	muon1->get_pz_fvtxmutr())-tracklet->get_fvtx_theta()) >0.001 && abs(TMath::ATan2(muon0->get_py_fvtxmutr(),
	muon0->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && abs(TMath::ATan2(muon1->get_py_fvtxmutr(),
	muon1->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && tracklet->get_fvtx_theta()+0!=0 && 
	x0y0 < _cut_tracklet_dcar){
	ntrklets++;
	}

	}

	return ntrklets;
        
}

void dnpfilter::init_tree()
{
	T = new TTree("T","DiMu_Simple_TTree");

	Evt_vtxchi2 = -9999;
	Evt_vtxoor = -9999;
	mass = -9999;
	mass_fvtx = -9999;
	mass_fvtxmutr = -9999;
	pT = -9999;
	pz = -9999;
	rapidity = -9999;
	charge = -9999;
	lvl1_trigscaled = 0;
	lvl1scaledown = 0;

	same_event = false;
	dca_r = -9999;
	dca_z = -9999;

	Tr0_fvtx_dca_z = -9999;
	Tr1_fvtx_dca_z = -9999;
	Tr0_idhits = -9999;
	Tr1_idhits = -9999;
	Tr0_nidhits = -9999;
	Tr1_nidhits = -9999;
	Tr0_ntrhits = -9999;
	Tr1_ntrhits = -9999;
	Tr0_lastgap = -9999;
	Tr1_lastgap = -9999;
	Tr0_DG0 = -9999;
	Tr1_DG0 = -9999;
	Tr0_DDG0 = -9999;
	Tr1_DDG0 = -9999;
	Tr0_trchi2 = -9999;
	Tr1_trchi2 = -9999;

	Tr0_chi2_fvtx = -9999;
	Tr1_chi2_fvtx = -9999;
	Tr0_chi2_fvtxmutr = -9999;
	Tr1_chi2_fvtxmutr = -9999;

	Tr0_x0 = -9999;
	Tr1_x0 = -9999;
	Tr0_y0 = -9999;
	Tr1_y0 = -9999;

	Tr0_px = -9999;
	Tr1_px = -9999;
	Tr0_py = -9999;
	Tr1_py = -9999;
	Tr0_pz = -9999;
	Tr1_pz = -9999;

	Tr0_x0_fvtxmutr   = -9999;
	Tr0_y0_fvtxmutr   = -9999;
	Tr0_z0_fvtxmutr   = -9999;
	Tr0_px_fvtxmutr   = -9999;
	Tr0_py_fvtxmutr   = -9999;
	Tr0_pz_fvtxmutr   = -9999;
	Tr0_x0_fvtx	  = -9999;
	Tr0_y0_fvtx       = -9999;
	Tr0_z0_fvtx       = -9999;
	Tr0_px_fvtx       = -9999;
	Tr0_py_fvtx       = -9999;
	Tr0_pz_fvtx       = -9999;

	Tr1_x0_fvtxmutr   = -9999;
	Tr1_y0_fvtxmutr   = -9999;
	Tr1_z0_fvtxmutr   = -9999;
	Tr1_px_fvtxmutr   = -9999;
	Tr1_py_fvtxmutr   = -9999;
	Tr1_pz_fvtxmutr   = -9999;
	Tr1_x0_fvtx	  = -9999;
	Tr1_y0_fvtx       = -9999;
	Tr1_z0_fvtx       = -9999;
	Tr1_px_fvtx       = -9999;
	Tr1_py_fvtx       = -9999;
	Tr1_pz_fvtx       = -9999;

	Tr0_dr_fvtx = -9999;
	Tr1_dr_fvtx = -9999;
	Tr0_dphi_fvtx = -9999;
	Tr1_dphi_fvtx = -9999;
	Tr0_dtheta_fvtx = -9999;
	Tr1_dtheta_fvtx = -9999;
	Evt_fvtxX = -9999;
	Evt_fvtxY = -9999;
	Evt_fvtxZ = -9999;
        Evt_bbcZ = -9999;
	from_same_vertex = false;
        Mutr_tr0_closeZ0=-9999;
	Mutr_tr0_charge =-9999;
	Mutr_tr1_charge =-9999;
        Mutr_tr1_closeZ0=-9999;
	nTracklets = -9999;	
	nTracklets_all = -9999;
	nTracklets_samearm= -9999;
	nTracklets_samearm_all= -9999;
	bbindex = 0;
	ccindex = 0;
	flavor_index = 0;

	T->Branch("Evt_vtxchi2",&Evt_vtxchi2,"Evt_vtxchi2/F");
	T->Branch("Evt_vtxoor",&Evt_vtxoor,"Evt_vtxoor/F");
	T->Branch("mass",&mass,"mass/F");
	T->Branch("mass_fvtx",&mass_fvtx,"mass_fvtx/F");
	T->Branch("mass_fvtxmutr",&mass_fvtxmutr,"mass_fvtxmutr/F");
	T->Branch("pT",&pT,"pT/F");
	T->Branch("pz",&pz,"pz/F");
	T->Branch("rapidity",&rapidity,"rapidity/F");
	T->Branch("charge",&charge,"charge/I");

	T->Branch("lvl1_trigscaled",&lvl1_trigscaled,"lvl1_trigscaled/i");
	T->Branch("lvl1scaledown",&lvl1scaledown,"lvl1scaledown/i");

	T->Branch("same_event",&same_event,"same_event/O");
	T->Branch("Tr0_fvtx_dca_z",&Tr0_fvtx_dca_z,"Tr0_fvtx_dca_z/F");
	T->Branch("Tr1_fvtx_dca_z",&Tr1_fvtx_dca_z,"Tr1_fvtx_dca_z/F");
	T->Branch("Tr0_idhits",&Tr0_idhits,"Tr0_idhits/S");
	T->Branch("Tr1_idhits",&Tr1_idhits,"Tr1_idhits/S");
	T->Branch("Tr0_nidhits",&Tr0_nidhits,"Tr0_nidhits/S");
	T->Branch("Tr1_nidhits",&Tr1_nidhits,"Tr1_nidhits/S");
	T->Branch("Tr0_ntrhits",&Tr0_ntrhits,"Tr0_ntrhits/S");
	T->Branch("Tr1_ntrhits",&Tr1_ntrhits,"Tr1_ntrhits/S");
	T->Branch("Tr0_lastgap",&Tr0_lastgap,"Tr0_lastgap/S");
	T->Branch("Tr1_lastgap",&Tr1_lastgap,"Tr1_lastgap/S");
	T->Branch("Tr0_DG0",&Tr0_DG0,"Tr0_DG0/F");
	T->Branch("Tr1_DG0",&Tr1_DG0,"Tr1_DG0/F");
	T->Branch("Tr0_DDG0",&Tr0_DDG0,"Tr0_DDG0/F");
	T->Branch("Tr1_DDG0",&Tr1_DDG0,"Tr1_DDG0/F");
	T->Branch("Tr0_trchi2",&Tr0_trchi2,"Tr0_trchi2/F");
	T->Branch("Tr1_trchi2",&Tr1_trchi2,"Tr1_trchi2/F");

	T->Branch("dca_r",&dca_r,"dca_r/F");
	T->Branch("dca_z",&dca_z,"dca_z/F");
	T->Branch("Tr0_chi2_fvtx",&Tr0_chi2_fvtx,"Tr0_chi2_fvtx/F");
	T->Branch("Tr1_chi2_fvtx",&Tr1_chi2_fvtx,"Tr1_chi2_fvtx/F");
	T->Branch("Tr0_chi2_fvtxmutr",&Tr0_chi2_fvtxmutr,"Tr0_chi2_fvtxmutr/F");
	T->Branch("Tr1_chi2_fvtxmutr",&Tr1_chi2_fvtxmutr,"Tr1_chi2_fvtxmutr/F");

	T->Branch("Tr0_x0",&Tr0_x0,"Tr0_x0/F");
	T->Branch("Tr0_y0",&Tr0_y0,"Tr0_y0/F");
	T->Branch("Tr1_x0",&Tr1_x0,"Tr1_x0/F");
	T->Branch("Tr1_y0",&Tr1_y0,"Tr1_y0/F");

	T->Branch("Tr0_px",&Tr0_px,"Tr0_px/F");
	T->Branch("Tr1_px",&Tr1_px,"Tr1_px/F");
	T->Branch("Tr0_py",&Tr0_py,"Tr0_py/F");
	T->Branch("Tr1_py",&Tr1_py,"Tr1_py/F");
	T->Branch("Tr0_pz",&Tr0_pz,"Tr0_pz/F");
	T->Branch("Tr1_pz",&Tr1_pz,"Tr1_pz/F");

	T->Branch("Tr0_x0_fvtx",&Tr0_x0_fvtx,"Tr0_x0_fvtx/F");
	T->Branch("Tr0_y0_fvtx",&Tr0_y0_fvtx,"Tr0_y0_fvtx/F");
	T->Branch("Tr0_z0_fvtx",&Tr0_z0_fvtx,"Tr0_z0_fvtx/F");
	T->Branch("Tr0_px_fvtx",&Tr0_px_fvtx,"Tr0_px_fvtx/F");
	T->Branch("Tr0_py_fvtx",&Tr0_py_fvtx,"Tr0_py_fvtx/F");
	T->Branch("Tr0_pz_fvtx",&Tr0_pz_fvtx,"Tr0_pz_fvtx/F");

	T->Branch("Tr0_x0_fvtxmutr",&Tr0_x0_fvtxmutr,"Tr0_x0_fvtxmutr/F");
	T->Branch("Tr0_y0_fvtxmutr",&Tr0_y0_fvtxmutr,"Tr0_y0_fvtxmutr/F");
	T->Branch("Tr0_z0_fvtxmutr",&Tr0_z0_fvtxmutr,"Tr0_z0_fvtxmutr/F");
	T->Branch("Tr0_px_fvtxmutr",&Tr0_px_fvtxmutr,"Tr0_px_fvtxmutr/F");
	T->Branch("Tr0_py_fvtxmutr",&Tr0_py_fvtxmutr,"Tr0_py_fvtxmutr/F");
	T->Branch("Tr0_pz_fvtxmutr",&Tr0_pz_fvtxmutr,"Tr0_pz_fvtxmutr/F");

	T->Branch("Tr1_x0_fvtx",&Tr1_x0_fvtx,"Tr1_x0_fvtx/F");
	T->Branch("Tr1_y0_fvtx",&Tr1_y0_fvtx,"Tr1_y0_fvtx/F");
	T->Branch("Tr1_z0_fvtx",&Tr1_z0_fvtx,"Tr1_z0_fvtx/F");
	T->Branch("Tr1_px_fvtx",&Tr1_px_fvtx,"Tr1_px_fvtx/F");
	T->Branch("Tr1_py_fvtx",&Tr1_py_fvtx,"Tr1_py_fvtx/F");
	T->Branch("Tr1_pz_fvtx",&Tr1_pz_fvtx,"Tr1_pz_fvtx/F");

	T->Branch("Tr1_x0_fvtxmutr",&Tr1_x0_fvtxmutr,"Tr1_x0_fvtxmutr/F");
	T->Branch("Tr1_y0_fvtxmutr",&Tr1_y0_fvtxmutr,"Tr1_y0_fvtxmutr/F");
	T->Branch("Tr1_z0_fvtxmutr",&Tr1_z0_fvtxmutr,"Tr1_z0_fvtxmutr/F");
	T->Branch("Tr1_px_fvtxmutr",&Tr1_px_fvtxmutr,"Tr1_px_fvtxmutr/F");
	T->Branch("Tr1_py_fvtxmutr",&Tr1_py_fvtxmutr,"Tr1_py_fvtxmutr/F");
	T->Branch("Tr1_pz_fvtxmutr",&Tr1_pz_fvtxmutr,"Tr1_pz_fvtxmutr/F");

	T->Branch("Tr0_xst1",&Tr0_xst1,"Tr0_xst1/F");
	T->Branch("Tr1_xst1",&Tr1_xst1,"Tr1_xst1/F");
	T->Branch("Tr0_yst1",&Tr0_yst1,"Tr0_yst1/F");
	T->Branch("Tr1_yst1",&Tr1_yst1,"Tr1_yst1/F");
	T->Branch("Tr0_dca_r",&Tr0_dca_r,"Tr0_dca_r/F");
	T->Branch("Tr1_dca_r",&Tr1_dca_r,"Tr1_dca_r/F");
	T->Branch("Tr0_dca_z",&Tr0_dca_z,"Tr0_dca_z/F");
	T->Branch("Tr1_dca_z",&Tr1_dca_z,"Tr1_dca_z/F");
	T->Branch("Tr0_dr_fvtx",&Tr0_dr_fvtx,"Tr0_dr_fvtx/F");
	T->Branch("Tr1_dr_fvtx",&Tr1_dr_fvtx,"Tr1_dr_fvtx/F");
	T->Branch("Tr0_dphi_fvtx",&Tr0_dphi_fvtx,"Tr0_dphi_fvtx/F");
	T->Branch("Tr1_dphi_fvtx",&Tr1_dphi_fvtx,"Tr1_dphi_fvtx/F");
	T->Branch("Tr0_dtheta_fvtx",&Tr0_dtheta_fvtx,"Tr0_dtheta_fvtx/F");
	T->Branch("Tr1_dtheta_fvtx",&Tr1_dtheta_fvtx,"Tr1_dtheta_fvtx/F");

	T->Branch("Evt_fvtxX",&Evt_fvtxX,"Evt_fvtxX/F");
	T->Branch("Evt_fvtxY",&Evt_fvtxY,"Evt_fvtxY/F");
	T->Branch("Evt_fvtxZ",&Evt_fvtxZ,"Evt_fvtxZ/F");
	T->Branch("Evt_bbcZ",&Evt_bbcZ,"Evt_bbcZ/F");

	T->Branch("from_same_vertex",&from_same_vertex,"from_same_vertex/O");
	T->Branch("Mutr_tr0_closeZ0",&Mutr_tr0_closeZ0,"Mutr_tr0_closeZ0/F");
	T->Branch("Mutr_tr1_closeZ0",&Mutr_tr1_closeZ0,"Mutr_tr1_closeZ0/F");
	T->Branch("Mutr_tr0_charge",&Mutr_tr0_charge,"Mutr_tr0_charge/I");
	T->Branch("Mutr_tr1_charge",&Mutr_tr1_charge,"Mutr_tr1_charge/I");
	T->Branch("nTracklets",&nTracklets,"nTracklets/I");
	T->Branch("nTracklets_all",&nTracklets_all,"nTracklets_all/I");
	T->Branch("nTracklets_samearm",&nTracklets_samearm,"nTracklets_samearm/I");
	T->Branch("nTracklets_samearm_all",&nTracklets_samearm_all,"nTracklets_samearm_all/I");
	T->Branch("bbindex",&bbindex,"bbindex/I");
	T->Branch("ccindex",&ccindex,"ccindex/I");
	T->Branch("flavor_index",&flavor_index,"flavor_index/I");

	}

void dnpfilter::init_histos()
{
//-------------------------------------- histograms ----------------------------------------------------
  sprintf (namest, "Pz_dist");
  Pz_dist    = new TH1F(namest, namest, 600, -150, 150.0);
}

bool dnpfilter::pass_dimuon_cut(const DiMuon *dimuon) const
{
	if(!dimuon) {return false;}

	else if(
			sqrt(dimuon->get_X0()*dimuon->get_X0() + dimuon->get_Y0()*dimuon->get_Y0()) <2.0
			&& dimuon->get_Evt_vtxchi2() < 5.0
			&& dimuon->get_Tr0_DG0() < 20.0
			&& dimuon->get_Tr1_DG0() < 20.0
			&& dimuon->get_Tr0_DDG0() < 10.0
			&& dimuon->get_Tr1_DDG0() < 10.0
			&& dimuon->get_Tr0_ntrhits() > 9
			&& dimuon->get_Tr1_ntrhits() > 9
			&& dimuon->get_Tr0_nidhits() > 5
			&& dimuon->get_Tr1_nidhits() > 5
			&& dimuon->get_Tr0_lastgap() > 3
			&& dimuon->get_Tr1_lastgap() > 3
			&& fabs(dimuon->get_Tr0_dca_r()) < 10.0
			&& fabs(dimuon->get_Tr1_dca_r()) < 10.0
			&& dimuon->get_pT() < 20.0
			&& fabs(dimuon->get_pz()) < 100.0
			&& fabs(dimuon->get_rapidity()) > 1.2
			&& fabs(dimuon->get_rapidity()) < 2.4
			&& dimuon->get_Tr0_pz() * dimuon->get_Tr1_pz() > 0
			&& dimuon->get_same_event() == 1
			&& dimu_container->get_nDiMuons() == 1
			&& dimuon->get_Tr0_trchi2() < 10
			&& dimuon->get_Tr1_trchi2() < 10
			) {return true;}
        else {return false;}
}


//----------------------------------------! global termination-------------------------------------
int dnpfilter::End(PHCompositeNode *topNode)
{

	_out_file->cd();

	T->Write();
	//Pz_dist->Write();

	return 0;
}


dnpfilter::~dnpfilter()
{

	delete T;
	T = NULL;

	_out_file->Close();
}
