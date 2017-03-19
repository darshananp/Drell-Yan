// $Id: saModuleSimpleDimuon.C created by jinhuang modified for Drell-Yan Analysis By Darshana Perera 
#include <Fun4AllReturnCodes.h>
#include <PHCompositeNode.h>
#include <getClass.h>
#include <MWGConsts.h>
#include <Tools.h>
#include <Fun4AllServer.h>
#include <Fun4AllHistoManager.h>
#include <TH1F.h>
#include <TriggerHelper.h>
#include <TrigRunLvl1.h>
#include <TrigLvl1.h>
#include <TH3.h>
#include <TH2.h>
#include <TF1.h>
#include <SyncObject.h>
#include <RunHeader.h>
#include <TMath.h>
#include <SingleMuonContainer.h>
#include <DiMuonContainer.h>
#include <DiMuon.h>
#include <MCDiMuonContainer.h>
#include <MCDiMuon.h>
#include <TFvtxCompactTrk.h>
#include <TFvtxCompactTrkMap.h>
#include "saFlag.h"
#include <SingleMuon.h>
#include <Fun4AllReturnCodes.h>
#include <PHNodeIterator.h>
#include <getClass.h>
#include <TClonesArray.h>
#include "saModuleSimpleDimuon.h"

using namespace std;

saModuleSimpleDimuon::saModuleSimpleDimuon(const std::string &name) :saModuleBase(name), _h_mass(NULL)
{

}

saModuleSimpleDimuon::~saModuleSimpleDimuon()
{
 	_use_cut_tracklet_chi2=false;
	_cut_tracklet_chi2=4.;
	_cut_tracklet_dcar=1.5;
	_use_2_hit_tracklet=false;
}

//! global initialization
int saModuleSimpleDimuon::init(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  // make a simple histogram //the template histogram, any TH1 and derivatives is accepted
  TH1F * h_mass = new TH1F("InvMass", "Invariant Mass;Invariant Mass (GeV)", 2,4, 8);
  _h_mass = new saHist(h_mass, saHist::DEFAULT_FLAG | saHist::RUN_PROPERTY | saHist::FILL_PROPERTY| saHist::CROSSING_CHECKS,Verbosity
  ());
  hm.get()->registerHisto(_h_mass);

  // make the flag node which is syncronized with cut on picodst_object
  PHCompositeNode* dstNode = NULL;
  PHNodeIterator nodeItr(topNode);
  dstNode = static_cast<PHCompositeNode*>(nodeItr.findFirst("PHCompositeNode",
      "DST"));
  if (!dstNode)
    {
      dstNode = new PHCompositeNode("DST");
      topNode->addNode(dstNode);
    }
  saFlagC *flags = new saFlagC();
  if (flags)
    {

      // make a new flag node called DST/SimpleDimuonFlag
      PHIODataNode<PHObject> * node = new PHIODataNode<PHObject>(flags,
          "SimpleDimuonFlag", "PHObject");

      if (!node)
        {

          cout
              << "saModuleSimpleDimuon::Init failed to create saEventProperty Node"
              << endl;
          return ABORTRUN;

        }
      else
        {

          dstNode->addNode(node);
          cout << "saFlag Node is added with version " << flags->ClassName()
              << " as " << node->getName() << endl;
        }
    }
  else
    {
      cout << "saModuleSimpleDimuon::Init failed to create saEventProperty"
          << endl;
      return ABORTRUN;
    }

  return EVENT_OK;
}

//! Run initialization
int saModuleSimpleDimuon::init_run(PHCompositeNode *topNode,sa_hist_mangager_ptr hm)
{
  return EVENT_OK;
}

//! event method
int saModuleSimpleDimuon::event(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{
//---------------------------------------------------------------------------------------------------------
   dimuons = findNode::getClass<DiMuonContainer>(topNode,"DiMuonContainer");
  if (!dimuons)
    {
      cout<< "saModuleSimpleDimuon:: DiMuonContainer - ERROR - not in Node Tree"<< endl;
      return ABORTRUN;
    }
//---------------------------------------------------------------------------------------------------------
  const saEventProperty * ep = findNode::getClass<saEventProperty>(topNode,"saEventProperty");
  if (!ep)
    {
      cout<< "saModuleSimpleDimuon::event - ERROR  - Cannot find EventProperty node in Top Node"<< endl;
      return ABORTRUN;
    }
//---------------------------------------------------------------------------------------------------------

	singlemuoncontainer = findNode::getClass<SingleMuonContainer>(topNode,"SingleMuonContainer");
	if (!singlemuoncontainer){
		cout << "dnpfilter:: SingleMuonContainer not in Node Tree" << endl;
		return ABORTRUN;
	}

//-------------------------------------------------------------------------------------------
	fvtx_trk_map = findNode::getClass<DSTReader>(topNode,"DSTReader");
	if(!fvtx_trk_map){
		cout << PHWHERE << "dnpfilter:: DSTReader not in Node Tree" << endl;
		return ABORTRUN;
	}
//------------------------------------------------------------------------------------------------
  if (Verbosity() >= 2)
    {
      cout << "saModuleSimpleDimuon::event - INFO  - ";
      ep->identify();
    }

//---------------------------------------------------------------------------------------------------------
  saFlagC *flag = findNode::getClass<saFlagC>(topNode, "SimpleDimuonFlag");
  if (!flag)
    {
      cout << "saModuleSimpleDimuon::event - SimpleDimuonFlag not in Node Tree"<< endl;
      return ABORTRUN;
    }
//---------------------------------------------------------------------------------------------------------
  // syncronizely set the cut flag with dimuon container
  flag->flags.Set(dimuons->get_nDiMuons());
  flag->flags.Reset(0);
//---------------------------------------------------------------------------------------------------------

  bool save_event = false;

  double Evt_bbcZ = dimuons->get_Evt_bbcZ();

  if ( fabs(Evt_bbcZ) >= 10 ) return DISCARDEVENT;
    
  for (unsigned int imuon = 0; imuon < dimuons->get_nDiMuons(); imuon++)
    {

      //DiMuon *dimuon = dynamic_cast<DiMuon*> (dimu_container->get_DiMuon(idimu));

      dimuon = dynamic_cast<DiMuon*> (dimuons->get_DiMuon(imuon));

      nTracklets = GetNumberofTracklets(fvtx_trk_map,dimuon);

      if ( !passesCuts() ) continue;

      if (nTracklets >= 10) continue;
     
      cout<< "*****************Number of tracklets = "<<nTracklets<<" **************Mass = "<<dimuon->get_mass()<<endl;
      
      save_event = true;
      flag->flags[imuon] = 1;
      _h_mass->Fill(ep, dimuon->get_mass());

    }

  return save_event ? EVENT_OK : DISCARDEVENT;
}


bool saModuleSimpleDimuon::passesCuts() const{
  // Cuts

  if ( dimuon->get_mass() >= 8 || dimuon->get_mass() <= 4 )
    return false;

  if ( !dimuon->get_same_event() ) 
    return false;

  int arm0 = (dimuon->get_Tr0_pz() > 0) ? 1 : 0;
  int arm1 = (dimuon->get_Tr1_pz() > 0) ? 1 : 0;

  if ( arm0 != arm1 ) 
    return false;
      
  if ( dimuon->get_Evt_vtxchi2() >= 4 ) 
    return false;

  if ( dimuon->get_charge() != 0 ) 
    return false;

  if ( fabs(dimuon->get_rapidity()) <= 1.2 || fabs(dimuon->get_rapidity()) >= 2.4 )
    return false;

  if ( fabs(dimuon->get_Tr0_pz()) <= 2 || fabs(dimuon->get_Tr1_pz()) <= 2.0 )
    return false;

  if ( dimuon->get_Tr0_DG0() >= 20 || dimuon->get_Tr1_DG0() >= 20 )
    return false;

  if ( dimuon->get_Tr0_DDG0() >= 8 || dimuon->get_Tr1_DDG0() >= 8 )
    return false;

  if ( dimuon->get_Tr0_ntrhits() <= 10 || dimuon->get_Tr1_ntrhits() <= 10 )
    return false;
  if ( dimuon->get_Tr0_nidhits() <= 6 || dimuon->get_Tr1_nidhits() <= 6 )
    return false;
      
  if ( dimuon->get_Tr0_lastgap() <= 3 || dimuon->get_Tr1_lastgap() <= 3 ) 
    return false;

  if ( dimuon->get_Tr0_trchi2() >= 10 || dimuon->get_Tr1_trchi2() >= 10 )
    return false;


  return true;
}


int saModuleSimpleDimuon::GetNumberofTracklets(DSTReader *fvtx_trk_map, const DiMuon *dimuon)
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
	//cout<<"No tracklet"<<__LINE__<<"size: "<<array->GetSize()<<endl;
	break;
	}
	
	//if(_use_cut_tracklet_chi2 && (tracklet->get_chi2_ndf() > _cut_tracklet_chi2)) continue;
	if(!_use_2_hit_tracklet && tracklet->get_nhits() <= 2) continue;

	SingleMuon *muon0 = singlemuoncontainer->get_SingleMuon(0);
	SingleMuon *muon1 = singlemuoncontainer->get_SingleMuon(1);

	float xx0 = tracklet->get_fvtx_vtx().getX()-(tracklet->get_fvtx_vtx().getZ()- dimuons->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*cos(tracklet->get_fvtx_phi());
	float yy0 = tracklet->get_fvtx_vtx().getY()-(tracklet->get_fvtx_vtx().getZ()- dimuons->get_Evt_fvtxZ())*
	tan(tracklet->get_fvtx_theta())*sin(tracklet->get_fvtx_phi());
	float x0y0=sqrt((xx0 - dimuons->get_Evt_fvtxX())*(xx0 - dimuons->get_Evt_fvtxX()) +
	(yy0 - dimuons->get_Evt_fvtxY())*(yy0 - dimuons->get_Evt_fvtxY()));

	if (fabs(TMath::ATan2(sqrt(muon0->get_px_fvtxmutr()*muon0->get_px_fvtxmutr()+
	muon0->get_py_fvtxmutr()*muon0->get_py_fvtxmutr()),muon0->get_pz_fvtxmutr())-
	tracklet->get_fvtx_theta()) >0.001 && fabs(TMath::ATan2(sqrt(muon1->get_px_fvtxmutr()*
	muon1->get_px_fvtxmutr()+muon1->get_py_fvtxmutr()*muon1->get_py_fvtxmutr()),
	muon1->get_pz_fvtxmutr())-tracklet->get_fvtx_theta()) >0.001 && fabs(TMath::ATan2(muon0->get_py_fvtxmutr(),
	muon0->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && fabs(TMath::ATan2(muon1->get_py_fvtxmutr(),
	muon1->get_px_fvtxmutr())-tracklet->get_fvtx_phi())>0.001 && tracklet->get_fvtx_theta()+0!=0 && 
	x0y0 < 1.5){
	ntrklets++;
	}

	}

	return ntrklets;
        
}


//! global termination
int saModuleSimpleDimuon::end(PHCompositeNode *topNode, sa_hist_mangager_ptr hm)
{

  //calculate asymmetry with relative lumi of BbcNoCutPHENIX
  _h_mass->CalcAsymmetry(hm.get()->getHisto_DefaultLumi(),saHist::BbcNoCutPHENIX);

  return EVENT_OK;
}

