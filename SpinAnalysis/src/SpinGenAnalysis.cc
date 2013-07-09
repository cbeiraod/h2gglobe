#include "../interface/SpinGenAnalysis.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#define PADEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
SpinGenAnalysis::SpinGenAnalysis()
{
	name_ = "SpinGenAnalysis";
}

// ----------------------------------------------------------------------------------------------------
SpinGenAnalysis::~SpinGenAnalysis()
{
}

// ----------------------------------------------------------------------------------------------------
void SpinGenAnalysis::Term(LoopAll& l)
{
}

// ----------------------------------------------------------------------------------------------------
void SpinGenAnalysis::Init(LoopAll& l)
{
	doSystematics = false;
	StatAnalysis::Init(l);
}

// ----------------------------------------------------------------------------------------------------
bool SpinGenAnalysis::SkimEvents(LoopAll& l, int jentry)
{
	return true;
}

// ----------------------------------------------------------------------------------------------------
void SpinGenAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s )
{
}

// ----------------------------------------------------------------------------------------------------
bool SpinGenAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight,
		int & category, int & diphoton_id,
		bool & isCorrectVertex,
		float &kinematic_bdtout,
		bool isSyst,
		float syst_shift, bool skipSelection,
		BaseGenLevelSmearer *genSys, BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys)
{
	assert( isSyst || !skipSelection );

	int cur_type = l.itype[l.current];
	float sampleweight = l.sampleContainer[l.current_sample_index].weight();
	
	TLorentzVector* part1 = (TLorentzVector*) l.gh_glu1_p4->At(0);
	TLorentzVector* part2 = (TLorentzVector*) l.gh_glu2_p4->At(0);
	TLorentzVector* pho1  = (TLorentzVector*) l.gh_pho1_p4->At(0);
	TLorentzVector* pho2  = (TLorentzVector*) l.gh_pho2_p4->At(0);
	
	TLorentzVector diphoton = *pho1 + *pho2;
	
	TLorentzVector *pho1_boosted, *pho2_boosted, *part1_boosted, *part2_boosted;
	pho1_boosted  = (TLorentzVector*) pho1->Clone();
	pho2_boosted  = (TLorentzVector*) pho2->Clone();
	part1_boosted = (TLorentzVector*)part1->Clone();
	part2_boosted = (TLorentzVector*)part2->Clone();
	
	TVector3 boost = diphoton.BoostVector();
	pho1_boosted->Boost(-boost);
	pho2_boosted->Boost(-boost);
	part1_boosted->Boost(-boost);
	part2_boosted->Boost(-boost);
	
	TVector3 gg_SQA = part1_boosted->Vect() - part2_boosted->Vect();
	gg_SQA = gg_SQA.Unit();
	
	{
		TVector3 q = pho1_boosted->Vect();
		Double_t ptot2 = gg_SQA.Mag2()*q.Mag2();
		if(ptot2 == 0)
			gg_costh = 0;
		else
		{
			Double_t arg = gg_SQA.Dot(q)/TMath::Sqrt(ptot2);
			if(arg >  1.0) arg =  1.0;
			if(arg < -1.0) arg = -1.0;
			gg_costh = arg;
		}
	}

	//gg_costh = gg_SQA.angle(pho1_boosted.Vect());

	//gg_costh should follow a distribution of the type: 5/32 (1 + 6 costh^2 + costh^4)
	//we want to reweight it to: 3/8 (1 + costh^2) && 5/12 (1 + costh^4)
	// weight1 = (3/8 (1 + costh^2)) / (5/32 (1 + 6 costh^2 + costh^4))
	// weight2 = (5/12 (1 + costh^4)) / (5/32 (1 + 6 costh^2 + costh^4))
	
	double costh2 = gg_costh*gg_costh;
	double costh4 = costh2*costh2;
	weight1 = 12./5. * (1. + costh2)/(1. + 6*costh2 + costh4);
	weight2 = 8./3. * (1. + costh4)/(1. + 6*costh2 + costh4);

	return true;
}

// ----------------------------------------------------------------------------------------------------
void SpinGenAnalysis::FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA,
		int category, float weight, bool isCorrectVertex, int diphoton_id)
{
	l.FillTree("run",l.run);
	l.FillTree("event",l.event);
	l.FillTree("lumis",l.lumis);
	l.FillTree("CMS_hgg_mass",mass);
	l.FillTree("weight",weight);
	l.FillTree("category",category);
	l.FillTree("diphotonMVA",diphotonMVA);

	l.FillTree("sampleType",cur_type);
	
	l.FillTree("costh", gg_costh);
	l.FillTree("spin_weight1", weight1);
	l.FillTree("spin_weight2", weight2);

	TLorentzVector lead_p4, sublead_p4, Higgs;
	float lead_r9 = 0., sublead_r9 = 0.;
	TVector3 * vtx;

	l.FillTree("leadPt",(float)lead_p4.Pt());
	l.FillTree("subleadPt",(float)sublead_p4.Pt());
	l.FillTree("leadR9",lead_r9);
	l.FillTree("subleadR9",sublead_r9);

	SpinGenAnalysis::fillControlPlots(lead_p4, sublead_p4, Higgs, category, weight, l);
}

// ----------------------------------------------------------------------------------------------------
void SpinGenAnalysis::fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs,
int category, float evweight, LoopAll & l )
{
	// control plots
	if( category>=0 )
	{
		fillControlPlots( lead_p4, sublead_p4, Higgs, -1, evweight, l );
	}

	l.FillHist("all_mass",category+1, Higgs.M(), evweight);

	if( Higgs.M()>=massMin && Higgs.M()<=massMax  )
	{
		l.FillHist("mass",category+1, Higgs.M(), evweight);
		l.FillHist("eta",category+1, Higgs.Eta(), evweight);
		l.FillHist("pt",category+1, Higgs.Pt(), evweight);

		l.FillHist("pho_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho1_pt",category+1,lead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,lead_p4.Eta(), evweight);
		l.FillHist("pho1_eta",category+1,lead_p4.Eta(), evweight);

		l.FillHist("pho_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho2_pt",category+1,sublead_p4.Pt(), evweight);
		l.FillHist("pho_eta",category+1,sublead_p4.Eta(), evweight);
		l.FillHist("pho2_eta",category+1,sublead_p4.Eta(), evweight);
	}
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
