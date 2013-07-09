#include "../interface/SpinGenAnalysis.h"

#include "../interface/JetHandler.h"

#include "Sorters.h"
#include "PhotonReducedInfo.h"
#include <iostream>
#include <algorithm>

#include "CMGTools/External/interface/PileupJetIdentifier.h"

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

	return true;
}

// ----------------------------------------------------------------------------------------------------
void SpinGenAnalysis::FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA,
int category, float weight, bool isCorrectVertex, int diphoton_id)
{
	l.FillTree("run",l.run);
	l.FillTree("lumis",l.lumis);
	l.FillTree("event",l.event);
	l.FillTree("mass",mass);
	l.FillTree("weight",weight);
	l.FillTree("category",category);
	l.FillTree("diphotonMVA",diphotonMVA);
	l.FillTree("vbfMVA",myVBF_MVA);
	l.FillTree("VBFevent", VBFevent);

	l.FillTree("sampleType",cur_type);

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

	float mass = Higgs.M();
	l.FillHist("all_mass",category+1, Higgs.M(), evweight);

	if( mass>=massMin && mass<=massMax  )
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

		if (VBFevent)
		{
			float myweight =  1;
			float sampleweight = l.sampleContainer[l.current_sample_index].weight();
			if(evweight*sampleweight!=0) myweight=evweight/sampleweight;

			l.FillCutPlots(category+1,1,"_sequential",evweight,myweight);
		}
	}
}

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
