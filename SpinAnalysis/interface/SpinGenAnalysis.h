#ifndef __SpinGenAnalysis__
#define __SpinGenAnalysis__

#include "BaseAnalysis.h"
#include "BaseSmearer.h"
#include "PhotonAnalysis/interface/StatAnalysis.h"
#include "RooContainer.h"
#include "VertexAnalysis/interface/HggVertexAnalyzer.h"

#include "EnergySmearer.h"
#include "EfficiencySmearer.h"
#include "DiPhoEfficiencySmearer.h"
#include "KFactorSmearer.h"
#include <iostream>
#include <fstream>
#include "math.h"

// ------------------------------------------------------------------------------------
class SpinGenAnalysis : public StatAnalysis
{
public:
	SpinGenAnalysis();
	virtual ~SpinGenAnalysis();

	virtual const std::string & name() const { return name_; };

	// LoopAll analysis interface implementation
	void GetBranches(TTree *t, std::set<TBranch *>& s );
	bool SkimEvents(LoopAll&, int);
	void Init(LoopAll&);
	void Term(LoopAll&);

	bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category,
	int & diphoton_id,
	bool & isCorrectVertex,
	float &kinematic_bdtout,
	bool isSyst=false,
	float syst_shift=0., bool skipSelection=false,
	BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0);

	void FillRooContainer(LoopAll& l, int cur_type, float mass, float diphotonMVA,
	int category, float weight, bool isCorrectVertex, int diphoton_id);

protected:
	void fillControlPlots(const TLorentzVector & lead_p4, const  TLorentzVector & sublead_p4, const TLorentzVector & Higgs,
	int category, float evweight, LoopAll & l );

};

#endif


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
