#include "TTargEventAction.hh"
#include "TTargRunAction.hh"
#include "TTargAnalysis.hh"
#include "TTargHit.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessTable.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"

TTargEventAction::TTargEventAction(TTargRunAction* runAction)
 : G4UserEventAction(),
   fRunAction(runAction),
   fEdep(0.),
   nBrem(0)
{
}

TTargEventAction::~TTargEventAction()
{
}

void TTargEventAction::BeginOfEventAction(const G4Event*)
{
  fEdep = 0.;
/*  G4Region* region;
  G4String regName;
  G4ProductionCuts* cut;
  regName = "Box2";
  region = G4RegionStore::GetInstance()->GetRegion(regName);
  cut = new G4ProductionCuts;
  cut->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("DarkPhoton"));
  cut->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("e-"));
  cut->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("e+"));
  cut->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("gamma"));
  region->SetProductionCuts(cut);
*/
     G4bool active = true;
     G4String pname = "eDBrem";
     G4ProcessTable* ptable = G4ProcessTable::GetProcessTable();
     ptable->SetProcessActivation(pname, active);
     rstnBrem();
}

void TTargEventAction::EndOfEventAction(const G4Event* event)
{
   auto anMan = G4AnalysisManager::Instance();
   auto hc = event->GetHCofThisEvent()->GetHC(1);
   auto hit = static_cast<TTargHit*>(hc->GetHit(0));
   G4double energy = hit->GetEnergy();
   anMan->FillH1(0,energy,1);
}


TTargRunAction* TTargEventAction::GetRunAction() const
{
   return fRunAction;
}

