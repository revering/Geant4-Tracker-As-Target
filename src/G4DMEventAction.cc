#include "G4DMEventAction.hh"
#include "G4DMRunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4ProcessTable.hh"

G4DMEventAction::G4DMEventAction(G4DMRunAction* runAction)
 : G4UserEventAction(),
   fRunAction(runAction),
   fEdep(0.),
   nBrem(0)
{
}

G4DMEventAction::~G4DMEventAction()
{
}

void G4DMEventAction::BeginOfEventAction(const G4Event*)
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

void G4DMEventAction::EndOfEventAction(const G4Event*)
{
}


G4DMRunAction* G4DMEventAction::GetRunAction() const
{
   return fRunAction;
}

