#include "TTargSteppingAction.hh"
#include "TTargEventAction.hh"
#include "TTargDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4DarkPhoton.hh"
#include "G4ProcessTable.hh"

TTargSteppingAction::TTargSteppingAction(TTargEventAction* eventAction)
 : G4UserSteppingAction(),
   fEventAction(eventAction),
   fScoringVolume(0),
   NEmissions(0)
{
}

TTargSteppingAction::~TTargSteppingAction()
{
}

void TTargSteppingAction::UserSteppingAction(const G4Step* step)
{
   if((fEventAction->GetnBrem())==0)
   {
      G4Track* track = step->GetTrack();
      const G4DynamicParticle* p = track->GetDynamicParticle();
      const G4ParticleDefinition* theParticleDefinition = p->GetParticleDefinition();
      if( theParticleDefinition == G4DarkPhoton::DarkPhoton())
      {
/*        G4Region* region;
	G4String regName;
	G4ProductionCuts* cut;
        regName = "Box2";
	region = G4RegionStore::GetInstance()->GetRegion(regName);
	cut = new G4ProductionCuts;
	cut->SetProductionCut(1000000*CLHEP::mm,G4ProductionCuts::GetIndex("DarkPhoton"));
        cut->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("e-"));
        cut->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("e+"));
        cut->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("gamma"));
	region->SetProductionCuts(cut);
*/
         G4bool active = false;
         G4String pname = "eDBrem";
         G4ProcessTable* ptable = G4ProcessTable::GetProcessTable();
         ptable->SetProcessActivation(pname, active); 
	 fEventAction->incnBrem();
      }
   }
/*   G4Track* track = step->GetTrack();
   const G4DynamicParticle* p = track->GetDynamicParticle();
   const G4ParticleDefinition* theParticleDefinition = p->GetParticleDefinition();
   if ( theParticleDefinition == G4Electron::ElectronDefinition() ||
        theParticleDefinition == G4Positron::PositronDefinition() )
   {
      G4double ekin = track->GetKineticEnergy();
      G4double DensityMat = track->GetMaterial()->GetDensity()/(g/cm3);
      if( fEventAction->GetRunAction()->GetDarkPhoton()->Emission(ekin, DensityMat, step->GetStepLength()) ) 
      {
         if(NEmissions)
         {
            G4cout << "Attempt to emit more than one A in the event, E = " << ekin << G4endl;
         }
         else
         {
            NEmissions++;
            G4double XAcc = fEventAction->GetRunAction()->GetDarkPhoton()->SimulateEmission(ekin);
            track->SetKineticEnergy(ekin*(1.-XAcc)*GeV);
         }
      }
   }
*/
}
