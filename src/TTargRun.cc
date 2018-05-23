#include "TTargRun.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

TTargRun::TTargRun() : G4Run()
{;}

TTargRun::~TTargRun()
{;}

void TTargRun::RecordEvent(const G4Event* evt)
{
   nEvent++;
//   G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
//   G4THitsMap<G4double>* f_pt = (G4THitsMap<G4double>*)(HCE->GetHC(f_pt_ID));
//   G4THitsMap<G4double>* i_pt = (G4THitsMap<G4double>*)(HCE->GetHC(i_pt_ID));
}



