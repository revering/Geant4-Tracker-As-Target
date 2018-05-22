#include "TTargPhysicsList.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DarkPhoton.hh"
#include "G4PhysicsListHelper.hh"
#include "G4Electron.hh"
#include "G4eBremsstrahlung.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessTable.hh"
#include "G4Gamma.hh"
#include "G4eDarkBremsstrahlung.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4muDarkBremsstrahlung.hh"


DMPhysicsList::DMPhysicsList() : G4VModularPhysicsList(),
   fEmPhysicsList(0)
{
   fEmPhysicsList = new G4EmStandardPhysics();
}

DMPhysicsList::~DMPhysicsList()
{
   delete fEmPhysicsList;
}

void DMPhysicsList::SetCuts()
{
   G4VUserPhysicsList::SetCuts();
   SetCutValue(100000000*CLHEP::mm,"DarkPhoton");
 

   G4Region* region;
   G4String regName;
   G4ProductionCuts* cuts;
 
   cuts = new G4ProductionCuts;
//   cuts->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("e-"));
//   cuts->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("e+"));
//   cuts->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("gamma"));

   regName = "Tracker";
   region = G4RegionStore::GetInstance()->GetRegion(regName);
   cuts->SetProductionCut(0.01*CLHEP::mm,G4ProductionCuts::GetIndex("DarkPhoton"));
   region->SetProductionCuts(cuts);

}

void DMPhysicsList::ConstructParticle()
{
   G4DarkPhoton::DarkPhoton();
   fEmPhysicsList->ConstructParticle();

}


void DMPhysicsList::ConstructProcess()
{
   AddTransportation();
   fEmPhysicsList->ConstructProcess();
   G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
   G4ParticleDefinition* particle = G4Electron::ElectronDefinition();
   ph->RegisterProcess(new G4eDarkBremsstrahlung(), particle);
   G4ParticleDefinition* muon = G4MuonMinus::MuonMinusDefinition();
   ph->RegisterProcess(new G4muDarkBremsstrahlung(), muon);
   G4ParticleDefinition* muplus = G4MuonPlus::MuonPlusDefinition();
   ph->RegisterProcess(new G4muDarkBremsstrahlung(), muplus);
//   ph->RegisterProcess(new G4eBremsstrahlung(), particle);
}

