#include "G4muDarkBremsstrahlung.hh"
#include "G4SystemOfUnits.hh"
#include "G4DarkPhoton.hh"
#include "G4muDarkBremsstrahlungModel.hh"
#include "G4UnitsTable.hh"
#include "G4LossTableManager.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"

using namespace std;

G4muDarkBremsstrahlung::G4muDarkBremsstrahlung(const G4String& name):
   G4VEnergyLossProcess(name),
   isInitialised(false)
{  
   G4int subtype = 64;   
   SetProcessSubType(subtype);
   SetSecondaryParticle(G4DarkPhoton::DarkPhoton());
   SetIonisation(false);
}

G4muDarkBremsstrahlung::~G4muDarkBremsstrahlung()
{}

G4bool G4muDarkBremsstrahlung::IsApplicable(const G4ParticleDefinition& p)
{
   return (&p == G4MuonPlus::MuonPlus() || &p == G4MuonMinus::MuonMinus());
}

void G4muDarkBremsstrahlung::InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                                    const G4ParticleDefinition*)
{
   if(!isInitialised)
   {  
      SetEmModel(new G4muDarkBremsstrahlungModel(),1);

      G4double energyLimit = 1*GeV;

      EmModel(1)->SetLowEnergyLimit(MinKinEnergy());
      EmModel(1)->SetHighEnergyLimit(energyLimit);
        
      G4VEmFluctuationModel* fm = 0;
      AddEmModel(1, EmModel(1), fm);

      isInitialised = true;
   }

//   G4LossTableManager* man = G4LossTableManager::Instance();
   G4double eth = 1*MeV;
   EmModel(1)->SetSecondaryThreshold(eth);
   EmModel(1)->SetLPMFlag(false);
}

void G4muDarkBremsstrahlung::PrintInfo()
{
   if(EmModel(1))
   {
//      G4LossTableManager* man = G4LossTableManager::Instance();
      G4cout << "    LPM flag: " << "false " << " for E > " << EmModel(1)->HighEnergyLimit()/GeV<< " GeV";
      G4cout << G4endl;
   }
}


