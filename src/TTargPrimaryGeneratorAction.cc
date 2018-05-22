#include "TTargPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

TTargPrimaryGeneratorAction::TTargPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
   G4int n_particle = 1;
   fParticleGun = new G4ParticleGun(n_particle);

   G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
   G4ParticleDefinition* particle = particleTable->FindParticle("e-");
   fParticleGun->SetParticleDefinition(particle);
   fParticleGun->SetParticlePosition(G4ThreeVector(-100.,0.,0.));
   fParticleGun->SetParticleEnergy(4*GeV);
   fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,1.,0.));
}

TTargPrimaryGeneratorAction::~TTargPrimaryGeneratorAction()
{
   delete fParticleGun;
}

void TTargPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

   fParticleGun->SetParticlePosition(G4ThreeVector(-50*cm,-110*cm,0));
   fParticleGun->GeneratePrimaryVertex(anEvent);

   fParticleGun->SetParticlePosition(G4ThreeVector(50*cm,-110*cm,0));
   fParticleGun->GeneratePrimaryVertex(anEvent);
}
