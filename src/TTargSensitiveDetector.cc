#include "G4VSensitiveDetector.hh"
#include "TTargSensitiveDetector.hh"
#include "TTargHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"

TTargSensitiveDetector::TTargSensitiveDetector(G4String name)
: G4VSensitiveDetector(name),
  fHitsCollection(nullptr)
{
   collectionName.insert("HitsColl");
}

TTargSensitiveDetector::~TTargSensitiveDetector()
{}

void TTargSensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
   fHitsCollection = new TTargHitsCollection(SensitiveDetectorName, collectionName[0]);

   hce->AddHitsCollection(1,fHitsCollection);
}

G4bool TTargSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
   const G4Track* track = step->GetTrack();
   if (track->GetDefinition()->GetParticleName() == "e-")
   {
      G4double e = track->GetKineticEnergy();
      G4ThreeVector m = track->GetMomentumDirection();
//      G4double pID = track->GetDefinition()->GetParticleID();
      auto hit = new TTargHit(e, m, 1);
      fHitsCollection->insert(hit);
   }

   return true;
}
