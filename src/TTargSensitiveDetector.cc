#include "G4VSensitiveDetector.hh"
#include "TTargSensitiveDetector.hh"
#include "TTargHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"

TTargSensitiveDetector::TTargSensitiveDetector(G4String name)
: G4VSensitiveDetector(name),
  fHitsCollection(nullptr), fHCID(-1)
{
   collectionName.insert("HitsColl");
}

TTargSensitiveDetector::~TTargSensitiveDetector()
{}

void TTargSensitiveDetector::Initialize(G4HCofThisEvent* hce)
{
   fHitsCollection = new TTargHitsCollection(SensitiveDetectorName, collectionName[0]);

   if (fHCID<0)
   {
      fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
   }

   hce->AddHitsCollection(fHCID,fHitsCollection);
   G4cout<<fHCID <<"\n";
}

G4bool TTargSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*)
{
   const G4Track* track = step->GetTrack();
   G4cout << "Found a hit, particle name is: " << track->GetDefinition()->GetParticleName() << "\n";
//   if (track->GetDefinition()->GetParticleName() == "e-")
//   {
      G4double e = track->GetKineticEnergy();
      G4ThreeVector m = track->GetMomentumDirection();
      G4String Id = track->GetDefinition()->GetParticleName();
//      G4double pID = track->GetDefinition()->GetParticleID();
      auto hit = new TTargHit(e, m, Id);
      fHitsCollection->insert(hit);
  // }


   return true;
}
