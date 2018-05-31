#ifndef TTargSensitiveDetector_h
#define TTargSensitiveDetector_h 

#include "G4VSensitiveDetector.hh"
#include "TTargHit.hh"


class TTargSensitiveDetector : public G4VSensitiveDetector
{
   public:
      TTargSensitiveDetector(G4String name);
      ~TTargSensitiveDetector();

   public:
      G4bool ProcessHits(G4Step* step, G4TouchableHistory* R0hist);
      void Initialize(G4HCofThisEvent* HCE);

   private:
      TTargHitsCollection* fHitsCollection;
      G4int fHCID;

};

#endif
