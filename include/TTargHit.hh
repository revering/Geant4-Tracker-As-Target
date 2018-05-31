#ifndef TTargHit_h
#define TTargHit_h

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"

class TTargHit : public G4VHit
{
   public:
      TTargHit();
      TTargHit(G4double e, G4ThreeVector m, G4String pID);
      virtual ~TTargHit();
      G4double GetEnergy() {return energy;}
      G4String GetParticleId() {return particleID;}

   private:
      G4ThreeVector mom;
      G4double energy;
      G4String particleID;

};

using TTargHitsCollection = G4THitsCollection<TTargHit>;

#endif
