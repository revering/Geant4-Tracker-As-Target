#include "TTargHit.hh"
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"

TTargHit::TTargHit() 
: G4VHit(), 
  mom(0), 
  energy(0.), 
  particleID("Nothing")
{}

TTargHit::TTargHit(G4double e, G4ThreeVector m, G4String pID) 
: G4VHit(), 
  mom(m), 
  energy(e),
  particleID(pID)
{}

TTargHit::~TTargHit()
{}

