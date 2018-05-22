#ifndef G4DarkPhoton_h
#define G4DarkPhoton_h

#include "G4ParticleDefinition.hh"
#include "globals.hh"

#include "CLHEP/Units/SystemOfUnits.h"

class G4DarkPhoton : public G4ParticleDefinition
{
   private:
      static G4DarkPhoton* theDarkPhoton;

      G4DarkPhoton(
         const G4String&      Name,
	 G4double             mass,
	 G4double             width,
	 G4double             charge,
	 G4int                iSpin,
	 G4int                iParity,
	 G4int                iConjugation,
	 G4int                iIsospin,
	 G4int                iIsospin3,
	 G4int                gParity,
	 const G4String&      pType,
	 G4int                lepton,
	 G4int                baryon,
	 G4int                encoding,
	 G4bool               stable,
	 G4double             lifetime,
	 G4DecayTable         *decaytable
	 );

      virtual ~G4DarkPhoton();

   public:

      static G4DarkPhoton* DarkPhoton();
};


#endif


