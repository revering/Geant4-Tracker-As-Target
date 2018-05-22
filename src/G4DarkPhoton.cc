#include "G4DarkPhoton.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


G4DarkPhoton* G4DarkPhoton::theDarkPhoton = 0;

G4DarkPhoton::G4DarkPhoton(
   const G4String&     aName, 
   G4double            mass, 
   G4double            width,       
   G4double            charge,
   G4int               iSpin, 
   G4int               iParity,   
   G4int               iConjugation, 
   G4int               iIsospin, 
   G4int               iIsospin3,
   G4int               gParity,
   const G4String&     pType,
   G4int               lepton,
   G4int               baryon,
   G4int               encoding,
   G4bool              stable, 
   G4double            lifetime,
   G4DecayTable        *decaytable)	
   : G4ParticleDefinition( aName, mass, width, charge, iSpin, iParity, 
                           iConjugation, iIsospin, iIsospin3, gParity, pType,    
			   lepton, baryon, encoding, stable, lifetime, decaytable ) 
{

}

G4DarkPhoton::~G4DarkPhoton()
{

}

G4DarkPhoton* G4DarkPhoton::DarkPhoton()
{
   if(!theDarkPhoton) {
     
      const G4String&     name = "DarkPhoton";
      G4double            mass = 100.*MeV; 
      G4double            width = 0.;       
      G4double            charge = 0;
      G4int               iSpin = 0;
      G4int               iParity = 0; 
      G4int               iConjugation = 0; 
      G4int               iIsospin = 0;
      G4int               iIsospin3 = 0;
      G4int               gParity = 0;
      const G4String&     pType = "boson";
      G4int               lepton = 0;
      G4int               baryon = 0;
      G4int               encoding = 0;
      G4bool              stable = true;
      G4double            lifetime = -1.0;
      G4DecayTable*        decaytable = 0;

      theDarkPhoton = new G4DarkPhoton(
         name, mass, width, charge, iSpin, iParity,
	 iConjugation, iIsospin, iIsospin3, gParity, pType,
	 lepton, baryon, encoding, stable, lifetime, decaytable);
   }
   return theDarkPhoton;
}


