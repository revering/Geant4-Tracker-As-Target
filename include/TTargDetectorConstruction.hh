#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4PhysicalVolume;
class G4LogicalVolume;

class DMDetectorConstruction : public G4VUserDetectorConstruction
{
   public:
      DMDetectorConstruction();
      virtual ~DMDetectorConstruction();

      virtual G4VPhysicalVolume* Construct();

      G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

    protected:
      G4LogicalVolume* fScoringVolume;
};


#endif
