#include "TTargDetectorConstruction.hh"
#include "TTargSensitiveDetector.hh"

#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SDManager.hh"

DMDetectorConstruction::DMDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{  }

DMDetectorConstruction::~DMDetectorConstruction()
{  }

G4VPhysicalVolume* DMDetectorConstruction::Construct()
{
   G4NistManager* nist = G4NistManager::Instance();

   G4double env_sizeXY = 3*cm, env_sizeZ = 3*cm;
   G4Material* env_mat = nist->FindOrBuildMaterial("G4_Si");
   G4double si_thickness = 600*um;
   G4double sp_thickness = 0.5*cm;
   G4bool checkOverlaps = true;

   G4double world_sizeXY = 1.2*env_sizeXY;
   G4double world_sizeZ = 1.2*env_sizeZ;
   G4double density = 0.001*g/cm3;
   G4double a =  14.007*g/mole;
   G4Material* world_mat = new G4Material("vacuum", 7., a, density);


   G4Box* SolidWorld = new G4Box("World", world_sizeXY, world_sizeXY, world_sizeZ);
   G4LogicalVolume* LogicWorld = new G4LogicalVolume(SolidWorld, world_mat, "World");
   G4VPhysicalVolume* physWorld = new G4PVPlacement( 0, G4ThreeVector(), LogicWorld, "World", 0, false, 0, checkOverlaps);

   G4Box* detector = new G4Box("Detector", env_sizeXY, env_sizeXY, 0.5*sp_thickness);
   G4LogicalVolume* ldetector = new G4LogicalVolume(detector, world_mat, "Detector");
   new G4PVPlacement(0,G4ThreeVector(0,0,si_thickness+1.1*sp_thickness), ldetector, "Detector", LogicWorld, false, 0, checkOverlaps);

   G4Box* predetector = new G4Box("Predetector", env_sizeXY, env_sizeXY, 0.5*sp_thickness);
   G4LogicalVolume* lpredetector = new G4LogicalVolume(predetector, world_mat, "Predetector");
   new G4PVPlacement(0,G4ThreeVector(0,0,-si_thickness-1.1*sp_thickness), lpredetector, "Predetector", LogicWorld, false, 0, checkOverlaps);

   G4Box* solidEnv = new G4Box("Target", env_sizeXY, env_sizeXY, si_thickness);
   G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv, env_mat, "Target");
   new G4PVPlacement(0,G4ThreeVector(0,0,0), logicEnv, "Target", LogicWorld, false, 0, checkOverlaps);
   const G4String& rName = "Tracker";
   G4Region* tracker = new G4Region(rName);
   logicEnv->SetRegion(tracker);
   tracker->AddRootLogicalVolume(logicEnv);

   auto sdManager = G4SDManager::GetSDMpointer();
   G4String SDname;
   auto scorer = new TTargSensitiveDetector(SDname="/scorer");
   sdManager->AddNewDetector(scorer);
   ldetector->SetSensitiveDetector(scorer);
   auto prescorer = new TTargSensitiveDetector(SDname="/prescorer");
   sdManager->AddNewDetector(prescorer);
   lpredetector->SetSensitiveDetector(prescorer);
   

   return physWorld;
}
