#include "TTargDetectorConstruction.hh"
#include "G4Region.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

DMDetectorConstruction::DMDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{  }

DMDetectorConstruction::~DMDetectorConstruction()
{  }

G4VPhysicalVolume* DMDetectorConstruction::Construct()
{
   G4NistManager* nist = G4NistManager::Instance();

   G4double env_sizeXY = 100*cm, env_sizeZ = 100*cm;
   G4Material* env_mat = nist->FindOrBuildMaterial("G4_W");

   G4bool checkOverlaps = true;

   G4double world_sizeXY = 1.2*env_sizeXY;
   G4double world_sizeZ = 1.2*env_sizeZ;
   G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");


   G4Box* SolidWorld = new G4Box("World", 1*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);

   G4LogicalVolume* LogicWorld = new G4LogicalVolume(SolidWorld, world_mat, "World");

   G4VPhysicalVolume* physWorld = new G4PVPlacement( 0, G4ThreeVector(), LogicWorld, "World", 0, false, 0, checkOverlaps);

   G4Box* detector = new G4Box("Detector", 1*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ);

   G4LogicalVolume* ldetector = new G4LogicalVolume(detector, world_mat, "Detector");

   new G4PVPlacement(0,G4ThreeVector(), ldetector, "Detector", LogicWorld, false, 0, checkOverlaps);

   G4Box* solidEnv = new G4Box("Target", 0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ);

   G4LogicalVolume* logicEnv = new G4LogicalVolume(solidEnv, env_mat, "Target");

   new G4PVPlacement(0,G4ThreeVector(50*cm,0,0), logicEnv, "Target", ldetector, false, 0, checkOverlaps);
    
   const G4String& rName = "Tracker";

   G4Region* tracker = new G4Region(rName);
   logicEnv->SetRegion(tracker);
   tracker->AddRootLogicalVolume(logicEnv);

   G4LogicalVolume* box2 = new G4LogicalVolume(solidEnv, env_mat, "Box2");

   new G4PVPlacement(0, G4ThreeVector(-50*cm,0,0), box2, "Box2", ldetector, false, 0, checkOverlaps);

   const G4String& r2Name = "Box2";

   G4Region* ltracker = new G4Region(r2Name);
   box2->SetRegion(ltracker);
   ltracker->AddRootLogicalVolume(box2);

   fScoringVolume = logicEnv;

   return physWorld;
}
