#include "G4DMRunAction.hh"
#include "G4DMPrimaryGeneratorAction.hh"
#include "G4DMDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParameterManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


G4DMRunAction::G4DMRunAction()
 : G4UserRunAction()
{
}

G4DMRunAction::~G4DMRunAction()
{
}

void G4DMRunAction::BeginOfRunAction(const G4Run*)
{

}

void G4DMRunAction::EndOfRunAction(const G4Run* run)
{
}

