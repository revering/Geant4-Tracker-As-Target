#include "TTargRunAction.hh"
#include "TTargPrimaryGeneratorAction.hh"
#include "TTargDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParameterManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


TTargRunAction::TTargRunAction()
 : G4UserRunAction()
{
}

TTargRunAction::~TTargRunAction()
{
}

void TTargRunAction::BeginOfRunAction(const G4Run*)
{

}

void TTargRunAction::EndOfRunAction(const G4Run* run)
{
}

