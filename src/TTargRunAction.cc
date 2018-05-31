#include "TTargRunAction.hh"
#include "TTargPrimaryGeneratorAction.hh"
#include "TTargDetectorConstruction.hh"
#include "TTargRun.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4ParameterManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "TTargAnalysis.hh"

TTargRunAction::TTargRunAction()
 : G4UserRunAction()
{
   G4AnalysisManager* anMan = G4AnalysisManager::Instance();
   anMan->SetFileName("TestOutput");

   anMan->CreateH1("Regular","Regular", 2, 0., 2);
   anMan->CreateH1("Dark","Dark", 2, 0., 2);

}

TTargRunAction::~TTargRunAction()
{
   delete G4AnalysisManager::Instance();
}

G4Run* TTargRunAction::GenerateRun()
{ return new TTargRun; }

void TTargRunAction::BeginOfRunAction(const G4Run*)
{
   G4AnalysisManager* anMan = G4AnalysisManager::Instance();
   anMan->OpenFile();
}

void TTargRunAction::EndOfRunAction(const G4Run* run)
{
   G4AnalysisManager* anMan = G4AnalysisManager::Instance();
   anMan->Write();
   anMan->CloseFile();
   
}


