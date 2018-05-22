#include "G4UImanager.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "TTargPhysicsList.hh"
#include "TTargDetectorConstruction.hh"
#include "TTargActionInitialization.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

int main( int argc, char** argv )
{
#ifdef G4MULTITHREADED
   G4MTRunManager* runManager = new G4MTRunManager;
#else
   G4RunManager* runManager = new G4RunManager;  
#endif

   runManager->SetUserInitialization(new DMDetectorConstruction);
   runManager->SetUserInitialization(new DMPhysicsList);
   runManager->SetUserInitialization(new TTargActionInitialization());


   runManager->Initialize();

   G4UImanager* UI = G4UImanager::GetUIpointer();	
   G4VisManager* visManager = new G4VisExecutive;
   visManager->Initialize();

   UI->ApplyCommand("/run/verbose 1");
   UI->ApplyCommand("/event/verbose 1");
   UI->ApplyCommand("/tracking/verbose 0");

   if( argc == 1)
   {
   #ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      UI->ApplyCommand("/control/execute init_vis.mac");
      ui->SessionStart();
      delete ui;
   #endif
   }
   else
   {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
   }

   

   delete visManager;
   delete runManager;
   return 0;
}
