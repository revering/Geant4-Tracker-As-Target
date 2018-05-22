#include "G4DMActionInitialization.hh"
#include "G4DMPrimaryGeneratorAction.hh"
#include "G4DMRunAction.hh"
#include "G4DMEventAction.hh"
#include "G4DMSteppingAction.hh"

G4DMActionInitialization::G4DMActionInitialization()
 : G4VUserActionInitialization()
 {}


G4DMActionInitialization::~G4DMActionInitialization()
{}

void G4DMActionInitialization::Build() const
{
   SetUserAction(new G4DMPrimaryGeneratorAction);
   
   G4DMRunAction* runAction = new G4DMRunAction;
   SetUserAction(runAction);
 
   G4DMEventAction* eventAction = new G4DMEventAction(runAction);
   SetUserAction(eventAction);

   G4DMSteppingAction* stepAction = new G4DMSteppingAction(eventAction);
   SetUserAction(stepAction);

}




