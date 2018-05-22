#include "TTargActionInitialization.hh"
#include "TTargPrimaryGeneratorAction.hh"
#include "TTargRunAction.hh"
#include "TTargEventAction.hh"
#include "TTargSteppingAction.hh"

TTargActionInitialization::TTargActionInitialization()
 : G4VUserActionInitialization()
 {}


TTargActionInitialization::~TTargActionInitialization()
{}

void TTargActionInitialization::Build() const
{
   SetUserAction(new TTargPrimaryGeneratorAction);
   
   TTargRunAction* runAction = new TTargRunAction;
   SetUserAction(runAction);
 
   TTargEventAction* eventAction = new TTargEventAction(runAction);
   SetUserAction(eventAction);

   TTargSteppingAction* stepAction = new TTargSteppingAction(eventAction);
   SetUserAction(stepAction);

}




