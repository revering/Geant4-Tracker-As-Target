#ifndef TTargSteppingAction_h
#define TTargSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class TTargEventAction;
class G4LogicalVolume;

class TTargSteppingAction : public G4UserSteppingAction
{
   public:
      TTargSteppingAction(TTargEventAction* eventAction);
      virtual ~TTargSteppingAction();

      virtual void UserSteppingAction(const G4Step*);

   private:
      TTargEventAction* fEventAction;
      G4LogicalVolume* fScoringVolume;
      int NEmissions;
};

#endif
