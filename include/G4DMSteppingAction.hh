#ifndef G4DMSteppingAction_h
#define G4DMSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4DMEventAction;
class G4LogicalVolume;

class G4DMSteppingAction : public G4UserSteppingAction
{
   public:
      G4DMSteppingAction(G4DMEventAction* eventAction);
      virtual ~G4DMSteppingAction();

      virtual void UserSteppingAction(const G4Step*);

   private:
      G4DMEventAction* fEventAction;
      G4LogicalVolume* fScoringVolume;
      int NEmissions;
};

#endif
