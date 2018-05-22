#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Parameter.hh"
#include "globals.hh"

class G4Run;

class G4DMRunAction : public G4UserRunAction
{
   public:
      G4DMRunAction();
      virtual ~G4DMRunAction();

      virtual void BeginOfRunAction(const G4Run*);
      virtual void EndOfRunAction(const G4Run*);
      
};

#endif

