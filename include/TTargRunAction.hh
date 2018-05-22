#ifndef B1RunAction_h
#define B1RunAction_h 1

#include "G4UserRunAction.hh"
#include "G4Parameter.hh"
#include "globals.hh"

class G4Run;

class TTargRunAction : public G4UserRunAction
{
   public:
      TTargRunAction();
      virtual ~TTargRunAction();

      virtual void BeginOfRunAction(const G4Run*);
      virtual void EndOfRunAction(const G4Run*);
      
};

#endif

