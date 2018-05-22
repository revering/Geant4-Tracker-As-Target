#ifndef TTargEventAction_h
#define TTargEventAction_h 1

#include "G4UserEventAction.hh"
#include "TTargRunAction.hh"
#include "globals.hh"


class TTargEventAction : public G4UserEventAction
{
   public:
      TTargEventAction(TTargRunAction* runAction);
      virtual ~TTargEventAction();
 
      virtual void BeginOfEventAction(const G4Event* event);
      virtual void EndOfEventAction(const G4Event* event);

      void AddEdep(G4double edep) { fEdep += edep; }
      TTargRunAction* GetRunAction() const;
      G4int GetnBrem() { return nBrem; }
      void incnBrem() { nBrem++;}
      void rstnBrem() { nBrem = 0;}

   private:
      TTargRunAction* fRunAction;
      G4double fEdep;
      G4int nBrem;
};

#endif

