#ifndef TTargActionInitialization_h
#define TTargActionInitialization_h 

#include "TTargPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"


class TTargActionInitialization : public G4VUserActionInitialization
{
   public:
     TTargActionInitialization();
     virtual ~TTargActionInitialization();

     virtual void Build() const;

};

#endif
