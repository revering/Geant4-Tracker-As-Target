#ifndef G4DMActionInitialization_h
#define G4DMActionInitialization_h 

#include "G4DMPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"


class G4DMActionInitialization : public G4VUserActionInitialization
{
   public:
     G4DMActionInitialization();
     virtual ~G4DMActionInitialization();

     virtual void Build() const;

};

#endif
