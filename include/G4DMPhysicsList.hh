#ifndef G4DMPhysicsList
#define G4DMPhysicsList 

#include "G4VModularPhysicsList.hh"

class DMPhysicsList : public G4VModularPhysicsList
{
  public:
     DMPhysicsList();
     virtual ~DMPhysicsList();
     
     virtual void SetCuts();
     virtual void ConstructParticle();
     virtual void ConstructProcess();

//     void AddPhysicsList(const G4String& name);

   private:

      G4VPhysicsConstructor* fEmPhysicsList;
};

#endif
