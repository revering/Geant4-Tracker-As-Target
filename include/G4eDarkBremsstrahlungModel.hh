#ifndef G4eDarkBremsstrahlungModel_h
#define G4eDarkBremsstrahlungModel_h

#include "G4VEmModel.hh"

struct ParamsForChi {double AA; double ZZ; double MMA; double EE0;};

class G4Element;
class G4ParticleChangeForLoss;

class G4eDarkBremsstrahlungModel : public G4VEmModel
{
   public:

      G4eDarkBremsstrahlungModel(const G4ParticleDefinition* p = 0,
                                 const G4String& nam = "eDBrem");

      virtual ~G4eDarkBremsstrahlungModel();

      virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

      
     virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                 G4double tkin, 
                                                 G4double Z,   G4double,
                                                 G4double cut,
                                                 G4double maxE = DBL_MAX);


     G4DataVector* ComputePartialSumSigma(const G4Material* material,
                                          G4double tkin, G4double cut);

     virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                    const G4MaterialCutsCouple*,
                                    const G4DynamicParticle*,
                                    G4double tmin,
                                    G4double maxEnergy);

   protected:

      const G4Element* SelectRandomAtom(const G4MaterialCutsCouple* couple);

   private:

      void SetParticle(const G4ParticleDefinition* p);
     
      static G4double chi(double t, void * pp);
 
      static G4double DsigmaDx(double x, void * pp);
      
      // hide assignment operator
      G4eDarkBremsstrahlungModel & operator=(const  G4eDarkBremsstrahlungModel &right);
      G4eDarkBremsstrahlungModel(const  G4eDarkBremsstrahlungModel&);

   protected:

      const G4ParticleDefinition*   particle;
      G4ParticleDefinition*         theDarkPhoton;
      G4ParticleChangeForLoss*      fParticleChange;
      G4double                      MA;
      G4bool                        isElectron;

   private:

      G4double highKinEnergy;
      G4double lowKinEnergy;
      G4double probsup;
      G4bool   isInitialised;

      std::vector<G4DataVector*> partialSumSigma;
  
	};


#endif







