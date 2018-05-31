#ifndef G4muDarkBremsstrahlungModel_h
#define G4muDarkBremsstrahlungModel_h

#include "G4VEmModel.hh"

struct ParamsForChi {double AA; double ZZ; double MMA; double EE0;};

class G4Element;
class G4ParticleChangeForLoss;

class G4muDarkBremsstrahlungModel : public G4VEmModel
{
   public:

      G4muDarkBremsstrahlungModel(const G4ParticleDefinition* p = 0,
                                 const G4String& nam = "muDBrem");

      virtual ~G4muDarkBremsstrahlungModel();

      virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

//      virtual G4double ComputeDEDXPerVolume(const G4Material*,
//                                            const G4ParticleDefinition*,
//                                           G4double kineticEnergy,
//                                            G4double cutEnergy);
      
     virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                 G4double tkin, 
                                                 G4double Z,   G4double,
                                                 G4double cut,
                                                 G4double maxE = DBL_MAX);
      
//     virtual G4double CrossSectionPerVolume(const G4Material*,
//                                            const G4ParticleDefinition*,
//                                            G4double kineticEnergy,
//                                            G4double cutEnergy,
//                                            G4double maxEnergy);
     
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

//      G4double ComputeBremLoss(G4double Z, G4double tkin, G4double cut);
      
//      G4double PositronCorrFactorLoss(G4double Z, G4double tkin, G4double cut);
      
//      G4double PositronCorrFactorSigma(G4double Z, G4double tkin, G4double cut);
      
      G4DataVector* ComputePartialSumSigma(const G4Material* material,
      G4double tkin, G4double cut);
      
//      G4double SupressionFunction(const G4Material* material, G4double tkin,
//      G4double gammaEnergy);
      
//      inline G4double ScreenFunction1(G4double ScreenVariable);
      
//      inline G4double ScreenFunction2(G4double ScreenVariable);
      
      // hide assignment operator
      G4muDarkBremsstrahlungModel & operator=(const  G4muDarkBremsstrahlungModel &right);
      G4muDarkBremsstrahlungModel(const  G4muDarkBremsstrahlungModel&);

   protected:

      const G4ParticleDefinition*   particle;
      G4ParticleDefinition*         theDarkPhoton;
      G4ParticleChangeForLoss*      fParticleChange;
      G4double                      MA;
      G4double                      MuonMass;
      G4bool                        isMuon;

   private:

      G4double highKinEnergy;
      G4double lowKinEnergy;
      G4double probsup;
//      G4double MigdalConstant;
//      G4double LPMconstant;
      G4bool   isInitialised;

   std::vector<G4DataVector*> partialSumSigma;

	};

#endif







