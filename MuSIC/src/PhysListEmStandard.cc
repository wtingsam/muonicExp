//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: PhysListEmStandard.cc 100289 2016-10-17 08:46:50Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmStandard.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4KleinNishinaModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCapture.hh"
#include "G4MuonMinusAtomicCapture.hh"
#include "G4MuMinusCapturePrecompound.hh"
#include "G4MuonicAtomDecay.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4MscStepLimitType.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"

#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4ProcessTable.hh"
#include "G4DecayWithSpin.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetMinEnergy(10*eV);
    param->SetMaxEnergy(10*TeV);
    param->SetNumberOfBinsPerDecade(10);
    param->SetMscStepLimitType(fUseSafetyPlus);
    SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
    G4cout << "test1" << G4endl;

  // Add standard EM Processes
  //
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      ////ph->RegisterProcess(new G4RayleighScattering, particle);      
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);      
      G4ComptonScattering* cs   = new G4ComptonScattering;
      cs->SetEmModel(new G4KleinNishinaModel());
      ph->RegisterProcess(cs, particle);
      ph->RegisterProcess(new G4GammaConversion, particle);
     
    } else if (particleName == "e-") {
    
      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      //            
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);      
      ph->RegisterProcess(eIoni, particle);
      //
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);      
            
    } else if (particleName == "e+") {
    
      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      //     
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);      
      ph->RegisterProcess(eIoni, particle);
      //
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      //
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);    
                  
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {

      ph->RegisterProcess(new G4MuMultipleScattering(), particle);
      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.1, 50*um);      
      ph->RegisterProcess(muIoni, particle);
      ph->RegisterProcess(new G4MuBremsstrahlung(), particle);
      ph->RegisterProcess(new G4MuPairProduction(), particle);
        
        if( particleName == "mu-" ){
            //ph->RegisterProcess(new G4MuonMinusCapture(),particle);
            //ph->RegisterProcess(new G4MuonMinusAtomicCapture(),particle);
            //ph->RegisterProcess(new G4MuonMinusCapture(new G4MuMinusCapturePrecompound()),particle);
            ph->RegisterProcess(new G4MuonicAtomDecay(),particle);
        }
      
    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);      
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.1, 20*um);
      ph->RegisterProcess(hIoni, particle);
      ph->RegisterProcess(new G4hBremsstrahlung(), particle);
      ph->RegisterProcess(new G4hPairProduction(), particle);            
     
    } else if( particleName == "alpha" || 
               particleName == "He3"    ) {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);           
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 1*um);
      ph->RegisterProcess(ionIoni, particle);
      ph->RegisterProcess(new G4NuclearStopping(), particle);      
            
    } else if( particleName == "GenericIon" ) {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);          
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 1*um);
      ph->RegisterProcess(ionIoni, particle);
      ph->RegisterProcess(new G4NuclearStopping(), particle);                   
      
    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
               
      //all others charged particles except geantino
      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4hIonisation(), particle);
    }
  }

    //
    /*
    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();
    G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();
    G4VProcess* decay;
    decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());
    G4ProcessManager* pmanager;
    pmanager = G4MuonMinus::MuonMinus()->GetProcessManager();
    
    if (pmanager) {
        if (decay) pmanager->RemoveProcess(decay);
        pmanager->AddProcess(decayWithSpin);
        // set ordering for PostStepDoIt and AtRestDoIt
        pmanager ->SetProcessOrdering(decayWithSpin, idxPostStep);
        pmanager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }
    
    G4VProcess* process = processTable->
    FindProcess("muMinusCaptureAtRest",G4MuonMinus::MuonMinus());
    
    if (pmanager) {
        if (process) pmanager->RemoveProcess(process);
        process = new G4MuonMinusCapture(new G4MuMinusCapturePrecompound());
        pmanager->AddRestProcess(process);
    }
     */
  //
    
  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  de->SetFluo(true);
  de->SetAuger(false);   
  de->SetPIXE(false);  
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

