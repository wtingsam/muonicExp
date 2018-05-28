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
// $Id: MWDCSD.cc 101036 2016-11-04 09:00:23Z gcosmo $
//
/// \file MWDCSD.cc
/// \brief Implementation of the MWDCSD class

#include "MWDCSD.hh"
#include "MWDCHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MWDCSD::MWDCSD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert( "MWDCColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MWDCSD::~MWDCSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MWDCSD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new MWDCHitsCollection
  (SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) { 
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool MWDCSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
  auto edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;
  
  auto preStepPoint = step->GetPreStepPoint();
  auto touchable = preStepPoint->GetTouchable();
  auto copyNo = touchable->GetVolume()->GetCopyNo();
  auto hitTime = preStepPoint->GetGlobalTime();
  
    auto worldPos = preStepPoint->GetPosition();
    auto localPos
    = touchable->GetHistory()->GetTopTransform().TransformPoint(worldPos);
    
  // check if this finger already has a hit
  auto ix = -1;
  for (auto i=0;i<fHitsCollection->entries();i++) {
    if ((*fHitsCollection)[i]->GetID()==copyNo) {
      ix = i;
      break;
    }
  }

  if (ix>=0) {
    // if it has, then take the earlier time
    if ((*fHitsCollection)[ix]->GetTime()>hitTime) { 
      (*fHitsCollection)[ix]->SetTime(hitTime); 
    }
  }
  else {
    // if not, create a new hit and set it to the collection
    auto hit = new MWDCHit(copyNo,hitTime);
    auto physical = touchable->GetVolume();
    hit->SetLogV(physical->GetLogicalVolume());
    auto transform = touchable->GetHistory()->GetTopTransform();
    transform.Invert();
    hit->SetRot(transform.NetRotation());
    hit->SetPos(transform.NetTranslation());
      
      hit->SetWorldPos(worldPos);
      hit->SetLocalPos(localPos);
      
      hit->SetEdep(edep);

    fHitsCollection->insert(hit);
  }    
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
