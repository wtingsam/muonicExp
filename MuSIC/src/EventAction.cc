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
/// \file electromagnetic/TestEm2/src/EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 98761 2016-08-09 14:07:11Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EventManager.hh"
#include "ScintHit.hh"
#include "MWDCHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
using std::array;
using std::vector;

namespace {
    // Utility function which finds a hit collection with the given Id
    // and print warnings if not found
    G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
        auto hce = event->GetHCofThisEvent();
        if (!hce) {
            G4ExceptionDescription msg;
            msg << "No hits collection of this event found." << G4endl;
            G4Exception("EventAction::EndOfEventAction()",
                        "Code001", JustWarning, msg);
            return nullptr;
        }
        
        auto hc = hce->GetHC(collId);
        if ( ! hc) {
            G4ExceptionDescription msg;
            msg << "Hits collection " << collId << " of this event not found." << G4endl;
            G4Exception("EventAction::EndOfEventAction()",
                        "Code001", JustWarning, msg);
        }
        return hc;
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction(),
    ScintHCID  {{ -1, -1 }},
    DriftHCID{{ -1, -1, -1, -1, -1, -1, -1, -1 }}

{
    G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
 //additional initializations 
 Run* run 
   = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 run->InitializePerEvent();
    
    if (ScintHCID[0] == -1) {
        auto sdManager = G4SDManager::GetSDMpointer();
        
        // hits collections names
        array<G4String, kDim> hHCName
        = {{ "Degrader/ScintColl", "Trigger/ScintColl" }};
        array<G4String, kDim2> dHCName
        = {{ "MWDC1/MWDCColl", "MWDC2/MWDCColl", "MWDC3/MWDCColl", "MWDC4/MWDCColl", "MWDC5/MWDCColl", "MWDC6/MWDCColl", "MWDC7/MWDCColl", "MWDC8/MWDCColl"}};

        for (G4int iDet = 0; iDet < kDim; ++iDet) {
            // hit collections IDs
            ScintHCID[iDet]   = sdManager->GetCollectionID(hHCName[iDet]);
        }
        for (G4int iDet = 0; iDet < kDim2; ++iDet) {
        DriftHCID[iDet] = sdManager->GetCollectionID(dHCName[iDet]);
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* event)
{  
 Run* run 
   = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 run->FillPerEvent();
    
    auto fAnalysisManager = G4AnalysisManager::Instance();
    
    //Scintillator hits
    for (G4int iDet = 0; iDet < kDim; ++iDet) {
        auto hc = GetHC(event, ScintHCID[iDet]);
        if ( ! hc ) return;
        
        for (unsigned int i = 0; i<hc->GetSize(); ++i) {
            auto hit = static_cast<ScintHit*>(hc->GetHit(i));
            auto LocalPos = hit->GetLocalPos();
            
            fAnalysisManager->FillNtupleDColumn(iDet*5, hit->GetTime());
            fAnalysisManager->FillNtupleDColumn(iDet*5+1, LocalPos.x());
            fAnalysisManager->FillNtupleDColumn(iDet*5+2, LocalPos.y());
            fAnalysisManager->FillNtupleDColumn(iDet*5+3, LocalPos.z());
            fAnalysisManager->FillNtupleDColumn(iDet*5+4, hit->GetEdep());
        }
    }
    
    // Drift chambers hits
    for (G4int iDet = 0; iDet < kDim2; ++iDet) {
        auto hc = GetHC(event, DriftHCID[iDet]);
        if ( ! hc ) return;
        auto nhit = hc->GetSize();
        for (unsigned long i = 0; i < nhit; ++i) {
            auto hit = static_cast<MWDCHit*>(hc->GetHit(i));
            auto LocalPos = hit->GetLocalPos();
            
            fAnalysisManager->FillNtupleDColumn(iDet*5+10, hit->GetTime());
            fAnalysisManager->FillNtupleDColumn(iDet*5+10+1, LocalPos.x());
            fAnalysisManager->FillNtupleDColumn(iDet*5+10+2, LocalPos.y());
            fAnalysisManager->FillNtupleDColumn(iDet*5+10+3, LocalPos.z());
            fAnalysisManager->FillNtupleDColumn(iDet*5+10+4, hit->GetEdep());
        }
    }
    
    fAnalysisManager->AddNtupleRow();
    /*
    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( printModulo == 0 || event->GetEventID() % printModulo != 0) return;

    
    auto primary = event->GetPrimaryVertex(0)->GetPrimary(0);
    G4cout
    << G4endl
    << ">>> Event " << event->GetEventID() << " >>> Simulation truth : "
    << primary->GetG4code()->GetParticleName()
    << " " << primary->GetMomentum() << G4endl;
    
    
    // Hodoscopes
    
    for (G4int iDet = 0; iDet < kDim; ++iDet) {
        auto hc = GetHC(event, ScintHCID[iDet]);
        if ( ! hc ) return;
        G4cout << "Hodoscope " << iDet + 1 << " has " << hc->GetSize()  << " hits." << G4endl;
        for (unsigned int i = 0; i<hc->GetSize(); ++i) {
            hc->GetHit(i)->Print();
        }
    }
    */
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

