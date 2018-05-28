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
/// \file electromagnetic/TestEm2/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// 
// $Id: DetectorConstruction.cc 98761 2016-08-09 14:07:11Z gcosmo $

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4GlobalMagFieldMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4AutoDelete.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "ScintSD.hh"
#include "MWDCSD.hh"
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fNLtot(40),fNRtot(50),fDLradl(0.5),fDRradl(0.1),
 fDLlength(0.),fDRlength(0.),
 fMaterial(0), LogicDeg(nullptr),
 fEcalLength(0.),fEcalRadius(0.),
 fSolidEcal(0),fLogicEcal(0),fPhysiEcal(0),
 fDetectorMessenger(0)
{
  DefineMaterials();
  SetMaterial("HydrogenGas");
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    auto nistManager = G4NistManager::Instance();

  //
  // define few Elements by hand
  //
  G4double a, z;
    
  G4Element* H  = new G4Element("Hydrogen",  "H", z= 1., a=   1.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,  "O", z= 8., a=  16.00*g/mole);
  G4Element* Ge = new G4Element("Germanium", "Ge",z=32., a=  72.59*g/mole);
  G4Element* Bi = new G4Element("Bismuth",   "Bi",z=83., a= 208.98*g/mole);

  //
  // define materials
  //
  G4double density;
  G4int ncomponents, natoms;

  // water with ionisation potential 78 eV
  G4Material* H2O = 
  new G4Material("Water", density= 1.00*g/cm3, ncomponents=2);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);

  // pure materails
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
  new G4Material("Aluminium",   z=13., a= 26.98*g/mole, density= 2.7*g/cm3);
  new G4Material("Iron",        z=26., a= 55.85*g/mole, density= 7.87*g/cm3);  
  new G4Material("Copper",      z=29., a= 63.55*g/mole, density= 8.960*g/cm3); 
  new G4Material("Tungsten",    z=74., a=183.84*g/mole, density=19.35*g/cm3); 
  new G4Material("Lead",        z=82., a=207.19*g/mole, density=11.35*g/cm3);  
  new G4Material("Uranium"    , z=92., a=238.03*g/mole, density= 18.95*g/cm3);

  // compound material
  G4Material* BGO = 
  new G4Material("BGO", density= 7.10*g/cm3, ncomponents=3);
  BGO->AddElement(O , natoms=12);
  BGO->AddElement(Ge, natoms= 3);
  BGO->AddElement(Bi, natoms= 4);
    
    new G4Material("HydrogenGas", z=1, a=2.016*g/mole, density= 0.8988*mg/cm3,
                   kStateGas, 300.0*kelvin, 1*atmosphere);

    new G4Material("ArgonGas", z=18, a=39.948*g/mole, density= 1.782*mg/cm3,
                   kStateGas, 300.0*kelvin, 1*atmosphere);
    
    nistManager->FindOrBuildMaterial("G4_AIR");
    nistManager->FindOrBuildMaterial("G4_Ag");
    nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                   kStateGas, 2.73*kelvin, 3.e-18*pascal);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  //G4double Radl = fMaterial->GetRadlen();
    G4double Radl = 0.1*m;
    fDLlength = fDLradl*Radl; fDRlength = fDRradl*Radl;
    fEcalLength = fNLtot*fDLlength;  fEcalRadius = fNRtot*fDRlength;

  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

    DefineMaterials();
    auto WorldMaterial = G4Material::GetMaterial("G4_AIR");
    auto MWDCMaterial = G4Material::GetMaterial("ArgonGas");
    auto DegraderMaterial = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

    auto ChamberMaterial = G4Material::GetMaterial("G4_Ag");
    auto GasMaterial = G4Material::GetMaterial("HydrogenGas");

    
    //World
    auto SolidWorld = new G4Box("World",2*m/2,2*m/2,2*m/2);
    auto LogicWorld = new G4LogicalVolume(SolidWorld,WorldMaterial,"World");
    auto PhysiWorld = new G4PVPlacement(0,G4ThreeVector(),LogicWorld,"World",0,false,0);
    
    //Beampipe
    auto SolidDuct = new G4Tubs("Duct",60.0*mm/2.0,65.0*mm/2.0,100*mm/2.0,0,360*deg);
    auto LogicDuct = new G4LogicalVolume(SolidDuct,ChamberMaterial,"Duct",0,0,0);
    auto PhysiDuct = new G4PVPlacement(0,G4ThreeVector(0,0,-570*mm),LogicDuct,"Duct",LogicWorld,false,0);
    
    //Degrador
    auto SolidDeg = new G4Box("Degrader",100*mm/2,100*mm/2,5*mm/2);
     LogicDeg = new G4LogicalVolume(SolidDeg,DegraderMaterial,"Degrader");
    auto PhysiDeg = new G4PVPlacement(0,G4ThreeVector(0,0,-500*mm),LogicDeg,"Degrader",LogicWorld,false,0);

    //Trigger
    auto SolidTrg = new G4Box("Trigger",240*mm/2,5*mm/2,200*mm/2);
     LogicTrg = new G4LogicalVolume(SolidTrg,DegraderMaterial,"Trigger");
    auto PhysiTrg = new G4PVPlacement(0,G4ThreeVector(0,500*mm,0),LogicTrg,"Trigger",LogicWorld,false,0);
    
    //Bottom Counter
    //Trigger
    auto SolidTrg2 = new G4Box("Trigger2",240*mm/2,5*mm/2,200*mm/2);
    auto LogicTrg2 = new G4LogicalVolume(SolidTrg2,DegraderMaterial,"Trigger2");
    auto PhysiTrg2 = new G4PVPlacement(0,G4ThreeVector(0,-500*mm,0),LogicTrg2,"Trigger2",LogicWorld,false,0);
    
    //Plate on the magnet
    auto SolidPlate = new G4Box("Plate",100*mm/2,10*mm/2,450*mm/2);
    auto LogicPlate = new G4LogicalVolume(SolidPlate,DegraderMaterial,"Plate");
    auto PhysiPlate = new G4PVPlacement(0,G4ThreeVector(110*mm,450*mm,0),LogicPlate,"Plate",LogicWorld,false,0);
    PhysiPlate = new G4PVPlacement(0,G4ThreeVector(-110*mm,450*mm,0),LogicPlate,"Plate",LogicWorld,false,0);
    
    //MWDC for beam tracking
    auto SolidMWDC1 = new G4Box("MWDC",200*mm/2,240*mm/2,5*mm/2);
     LogicMWDC1 = new G4LogicalVolume(SolidMWDC1,MWDCMaterial,"MWDC");
    auto PhysiMWDC1 = new G4PVPlacement(0,G4ThreeVector(0,0,-450*mm),LogicMWDC1,"MWDC",LogicWorld,false,0);
    
    auto SolidMWDC2 = new G4Box("MWDC",200*mm/2,240*mm/2,5*mm/2);
     LogicMWDC2 = new G4LogicalVolume(SolidMWDC2,MWDCMaterial,"MWDC");
    auto PhysiMWDC2 = new G4PVPlacement(0,G4ThreeVector(0,0,-460*mm),LogicMWDC2,"MWDC",LogicWorld,false,0);
    
    auto SolidMWDC3 = new G4Box("MWDC",200*mm/2,240*mm/2,5*mm/2);
     LogicMWDC3 = new G4LogicalVolume(SolidMWDC3,MWDCMaterial,"MWDC");
    auto PhysiMWDC3 = new G4PVPlacement(0,G4ThreeVector(0,0,-470*mm),LogicMWDC3,"MWDC",LogicWorld,false,0);
    
    auto SolidMWDC4 = new G4Box("MWDC",200*mm/2,240*mm/2,5*mm/2);
     LogicMWDC4 = new G4LogicalVolume(SolidMWDC4,MWDCMaterial,"MWDC");
    auto PhysiMWDC4 = new G4PVPlacement(0,G4ThreeVector(0,0,-480*mm),LogicMWDC4,"MWDC",LogicWorld,false,0);
    
    //MWDC for electron tracking
    auto SolidMWDC5 = new G4Box("MWDC",240*mm/2,5*mm/2,200*mm/2);
     LogicMWDC5 = new G4LogicalVolume(SolidMWDC5,MWDCMaterial,"MWDC");
    auto PhysiMWDC5 = new G4PVPlacement(0,G4ThreeVector(0,460*mm,0),LogicMWDC5,"MWDC",LogicWorld,false,0);
    
    auto SolidMWDC6 = new G4Box("MWDC",240*mm/2,5*mm/2,200*mm/2);
     LogicMWDC6 = new G4LogicalVolume(SolidMWDC6,MWDCMaterial,"MWDC");
    auto PhysiMWDC6 = new G4PVPlacement(0,G4ThreeVector(0,470*mm,0),LogicMWDC6,"MWDC",LogicWorld,false,0);
    
    auto SolidMWDC7 = new G4Box("MWDC",240*mm/2,5*mm/2,200*mm/2);
     LogicMWDC7 = new G4LogicalVolume(SolidMWDC7,MWDCMaterial,"MWDC");
    auto PhysiMWDC7 = new G4PVPlacement(0,G4ThreeVector(0,480*mm,0),LogicMWDC7,"MWDC",LogicWorld,false,0);
    
    auto SolidMWDC8 = new G4Box("MWDC",240*mm/2,5*mm/2,200*mm/2);
     LogicMWDC8 = new G4LogicalVolume(SolidMWDC8,MWDCMaterial,"MWDC");
    auto PhysiMWDC8 = new G4PVPlacement(0,G4ThreeVector(0,490*mm,0),LogicMWDC8,"MWDC",LogicWorld,false,0);
    
    //
    //Chamber
    auto SolidChamber = new G4Tubs("Chamber",420.0*mm/2.0,425.0*mm/2.0,360*mm/2.0,0,340*deg);
    auto LogicChamber = new G4LogicalVolume(SolidChamber,ChamberMaterial,"Chamber",0,0,0);
    
    G4RotationMatrix *RotChamber = new G4RotationMatrix();  // make it unit vector
    RotChamber->rotateY(90*deg);
    RotChamber->rotateZ(170*deg);

    auto PhysiChamber = new G4PVPlacement(RotChamber,G4ThreeVector(),LogicChamber,"Chamber",LogicWorld,false,0);
    
    //Gas
    auto SolidGas = new G4Tubs("Gas",0.,420.0*mm/2.0,360*mm/2.0,0.,360*deg);
    auto LogicGas = new G4LogicalVolume(SolidGas,GasMaterial,"Gas",0,0,0);
    
    G4RotationMatrix *RotGas = new G4RotationMatrix();  // make it unit vector
    RotGas->rotateY(90*deg);
    
    auto PhysiGas = new G4PVPlacement(RotGas,G4ThreeVector(),LogicGas,"Gas",LogicWorld,false,0);

    //Chamber side
    auto SolidChamberSide = new G4Tubs("Chamber",0.,425.0*mm/2.0,10*mm/2.0,0,360*deg);
    auto LogicChamberSide = new G4LogicalVolume(SolidChamberSide,ChamberMaterial,"Chamber",0,0,0);
    auto PhysiChamberSide = new G4PVPlacement(RotGas,G4ThreeVector(185*mm,0,0),LogicChamberSide,"Chamber",LogicWorld,false,0);
    PhysiChamberSide = new G4PVPlacement(RotGas,G4ThreeVector(-185*mm,0,0),LogicChamberSide,"Chamber",LogicWorld,false,0);

    //Coils
    auto SolidCoil = new G4Tubs("Coil",440.0*mm/2.0,880.0*mm/2.0,100*mm/2.0,0,360*deg);
    auto LogicCoil = new G4LogicalVolume(SolidCoil,ChamberMaterial,"Coil",0,0,0);
    auto PhysiCoil = new G4PVPlacement(RotGas,G4ThreeVector(110*mm,0,0),LogicCoil,"Coil",LogicWorld,false,0);
    PhysiCoil = new G4PVPlacement(RotGas,G4ThreeVector(-110*mm,0,0),LogicCoil,"Coil",LogicWorld,false,0);
    
    //Visualization
    auto VisWorld = new G4VisAttributes(G4Colour(0,0,0));
    VisWorld->SetVisibility(false);
    LogicWorld->SetVisAttributes(VisWorld);
    
    auto VisDegrader = new G4VisAttributes(G4Colour(255.0/255.0,255.0/255.0,204.0/255.0));
    LogicDeg->SetVisAttributes(VisDegrader);
    
    auto VisTrigger = new G4VisAttributes(G4Colour(255.0/255.0,255.0/255.0,204.0/255.0));
    LogicTrg->SetVisAttributes(VisDegrader);
    
    auto VisMWDC = new G4VisAttributes(G4Colour(0.8,0.2,0.1));
    LogicMWDC1->SetVisAttributes(VisMWDC);
    LogicMWDC2->SetVisAttributes(VisMWDC);
    LogicMWDC3->SetVisAttributes(VisMWDC);
    LogicMWDC4->SetVisAttributes(VisMWDC);
    LogicMWDC5->SetVisAttributes(VisMWDC);
    LogicMWDC6->SetVisAttributes(VisMWDC);
    LogicMWDC7->SetVisAttributes(VisMWDC);
    LogicMWDC8->SetVisAttributes(VisMWDC);

    auto VisChamber = new G4VisAttributes(G4Colour(128.0/255.0,128.0/255.0,128.0/255.0));
    LogicChamber->SetVisAttributes(VisChamber);
    
    auto VisChamberSide = new G4VisAttributes(G4Colour(128.0/255.0,128.0/255.0,128.0/255.0));
    LogicChamberSide->SetVisAttributes(VisChamberSide);
    
    auto VisGas = new G4VisAttributes(G4Colour(128.0/255.0,128.0/255.0,128.0/255.0));
    LogicGas->SetVisAttributes(VisGas);
    
    auto VisCoil = new G4VisAttributes(G4Colour(0.8,0.2,0.2));
    LogicCoil->SetVisAttributes(VisCoil);

  //
  // Ecal
  //
  fSolidEcal = new G4Tubs("Ecal",0.,120*mm/2.0,500*mm/2.0,0.,360*deg);
  fLogicEcal = new G4LogicalVolume( fSolidEcal,fMaterial,"Ecal",0,0,0);
  //fPhysiEcal = new G4PVPlacement(0,G4ThreeVector(),fLogicEcal,"Ecal",LogicWorld,false,0);

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4cout << "Absorber is " << G4BestUnit(fEcalLength,"Length")
         << " of " << fMaterial->GetName() << G4endl;
  //
  //always return the physical World
  //
  return PhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(const G4String& materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if(pttoMaterial &&  fMaterial != pttoMaterial) {
    fMaterial = pttoMaterial;
    if(fLogicEcal) { fLogicEcal->SetMaterial(fMaterial); }
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetLBining(G4ThreeVector Value)
{
  fNLtot = (G4int)Value(0);
  if (fNLtot > kMaxBin) {
    G4cout << "\n ---> warning from SetLBining: "
           << fNLtot << " truncated to " << kMaxBin << G4endl;
    fNLtot = kMaxBin;
  }  
  fDLradl = Value(1);
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetRBining(G4ThreeVector Value)
{
  fNRtot = (G4int)Value(0);
  if (fNRtot > kMaxBin) {
    G4cout << "\n ---> warning from SetRBining: "
           << fNRtot << " truncated to " << kMaxBin << G4endl;
    fNRtot = kMaxBin;
  }    
  fDRradl = Value(1);
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{

    auto sdManager = G4SDManager::GetSDMpointer();
    G4String SDname;
    
    auto degSD = new ScintSD(SDname="Degrader");
    sdManager->AddNewDetector(degSD);
    LogicDeg->SetSensitiveDetector(degSD);
  
    auto trgSD = new ScintSD(SDname="Trigger");
    sdManager->AddNewDetector(trgSD);
    LogicTrg->SetSensitiveDetector(trgSD);
    
    auto MWDC1SD = new MWDCSD(SDname="MWDC1");
    sdManager->AddNewDetector(MWDC1SD);
    LogicMWDC1->SetSensitiveDetector(MWDC1SD);
    
    auto MWDC2SD = new MWDCSD(SDname="MWDC2");
    sdManager->AddNewDetector(MWDC2SD);
    LogicMWDC2->SetSensitiveDetector(MWDC2SD);
    
    auto MWDC3SD = new MWDCSD(SDname="MWDC3");
    sdManager->AddNewDetector(MWDC3SD);
    LogicMWDC3->SetSensitiveDetector(MWDC3SD);
    
    auto MWDC4SD = new MWDCSD(SDname="MWDC4");
    sdManager->AddNewDetector(MWDC4SD);
    LogicMWDC4->SetSensitiveDetector(MWDC4SD);
    
    auto MWDC5SD = new MWDCSD(SDname="MWDC5");
    sdManager->AddNewDetector(MWDC5SD);
    LogicMWDC5->SetSensitiveDetector(MWDC5SD);
    
    auto MWDC6SD = new MWDCSD(SDname="MWDC6");
    sdManager->AddNewDetector(MWDC6SD);
    LogicMWDC6->SetSensitiveDetector(MWDC6SD);
    
    auto MWDC7SD = new MWDCSD(SDname="MWDC7");
    sdManager->AddNewDetector(MWDC7SD);
    LogicMWDC7->SetSensitiveDetector(MWDC7SD);
    
    auto MWDC8SD = new MWDCSD(SDname="MWDC8");
    sdManager->AddNewDetector(MWDC8SD);
    LogicMWDC8->SetSensitiveDetector(MWDC8SD);
    
    if ( fFieldMessenger.Get() == 0 ) {
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue = G4ThreeVector();
        G4GlobalMagFieldMessenger* msg =
            new G4GlobalMagFieldMessenger(fieldValue);
        //msg->SetVerboseLevel(1);
        G4AutoDelete::Register(msg);
        fFieldMessenger.Put( msg );

    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
