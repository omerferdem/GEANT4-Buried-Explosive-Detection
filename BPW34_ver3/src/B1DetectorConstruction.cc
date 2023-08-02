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
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 2*cm, env_sizeZ = 2*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.*env_sizeXY;
  G4double world_sizeZ  = 1.*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //     
  // Shape 2
  //

  G4String name, symbol;
  G4double iz, z, a, n, abundance, density, numberOfatoms, ncomponents, numisotopes;
  
  //G4Material* shell_mat = nist->FindOrBuildMaterial("G4_BORON_CARBIDE");

  //G4Element* Li = new G4Element(name = "Lithium", symbol = "Li", z = 3.00, a = 6.941 * g / mole);
  //G4Element* F = new G4Element(name = "Flourine", symbol = "F", z = 9.00, a = 18.998 * g / mole);
  //G4Material* LiF_mat = new G4Material(name = "LiF", density = 2.64 * g / cm3, ncomponents = 2);
  //LiF_mat->AddElement(Li, 1);
  //LiF_mat->AddElement(F, 1);

  //G4Isotope* B10 = new G4Isotope(name = "B10", iz = 4, n = 10);
  //G4Isotope* B11 = new G4Isotope(name = "B11", iz = 4, n = 11);
  //G4Element* Enr_Boron = new G4Element(name="B", symbol="B", numisotopes = 2);
  //Enr_Boron->AddIsotope(B10, abundance = 90. * perCent);
  //Enr_Boron->AddIsotope(B11, abundance = 10. * perCent);
  //G4Element* C = new G4Element(name = "Carbon", symbol = "C", z = 6.00, a = 12.011 * g / mole);
  //G4Material* Enr_B4C = new G4Material(name = "B4C", density = 2.52 * g / cm3, ncomponents = 5);
  //Enr_B4C->AddElement(Enr_Boron, numberOfatoms = 1);
  //Enr_B4C->AddElement(Enr_Boron, numberOfatoms = 1);
  //Enr_B4C->AddElement(Enr_Boron, numberOfatoms = 1);
  //Enr_B4C->AddElement(Enr_Boron, numberOfatoms = 1);
  //Enr_B4C->AddElement(C, numberOfatoms = 1);
  
  G4Isotope* Li6 = new G4Isotope(name = "Li6", iz = 3, n = 6);
  G4Isotope* Li7 = new G4Isotope(name = "Li7", iz = 3, n = 7);
  G4Element* Enr_Lithium = new G4Element(name = "Enriched_Lithium", symbol = "Li", numisotopes = 2);
  Enr_Lithium->AddIsotope(Li6, abundance = 90. * perCent);
  Enr_Lithium->AddIsotope(Li7, abundance = 10. * perCent);
  G4Element* F = new G4Element(name = "Flourine", symbol = "F", z = 9.00, a = 18.998 * g / mole);
  G4Material* Enr_LiF = new G4Material(name = "Enriched_LiF", density = 2.64 * g / cm3, ncomponents = 2);
  Enr_LiF->AddElement(Enr_Lithium, numberOfatoms = 1);
  Enr_LiF->AddElement(F, numberOfatoms = 1);

  G4double thickness_param = 0.014 * mm;

  G4double COUNT_X = 1;
  G4double COUNT_Y = 1;
  G4double COUNT_Z = 1;

  G4double sens_hx = 0.1 * mm;
  G4double sens_hy = 1.5 * mm;
  G4double sens_hz = 1.25 * mm;

  G4double shell_hx = thickness_param+sens_hx * mm;
  G4double shell_hy = thickness_param + 1.5 * mm;
  G4double shell_hz = thickness_param + 1.25 * mm;

  G4double pos_shell_x = 0.0 * mm;
  G4double pos_shell_y = 0.0 * mm;
  G4double pos_shell_z = 0.0 * mm;

  G4Box* shellBox = new G4Box("Shell", shell_hx, shell_hy, shell_hz);
  G4LogicalVolume* shellLog = new G4LogicalVolume(shellBox, Enr_LiF, "Shell");

  //   ***   Placing Shells of Photodiodes   ***

  for (G4int kshell = 0; kshell < COUNT_X; kshell++) {
      pos_shell_x = shell_hx + shell_hx * 2 * kshell;
      for (G4int jshell = 0; jshell < COUNT_Z; jshell++) {
          pos_shell_z = shell_hz + shell_hz * 2 * jshell;
          for (G4int ishell = 0; ishell < COUNT_Y; ishell++) {
              pos_shell_y = shell_hy + shell_hy * 2 * ishell;
              new G4PVPlacement(0, G4ThreeVector(pos_shell_x, pos_shell_y, pos_shell_z), shellLog, "Shell", logicEnv, false, kshell * COUNT_Z * COUNT_Y + jshell * COUNT_Y + ishell, checkOverlaps);
          }
      }
  }

  //   ***   Sensitive Silicon of BPW34   ***

  G4Material* sensitive_mat = nist->FindOrBuildMaterial("G4_Si");

  // sensitive x position will be -1.6mm deep in PE (1.7*mm if middle is referenced)
  G4double pos_sens_x = -thickness_param-2*sens_hx*mm+sens_hx+shell_hx;
  G4double pos_sens_y = 0 * mm;
  G4double pos_sens_z = 0 * mm;

  G4Box* sensitiveBox = new G4Box("Sensitive", sens_hx, sens_hy, sens_hz);
  G4LogicalVolume* sensitiveLog = new G4LogicalVolume(sensitiveBox, sensitive_mat, "SensitiveLV");

  //   ***   Placing Sensitives of Photodiodes   ***

  new G4PVPlacement(0, G4ThreeVector(pos_sens_x, pos_sens_y, pos_sens_z), sensitiveLog, "Sensitive", shellLog, false, 0, checkOverlaps);


  //G4Box* extraBox = new G4Box("Extra", 0.5 * mm, shell_hy, shell_hz);
  //G4LogicalVolume* extraLog = new G4LogicalVolume(extraBox, LiF_mat, "Extra");
  //new G4PVPlacement(0, G4ThreeVector(-0.5 * mm, shell_hy, shell_hz), extraLog, "Extra", logicEnv, false, 0, checkOverlaps);

  // Set silicon as scoring volume
  //
  fScoringVolume = sensitiveLog;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
