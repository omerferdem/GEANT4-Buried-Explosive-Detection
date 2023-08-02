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
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = false;

  //     
  // World
  //
  G4double world_sizeX = 50 * cm;
  G4double world_sizeY = 160 * cm;
  G4double world_sizeZ = 120 * cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeX, world_sizeY, world_sizeZ);     //its size
      
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
  // Water_Shield
  //
  G4Material* h2o_mat = nist->FindOrBuildMaterial("G4_WATER");

  G4double h2o_hx = 50 * cm;
  G4double h2o_hy = 30 * cm;
  G4double h2o_hz = 20 * cm;

  G4double h2o_posx = 0.0 * cm;
  G4double h2o_posy = 70.0 * cm;
  G4double h2o_posz = -45.0 * cm;

  G4Box* h2oBox = new G4Box("H2Oshield", h2o_hx, h2o_hy, h2o_hz);
  G4LogicalVolume* h2oLog = new G4LogicalVolume(h2oBox, h2o_mat, "H2Oshield");
  new G4PVPlacement(0, G4ThreeVector(h2o_posx, h2o_posy, h2o_posz), h2oLog, "H2Oshield", logicWorld, false, 0, checkOverlaps);

  //     
  // Boron_Shield
  //
  G4Material* b4c_mat = nist->FindOrBuildMaterial("G4_BORON_CARBIDE");

  G4double b4c_hx = 50 * cm;
  G4double b4c_hy = 30 * cm;
  G4double b4c_hz = 2.5 * cm;

  G4double b4c_posx = 0.0 * cm;
  G4double b4c_posy = 70.0 * cm;
  G4double b4c_posz = -22.5 * cm;

  G4Box* b4cBox = new G4Box("B4Cshield", b4c_hx, b4c_hy, b4c_hz);
  G4LogicalVolume* b4cLog = new G4LogicalVolume(b4cBox, b4c_mat, "B4Cshield");
  new G4PVPlacement(0, G4ThreeVector(b4c_posx, b4c_posy, b4c_posz), b4cLog, "B4Cshield", logicWorld, false, 0, checkOverlaps);

  //     
  // Pb_Shield
  //
  G4Material* pb_mat = nist->FindOrBuildMaterial("G4_Pb");

  G4double Pb_hx = 50 * cm;
  G4double Pb_hy = 30 * cm;
  G4double Pb_hz = 10 * cm;

  G4double Pb_posx = 0.0 * cm;
  G4double Pb_posy = 70.0 * cm;
  G4double Pb_posz = -10.0 * cm;

  G4Box* pbBox = new G4Box("Pbshield", Pb_hx, Pb_hy, Pb_hz);
  G4LogicalVolume* pbLog = new G4LogicalVolume(pbBox, pb_mat, "Pbshield");
  new G4PVPlacement(0, G4ThreeVector(Pb_posx, Pb_posy, Pb_posz), pbLog, "Pbshield", logicWorld, false, 0, checkOverlaps);

  //
  // Soil
  //
  G4Material* sio2_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4double air_density = 0.001225 * g / cm3;
  G4double water_density = 0.997 * g / cm3;
  G4double sio2_density = 2.65 * g / cm3;
  G4double fractionmass, ncomponents;

  G4Material* soil_0water = new G4Material("soil_0water", air_density * 0.5 + sio2_density * 0.5, ncomponents = 2);
  soil_0water->AddMaterial(sio2_mat, fractionmass = ((sio2_density * 0.5) / (air_density * 0.5 + sio2_density * 0.5)) * perCent);
  soil_0water->AddMaterial(world_mat, fractionmass = ((air_density * 0.5) / (air_density * 0.5 + sio2_density * 0.5)) * perCent);

  G4Material* soil_10water = new G4Material("soil_10water", water_density * 0.1 + air_density * 0.4 + sio2_density * 0.5, ncomponents = 3);
  soil_10water->AddMaterial(sio2_mat, fractionmass = ((sio2_density * 0.5) / (water_density * 0.1 + air_density * 0.4 + sio2_density * 0.5)) * perCent);
  soil_10water->AddMaterial(world_mat, fractionmass = ((air_density * 0.4) / (water_density * 0.1 + air_density * 0.4 + sio2_density * 0.5)) * perCent);
  soil_10water->AddMaterial(h2o_mat, fractionmass = ((water_density * 0.1) / (water_density * 0.1 + air_density * 0.4 + sio2_density * 0.5)) * perCent);

  G4Material* soil_20water = new G4Material("soil_20water", water_density * 0.2 + air_density * 0.3 + sio2_density * 0.5, ncomponents = 3);
  soil_20water->AddMaterial(sio2_mat, fractionmass = ((sio2_density * 0.5) / (water_density * 0.2 + air_density * 0.3 + sio2_density * 0.5)) * perCent);
  soil_20water->AddMaterial(world_mat, fractionmass = ((air_density * 0.3) / (water_density * 0.2 + air_density * 0.3 + sio2_density * 0.5)) * perCent);
  soil_20water->AddMaterial(h2o_mat, fractionmass = ((water_density * 0.2) / (water_density * 0.2 + air_density * 0.3 + sio2_density * 0.5)) * perCent);

  G4double soil_hx = 50 * cm;
  G4double soil_hy = 80 * cm;
  G4double soil_hz = 120 * cm;

  G4double soil_posx = 0.0 * cm;
  G4double soil_posy = -80.0 * cm;
  G4double soil_posz = 0.0 * cm;

  G4Box* soilBox = new G4Box("Soil", soil_hx, soil_hy, soil_hz);
  G4LogicalVolume* soilLog = new G4LogicalVolume(soilBox, soil_20water, "Soil");
  new G4PVPlacement(0, G4ThreeVector(soil_posx, soil_posy, soil_posz), soilLog, "Soil", logicWorld, false, 0, checkOverlaps);

  //
  // Explosive
  //
  
  G4double exp_depth = 70. * cm;

  G4String name, symbol;
  G4double z, a, density, numberOfatoms, n, iz, abundance, numisotopes;
  
  G4Element* Nitrogen = new G4Element(name = "N", symbol = "N", z = 7, a = 14.0067 * g / mole);
  G4Element* Hidrogen = new G4Element(name = "Hidrogen", symbol = "H", z = 1, a = 1.00784 * g / mole);
  G4Element* Oxygen = new G4Element(name = "Oxygen", symbol = "O", z = 8, a = 15.999 * g / mole);
  G4Material* NH4NO3_mat = new G4Material(name = "AN", density = 1.72 * g / cm3, ncomponents = 3);
  NH4NO3_mat->AddElement(Nitrogen, fractionmass = 0.35);
  NH4NO3_mat->AddElement(Hidrogen, fractionmass = 0.05);
  NH4NO3_mat->AddElement(Oxygen, fractionmass = 0.60);

  G4double kg500_exp_xyz = 66.3 * cm;
  G4double ton1_exp_xyz = 83.5 * cm;
  G4double ton2_exp_xyz = 105. * cm;
  G4Box* expBox = new G4Box("Explosive", kg500_exp_xyz /2, kg500_exp_xyz /2, kg500_exp_xyz /2);
  G4LogicalVolume* expLog = new G4LogicalVolume(expBox, NH4NO3_mat, "Explosive");
  new G4PVPlacement(0, G4ThreeVector(0*cm, -soil_posy- kg500_exp_xyz /2- exp_depth, 0 * cm), expLog, "Explosive", soilLog, false, 0, checkOverlaps);
  
  //
  //CONVERTER
  //

  G4Isotope* Li6 = new G4Isotope(name = "Li6", iz = 3, n = 6);
  G4Isotope* Li7 = new G4Isotope(name = "Li7", iz = 3, n = 7);
  G4Element* Enr_Lithium = new G4Element(name = "Enriched_Lithium", symbol = "Li", numisotopes = 2);
  Enr_Lithium->AddIsotope(Li6, abundance = 90. * perCent);
  Enr_Lithium->AddIsotope(Li7, abundance = 10. * perCent);
  G4Element* F = new G4Element(name = "Flourine", symbol = "F", z = 9.00, a = 18.998 * g / mole);
  G4Material* Enr_LiF = new G4Material(name = "Enriched_LiF", density = 2.64 * g / cm3, ncomponents = 2);
  Enr_LiF->AddElement(Enr_Lithium, numberOfatoms = 1);
  Enr_LiF->AddElement(F, numberOfatoms = 1);

  G4double COUNT_X = 60;
  G4double COUNT_Y = 1;
  G4double COUNT_Z = 72;

  G4double shell_hx = 1.514 * mm;
  G4double shell_hy = 0.114 * mm;
  G4double shell_hz = 1.264 * mm;

  G4double init_x = 0.0 * cm - COUNT_X * shell_hx;
  G4double init_y = 50.0 * cm - COUNT_Y * shell_hy;
  G4double init_z = 30.0 * cm - COUNT_Z * shell_hz;
  G4double pos_shell_x = 0.0 * cm;
  G4double pos_shell_y = 0.0 * cm;
  G4double pos_shell_z = 0.0 * cm;

  G4Box* shellBox = new G4Box("Shell", shell_hx, shell_hy, shell_hz);
  G4LogicalVolume* shellLog = new G4LogicalVolume(shellBox, Enr_LiF, "Shell");

  //   ***   Placing Shells of Photodiodes   ***

  for (G4int kshell = 0; kshell < COUNT_X; kshell++) {
      pos_shell_x = init_x + shell_hx + shell_hx * 2 * kshell;
      for (G4int jshell = 0; jshell < COUNT_Z; jshell++) {
          pos_shell_z = init_y + shell_hz + shell_hz * 2 * jshell;
          for (G4int ishell = 0; ishell < COUNT_Y; ishell++) {
              pos_shell_y = init_z + shell_hy + shell_hy * 2 * ishell;
              new G4PVPlacement(0, G4ThreeVector(pos_shell_x, pos_shell_y, pos_shell_z), shellLog, "Shell", logicWorld, false, kshell * COUNT_Z * COUNT_Y + jshell * COUNT_Y + ishell, checkOverlaps);
          }
      }
  }

  //   ***   Sensitive Silicon of BPW34   ***

  G4Material* sensitive_mat = nist->FindOrBuildMaterial("G4_Si");

  G4double sens_hx = 1.5 * mm;
  G4double sens_hy = 0.1 * mm;
  G4double sens_hz = 1.25 * mm;
  // sensitive x position will be -1.6mm deep in PE (1.7*mm if middle is referenced)
  G4double pos_sens_x = 0;
  G4double pos_sens_y = 0;
  G4double pos_sens_z = 0;

  G4Box* sensitiveBox = new G4Box("Sensitive", sens_hx, sens_hy, sens_hz);
  G4LogicalVolume* sensitiveLog = new G4LogicalVolume(sensitiveBox, sensitive_mat, "Sensitive");
  new G4PVPlacement(0, G4ThreeVector(pos_sens_x, pos_sens_y, pos_sens_z), sensitiveLog, "Sensitive", shellLog, false, 0, checkOverlaps);

  // Set silicon as scoring volume
  //
  fScoringVolume = sensitiveLog;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
