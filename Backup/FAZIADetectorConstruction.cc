#include "FAZIADetectorConstruction.hh"
#include "G4SubtractionSolid.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformElectricField.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include <iostream>
#include "TMath.h"
#define Blktheta 0.069886

using std::stringstream;
using namespace std;
using namespace CLHEP;
  FAZIADetectorConstruction::FAZIADetectorConstruction()
: G4VUserDetectorConstruction()
{
}

FAZIADetectorConstruction::~FAZIADetectorConstruction()
{
}

G4VPhysicalVolume* FAZIADetectorConstruction::Construct()
{
  //----------------Material --------------------
  G4NistManager* nist = G4NistManager::Instance();

  const G4double labTemp = CLHEP::STP_Temperature + 20.*kelvin;

  G4Element*  elCs  = new G4Element("Caesium","Cs",  55,   132.90545*g/mole);
  G4Element*  elI   = new G4Element("Iodine", "I" , 53,   126.90447*g/mole);

  G4Material* CsI = new G4Material("CsI", 4.51*g/CLHEP::cm3, 2, kStateSolid, labTemp);
  CsI -> AddElement(elCs, 1);
  CsI -> AddElement(elI,  1);

  // ------------ Generate & Add Material Properties Table ------------
  //
  G4double photonEnergy[] =
  {	1.90738*eV, 1.98368*eV, 2.06633*eV,
	2.15617*eV, 2.25418*eV, 2.36152*eV,
	2.47960*eV, 2.61011*eV, 2.75511*eV,
	2.91718*eV, 3.0995*eV};
  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

  // CsI
  //fixed
  G4double refractiveIndex1[] =
  { 1.7789, 1.7820, 1.7856, 1.7897, 1.7945,
	1.8001, 1.8067, 1.8146, 1.8242, 1.8359,
	1.8506};

  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

  G4double absorption[] =
  {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
	30.000*m, 28.500*m, 27.000*m,17.500*m, 14.500*m };

  assert(sizeof(absorption) == sizeof(photonEnergy));

  G4double scintilFast[] =
  { 1.00, 1.00, 1.00, 
	1.00, 1.00, 1.00, 
	1.00,1.00, 1.00,
	1.00, 1.00 };

  assert(sizeof(scintilFast) == sizeof(photonEnergy));

  G4double scintilSlow[] =
  { 0.01, 1.00, 2.00,
	3.00, 4.00, 5.00,
	6.00,7.00, 6.00,
	5.00, 4.00 };

  assert(sizeof(scintilSlow) == sizeof(photonEnergy));

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
	->SetSpline(true);
  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
	->SetSpline(true);
  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
	->SetSpline(true);
  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
	->SetSpline(true);

  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
  myMPT1->AddConstProperty("YIELDRATIO",0.8);

  CsI->SetMaterialPropertiesTable(myMPT1);

  // Set the Birks Constant for the Water scintillator

  CsI->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  // -----------------------------------------------------
  // World

  G4Material* world_mat = nist -> FindOrBuildMaterial("G4_Galactic");
  //G4Material* world_mat = nist -> FindOrBuildMaterial("G4_AIR");		
  G4double world_size = 10*m;

  G4Box* solidWorld =
	new G4Box("World",                       // its name
		3*world_size,                // half x
		3*world_size,                // half y
		3*world_size);               // half zNDDetectorConstruction::LANDDetectorConstruction()


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
		true);                 //overlaps checking


  // -----------------------------------------------------
  // Detector
  /*Double_t side_si = 2cm;

	Double_t side_csi_back = 2.272cm;

	Double_t inter_si = 0.24cm;

	Double_t thick_si1 = 300 *um;
	Double_t thick_si2 = 500 *um;
	Double_t thick_csi = 10cm;*/
  const int qua = 1;
  const int tel = 4;

  G4Material* scintillator_mat = nist -> FindOrBuildMaterial("CsI");
  G4Material* Silicon_mat = nist -> FindOrBuildMaterial("G4_Si");
  G4Material* randomizer_mat = nist -> FindOrBuildMaterial("G4_Al");

  G4Box*randomizer = new G4Box("target",(20*mm)/2,(20*mm)/2,(0.1*mm)/2);
  G4LogicalVolume* logicalDetector0 
  = new G4LogicalVolume(randomizer,randomizer_mat,"target");

  G4VPhysicalVolume* PhysicalDetector0 
  = new G4PVPlacement(0,
  G4ThreeVector(0*mm,0*mm,0*mm),
  logicalDetector0,
  "target",
  logicWorld,
  false,
  1,
  true);
  /*string silicon_1[1][4][4] = { {"S1Q1T1","S1Q1T2","S1Q1T4","S1Q1T3"},
	{"S1Q2T1","S1Q2T2","S1Q2T4","S1Q2T3"},
	{"S1Q4T1","S1Q4T2","S1Q4T4","S1Q4T3"},
	{"S1Q3T1","S1Q3T2","S1Q3T4","S1Q3T3"} }
	string silicon_2[1][4][4] = { {"S2Q1T1","S2Q1T2","S2Q1T4","S2Q1T3"},
	{"S2Q2T1","S2Q2T2","S2Q2T4","S2Q2T3"},
	{"S2Q4T1","S2Q4T2","S2Q4T4","S2Q4T3"},
	{"S2Q3T1","S2Q3T2","S2Q3T4","S2Q3T3"} }
	string csicrysial[1][4][4] = { {"CsIQ1T1","CsIQ1T2","CsIQ1T4","CsIQ1T3"},
	{"CsIQ2T1","CsIQ2T2","CsIQ2T4","CsIQ2T3"},
	{"CsIQ4T1","CsIQ4T2","CsIQ4T4","CsIQ4T3"},
	{"CsIQ3T1","CsIQ3T2","CsIQ3T4","CsIQ3T3"} }*/
  /*G4int silicon_1Num[1][qua][tel] = {{{11,12,13,14}}};
  G4int silicon_2Num[1][qua][tel] = {{{2011,2012,2013,2014}}};
  G4int csicrysialNum[1][qua][tel] = {{{3011,3012,3013,3014}}};

  string silicon_1[1][qua][tel] = {{{"S1Q1T1","S1Q1T2","S1Q1T4","S1Q1T3"}}};
  string silicon_2[1][qua][tel] = {{{"S2Q1T1","S2Q1T2","S2Q1T4","S2Q1T3"}}};
  string csicrystal[1][qua][tel] = {{{"CsIQ1T1","CsIQ1T2","CsIQ1T4","CsIQ1T3"}}};
  string TubeSurface[1][qua][tel] = {{{"CsI11","CsI12","CsI13","CsI14"}}};

  G4Box* Silicon_1[1][qua][tel];
  G4Box* Silicon_2[1][qua][tel];
  G4Box* CsIcrystal[1][qua][tel];
  //G4Trd* CsIcrystal[1][qua][tel];
  
  G4LogicalVolume* logicalDetector1[1][qua][tel];
  G4LogicalVolume* logicalDetector2[1][qua][tel];
  G4LogicalVolume* logicalDetector3[1][qua][tel];
  
  G4VPhysicalVolume* PhysicalDetector1[1][qua][tel];
  G4VPhysicalVolume* PhysicalDetector2[1][qua][tel];
  G4VPhysicalVolume* PhysicalDetector3[1][qua][tel];
  
  G4OpticalSurface* opTubeSurface[1][qua][tel];*/
/*new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*degree, Angle1, 0*degree), G4ThreeVector(0*mm, 0*mm, 400*mm).rotateX(Angle1).rotateZ((360/CsIPhiN_1)*copyN*degree),"Si_PV", Si_LV_1, labPV, false, 111+copyN);
1000*(TMath::Sin(0.00142*(Blk+1)))*(TMath::Cos((pi/4.0)+Tel*(pi/2.0))*mm) , 1000*(TMath::Sin(0.00142*(Blk+1)))*(TMath::Sin((pi/4.0)+Tel*(pi/2.0)))*mm , 1000*(TMath::Cos(0.00142*(Blk+1)))*mm)*/
	//for(G4int Tel = 0; Tel < 4 ; Tel++){
	//}
  
  G4int silicon_11Num[tel] = {11,22,33,44};
  G4int silicon_12Num[tel] = {2011,2022,2033,2044};
  G4int csicrysial1Num[tel] = {3011,3022,3033,3044};
  
  G4int silicon_21Num[tel] = {12,23,34,41};
  G4int silicon_22Num[tel] = {2012,2023,2034,2041};
  G4int csicrysial2Num[tel] = {3012,3023,3034,3041};
  
  G4int silicon_31Num[tel] = {14,21,32,43};
  G4int silicon_32Num[tel] = {2014,2021,2032,2043};
  G4int csicrysial3Num[tel] = {3014,3021,3032,3043};
  
  G4int silicon_41Num[tel] = {13,24,31,42};
  G4int silicon_42Num[tel] = {2013,2024,2031,2042};
  G4int csicrysial4Num[tel] = {3013,3024,3031,3042};

  string silicon_11[tel] = {"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"};
  string silicon_12[tel] = {"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"};
  string csicrystal1[tel] = {"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"};
  
  string silicon_21[tel] = {"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"};
  string silicon_22[tel] = {"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"};
  string csicrystal2[tel] = {"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"};
  
  string silicon_31[tel] = {"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"};
  string silicon_32[tel] = {"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"};
  string csicrystal3[tel] = {"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"};
  
  string silicon_41[tel] = {"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"};
  string silicon_42[tel] = {"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"};
  string csicrystal4[tel] = {"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"};
  string TubeSurface[tel] = {"CsI11","CsI12","CsI13","CsI14"};

  G4Box* Silicon_11[tel];
  G4Box* Silicon_12[tel];
  G4Trd* CsIcrystal1[tel];
  
  G4Box* Silicon_21[tel];
  G4Box* Silicon_22[tel];
  G4Trd* CsIcrystal2[tel];
  
  G4Box* Silicon_31[tel];
  G4Box* Silicon_32[tel];
  G4Trd* CsIcrystal3[tel];
  
  G4Box* Silicon_41[tel];
  G4Box* Silicon_42[tel];
  G4Trd* CsIcrystal4[tel];
  //G4Trd* CsIcrystal[1][qua][tel];
  
  G4LogicalVolume* logicalDetector11[tel];
  G4LogicalVolume* logicalDetector12[tel];
  G4LogicalVolume* logicalDetector13[tel];
  G4VPhysicalVolume* PhysicalDetector11[tel];
  G4VPhysicalVolume* PhysicalDetector12[tel];
  G4VPhysicalVolume* PhysicalDetector13[tel];
  
  G4LogicalVolume* logicalDetector21[tel];
  G4LogicalVolume* logicalDetector22[tel];
  G4LogicalVolume* logicalDetector23[tel];
  G4VPhysicalVolume* PhysicalDetector21[tel];
  G4VPhysicalVolume* PhysicalDetector22[tel];
  G4VPhysicalVolume* PhysicalDetector23[tel];
  
  G4LogicalVolume* logicalDetector31[tel];
  G4LogicalVolume* logicalDetector32[tel];
  G4LogicalVolume* logicalDetector33[tel];
  G4VPhysicalVolume* PhysicalDetector31[tel];
  G4VPhysicalVolume* PhysicalDetector32[tel];
  G4VPhysicalVolume* PhysicalDetector33[tel];
  
  G4LogicalVolume* logicalDetector41[tel];
  G4LogicalVolume* logicalDetector42[tel];
  G4LogicalVolume* logicalDetector43[tel];
  G4VPhysicalVolume* PhysicalDetector41[tel];
  G4VPhysicalVolume* PhysicalDetector42[tel];
  G4VPhysicalVolume* PhysicalDetector43[tel];
  
  G4OpticalSurface* opTubeSurface[tel];

 // for(G4int Blk = 0 ; Blk < 1 ; Blk++)
  //{
	//for(G4int Qua = 0 ; Qua < 1 ; Qua++)
	//{
    G4double AngleN = -(0.810+4)*CLHEP::degree;
    G4double AngleNN = -(0.80+4)*CLHEP::degree;
  	G4RotationMatrix* Rot_1Si[tel];
  	G4RotationMatrix* Rot_2Si[tel];
  	G4RotationMatrix* Rot_3Si[tel];
  	G4RotationMatrix* Rot_4Si[tel];
	  for(G4int Tel = 0 ; Tel < 4 ; Tel++)
	  {
	Rot_1Si[Tel] = new G4RotationMatrix;
	Rot_2Si[Tel] = new G4RotationMatrix;
	Rot_3Si[Tel] = new G4RotationMatrix;
	Rot_4Si[Tel] = new G4RotationMatrix;
	if(Tel == 0){
  	Rot_1Si[Tel] -> rotateY(AngleNN*2);
  	Rot_1Si[Tel] -> rotateX(AngleNN*-2);}
	if(Tel == 1){
  	Rot_1Si[Tel] -> rotateY(AngleNN*-2);
  	Rot_1Si[Tel] -> rotateX(AngleNN*-2);}
	if(Tel == 2){
  	Rot_1Si[Tel] -> rotateY(AngleNN*-2);
  	Rot_1Si[Tel] -> rotateX(AngleNN*2);}
	if(Tel == 3){
  	Rot_1Si[Tel] -> rotateY(AngleNN*2);
  	Rot_1Si[Tel] -> rotateX(AngleNN*2);}
	
	if(Tel == 0){
  	Rot_2Si[Tel] -> rotateY(AngleN*1);
  	Rot_2Si[Tel] -> rotateX(AngleNN*-2);}
	if(Tel == 1){
  	Rot_2Si[Tel] -> rotateY(AngleNN*-2);
  	Rot_2Si[Tel] -> rotateX(AngleN*-1);}
	if(Tel == 2){
  	Rot_2Si[Tel] -> rotateY(AngleN*-1);
  	Rot_2Si[Tel] -> rotateX(AngleNN*2);}
	if(Tel == 3){
  	Rot_2Si[Tel] -> rotateY(AngleNN*2);
  	Rot_2Si[Tel] -> rotateX(AngleN*1);}
	
	if(Tel == 0){
  	Rot_3Si[Tel] -> rotateY(AngleNN*2);
  	Rot_3Si[Tel] -> rotateX(AngleN*-1);}
	if(Tel == 1){
  	Rot_3Si[Tel] -> rotateY(AngleN*-1);
  	Rot_3Si[Tel] -> rotateX(AngleNN*-2);}
	if(Tel == 2){
  	Rot_3Si[Tel] -> rotateY(AngleNN*-2);
  	Rot_3Si[Tel] -> rotateX(AngleN*1);}
	if(Tel == 3){
  	Rot_3Si[Tel] -> rotateY(AngleN*1);
  	Rot_3Si[Tel] -> rotateX(AngleNN*2);}
	
	if(Tel == 0){
  	Rot_4Si[Tel] -> rotateY(AngleN*1);
  	Rot_4Si[Tel] -> rotateX(AngleN*-1);}
	if(Tel == 1){
  	Rot_4Si[Tel] -> rotateY(AngleN*-1);
  	Rot_4Si[Tel] -> rotateX(AngleN*-1);}
	if(Tel == 2){
  	Rot_4Si[Tel] -> rotateY(AngleN*-1);
  	Rot_4Si[Tel] -> rotateX(AngleN*1);}
	if(Tel == 3){
  	Rot_4Si[Tel] -> rotateY(AngleN*1);
  	Rot_4Si[Tel] -> rotateX(AngleN*1);}
//one block's geometry----------------------------------------------------------------------------------------------------------		
		{//1
		//Si1
		Silicon_11[Tel] = new G4Box(silicon_11[Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector11[Tel] = new G4LogicalVolume(Silicon_11[Tel],Silicon_mat,silicon_11[Tel]);
		  PhysicalDetector11[Tel] = new G4PVPlacement(Rot_1Si[Tel],
			  G4ThreeVector(1000*(TMath::Sin(0.0450+Blktheta))*TMath::Cos((3.141592/4.0))*mm, 1000*(TMath::Sin(0.0450+Blktheta))*TMath::Sin((3.141592/4.0))*mm, 1000*(TMath::Cos(0.0450+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector11[Tel],//theta == 0.0144
			  silicon_11[Tel],
			  logicWorld,
			  false,
			  silicon_11Num[Tel],
			  true);

		  //Si2		 
		 Silicon_12[Tel] = new G4Box(silicon_12[Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		  logicalDetector12[Tel] = new G4LogicalVolume(Silicon_12[Tel] ,Silicon_mat ,silicon_12[Tel]);

		  PhysicalDetector12[Tel] = new G4PVPlacement(Rot_1Si[Tel],
			  G4ThreeVector(1000.625*(TMath::Sin(0.0450+Blktheta))*TMath::Cos((3.141592/4.0))*mm, 1000.625*(TMath::Sin(0.0450+Blktheta))*TMath::Sin((3.141592/4.0))*mm, 1000.625*(TMath::Cos(0.0450+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector12[Tel],
			  silicon_12[Tel],
			  logicWorld,
			  false,
			  silicon_12Num[Tel],
			  true);

		  //CsI crystal		  
		  CsIcrystal1[Tel] = new G4Trd(csicrystal1[Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		  logicalDetector13[Tel] = new G4LogicalVolume(CsIcrystal1[Tel] ,scintillator_mat ,csicrystal1[Tel]);

		  PhysicalDetector13[Tel] = new G4PVPlacement(Rot_1Si[Tel],
			  G4ThreeVector(1051.058*(TMath::Sin(0.0450+Blktheta))*TMath::Cos((3.141592/4.0))*mm, 1051.058*(TMath::Sin(0.0450+Blktheta))*TMath::Sin((3.141592/4.0))*mm, 1051.058*(TMath::Cos(0.0450+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector13[Tel],
			  csicrystal1[Tel],
			  logicWorld,
			  false,
			  csicrysial1Num[Tel],
			  true);
		}
		{//2
		Silicon_21[Tel] = new G4Box(silicon_21[Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector21[Tel] = new G4LogicalVolume(Silicon_21[Tel],Silicon_mat,silicon_21[Tel]);
		  PhysicalDetector21[Tel] = new G4PVPlacement(Rot_2Si[Tel],
			  G4ThreeVector(1000*(TMath::Sin(0.0330+Blktheta))*TMath::Cos((1.249))*mm, 1000*(TMath::Sin(0.0330+Blktheta))*TMath::Sin((1.249))*mm, 1000*(TMath::Cos(0.0330+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector21[Tel],//theta == 0.0316
			  silicon_21[Tel],
			  logicWorld,
			  false,
			  silicon_21Num[Tel],
			  true);

		  //Si2		 
		 Silicon_22[Tel] = new G4Box(silicon_22[Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		  logicalDetector22[Tel] = new G4LogicalVolume(Silicon_22[Tel] ,Silicon_mat ,silicon_22[Tel]);

		  PhysicalDetector22[Tel] = new G4PVPlacement(Rot_2Si[Tel],
			  G4ThreeVector(1000.625*(TMath::Sin(0.0330+Blktheta))*TMath::Cos((1.249))*mm, 1000.625*(TMath::Sin(0.0330+Blktheta))*TMath::Sin((1.249))*mm, 1000.625*(TMath::Cos(0.0330+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector22[Tel],
			  silicon_22[Tel],
			  logicWorld,
			  false,
			  silicon_22Num[Tel],
			  true);

		  //CsI crystal		  
		  CsIcrystal2[Tel] = new G4Trd(csicrystal2[Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		  logicalDetector23[Tel] = new G4LogicalVolume(CsIcrystal2[Tel] ,scintillator_mat ,csicrystal2[Tel]);

		  PhysicalDetector23[Tel] = new G4PVPlacement(Rot_2Si[Tel],
			  G4ThreeVector(1051.058*(TMath::Sin(0.0330+Blktheta))*TMath::Cos((1.249))*mm, 1051.058*(TMath::Sin(0.0330+Blktheta))*TMath::Sin((1.249))*mm, 1051.058*(TMath::Cos(0.0330+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector23[Tel],
			  csicrystal2[Tel],
			  logicWorld,
			  false,
			  csicrysial2Num[Tel],
			  true);
		}
		{//3
		Silicon_31[Tel] = new G4Box(silicon_31[Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector31[Tel] = new G4LogicalVolume(Silicon_31[Tel],Silicon_mat,silicon_31[Tel]);
		  PhysicalDetector31[Tel] = new G4PVPlacement(Rot_3Si[Tel],
			  G4ThreeVector(1000*(TMath::Sin(0.0330+Blktheta))*TMath::Cos((0.321))*mm, 1000*(TMath::Sin(0.0330+Blktheta))*TMath::Sin((0.321))*mm, 1000*(TMath::Cos(0.0330+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector31[Tel],//theta == 0.0147
			  silicon_31[Tel],
			  logicWorld,
			  false,
			  silicon_31Num[Tel],
			  true);

		  //Si2		 
		 Silicon_32[Tel] = new G4Box(silicon_32[Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		  logicalDetector32[Tel] = new G4LogicalVolume(Silicon_32[Tel] ,Silicon_mat ,silicon_32[Tel]);

		  PhysicalDetector32[Tel] = new G4PVPlacement(Rot_3Si[Tel],
			  G4ThreeVector(1000.625*(TMath::Sin(0.0330+Blktheta))*TMath::Cos((0.321))*mm, 1000.625*(TMath::Sin(0.0330+Blktheta))*TMath::Sin((0.321))*mm, 1000.625*(TMath::Cos(0.0330+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector32[Tel],
			  silicon_32[Tel],
			  logicWorld,
			  false,
			  silicon_32Num[Tel],
			  true);

		  //CsI crystal		  
		  CsIcrystal3[Tel] = new G4Trd(csicrystal3[Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		  logicalDetector33[Tel] = new G4LogicalVolume(CsIcrystal3[Tel] ,scintillator_mat ,csicrystal3[Tel]);

		  PhysicalDetector33[Tel] = new G4PVPlacement(Rot_3Si[Tel],
			  G4ThreeVector(1051.058*(TMath::Sin(0.0330+Blktheta))*TMath::Cos((0.321))*mm, 1051.058*(TMath::Sin(0.0330+Blktheta))*TMath::Sin((0.321))*mm, 1051.058*(TMath::Cos(0.0330+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector33[Tel],
			  csicrystal3[Tel],
			  logicWorld,
			  false,
			  csicrysial3Num[Tel],
			  true);
		}
		{//4
		Silicon_41[Tel] = new G4Box(silicon_41[Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector41[Tel] = new G4LogicalVolume(Silicon_41[Tel],Silicon_mat,silicon_41[Tel]);
		  PhysicalDetector41[Tel] = new G4PVPlacement(Rot_4Si[Tel],
			  G4ThreeVector(1000*(TMath::Sin(0.0146+Blktheta))*TMath::Cos((3.141592/4.0))*mm, 1000*(TMath::Sin(0.0146+Blktheta))*TMath::Sin((3.141592/4.0))*mm, 1000*(TMath::Cos(0.0146+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector41[Tel],//theta == 0.0147
			  silicon_41[Tel],
			  logicWorld,
			  false,
			  silicon_41Num[Tel],
			  true);

		  //Si2		 
		 Silicon_42[Tel] = new G4Box(silicon_42[Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		  logicalDetector42[Tel] = new G4LogicalVolume(Silicon_42[Tel] ,Silicon_mat ,silicon_42[Tel]);

		  PhysicalDetector42[Tel] = new G4PVPlacement(Rot_4Si[Tel],
			  G4ThreeVector(1000.625*(TMath::Sin(0.0146+Blktheta))*TMath::Cos((3.141592/4.0))*mm, 1000.625*(TMath::Sin(0.0146+Blktheta))*TMath::Sin((3.141592/4.0))*mm, 1000.625*(TMath::Cos(0.0146+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector42[Tel],
			  silicon_42[Tel],
			  logicWorld,
			  false,
			  silicon_42Num[Tel],
			  true);

		  //CsI crystal		  
		  CsIcrystal4[Tel] = new G4Trd(csicrystal4[Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		  logicalDetector43[Tel] = new G4LogicalVolume(CsIcrystal4[Tel] ,scintillator_mat ,csicrystal4[Tel]);

		  PhysicalDetector43[Tel] = new G4PVPlacement(Rot_4Si[Tel],
			  G4ThreeVector(1051.058*(TMath::Sin(0.0146+Blktheta))*TMath::Cos((3.141592/4.0))*mm, 1051.058*(TMath::Sin(0.0146+Blktheta))*TMath::Sin((3.141592/4.0))*mm, 1051.058*(TMath::Cos(0.0146+Blktheta))*mm).rotateZ(90*(Tel)*degree),
			  logicalDetector43[Tel],
			  csicrystal4[Tel],
			  logicWorld,
			  false,
			  csicrysial4Num[Tel],
			  true);
		}
//-------------------------------------------------------------------------------------------------------------------------------		
		opTubeSurface[Tel] = new G4OpticalSurface(TubeSurface[Tel]);
		opTubeSurface[Tel]->SetType(dielectric_metal);
		opTubeSurface[Tel]->SetFinish(polished);
		opTubeSurface[Tel]->SetModel(unified);

		new G4LogicalBorderSurface(TubeSurface[Tel],PhysicalDetector13[Tel],physWorld,opTubeSurface[Tel]);
		new G4LogicalBorderSurface(TubeSurface[Tel],PhysicalDetector23[Tel],physWorld,opTubeSurface[Tel]);
		new G4LogicalBorderSurface(TubeSurface[Tel],PhysicalDetector33[Tel],physWorld,opTubeSurface[Tel]);
		new G4LogicalBorderSurface(TubeSurface[Tel],PhysicalDetector43[Tel],physWorld,opTubeSurface[Tel]);
  const G4int num = 2;
  G4double ephoton[num] = {2.034*eV, 4.136*eV};

  //Optical alu rap Surface
  G4double refractiveIndex[num] = {0.48, 1.55};
  G4double specularLobe[num]    = {0.3, 0.3};
  G4double specularSpike[num]   = {0.2, 0.2};
  G4double backScatter[num]     = {0.2, 0.2};
  G4double reflecitivity[num]   = {0.90, 0.92};

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

  myST1->AddProperty("RINDEX",                ephoton, refractiveIndex, num);
  myST1->AddProperty("REFLECITIVITY",         ephoton, reflecitivity,   num);
  myST1->AddProperty("SPECULARLOBECONSTANT",  ephoton, specularLobe,    num);
  myST1->AddProperty("SPECULARSPIKECONSTANT", ephoton, specularSpike,   num);
  myST1->AddProperty("BACKSCATTERCONSTANT",   ephoton, backScatter,     num);

  G4cout << "Tube Surface G4MaterialPropertiesTable" << G4endl;
  myST1->DumpTable();

  opTubeSurface[Tel]->SetMaterialPropertiesTable(myST1);

	  }
	//}
  //}

  // ------------- Surfaces --------------
  //
  // Csi crystal
  //

  //G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>


  //if (opticalSurface) opticalSurface->DumpInfo();

  // Generate & Add Material Properties Table attached to the optical surfaces
  //

  //runManager -> SetSensitiveDetector(pvp);

  //Set detector's Color

  logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());

  G4VisAttributes* Silicon_11VisAtt[4];
  G4VisAttributes* Silicon_12VisAtt[4];
  G4VisAttributes* CsIcrystal1VisAtt[4];
  G4VisAttributes* Silicon_21VisAtt[4];
  G4VisAttributes* Silicon_22VisAtt[4];
  G4VisAttributes* CsIcrystal2VisAtt[4];
  G4VisAttributes* Silicon_31VisAtt[4];
  G4VisAttributes* Silicon_32VisAtt[4];
  G4VisAttributes* CsIcrystal3VisAtt[4];
  G4VisAttributes* Silicon_41VisAtt[4];
  G4VisAttributes* Silicon_42VisAtt[4];
  G4VisAttributes* CsIcrystal4VisAtt[4];
  //for(int i = 0 ;i < 1 ;i++)
  //{
	//for(int j = 0 ;j < 4 ;j++)
	//{
	  for(int k = 0 ;k < 4 ;k++)
	  {

		Silicon_11VisAtt[k] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
		Silicon_11VisAtt[k]-> SetForceWireframe(true);
		logicalDetector11[k] -> SetVisAttributes (Silicon_11VisAtt[k]);

		Silicon_12VisAtt[k] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
		Silicon_12VisAtt[k] -> SetForceWireframe(true);
		logicalDetector12[k] -> SetVisAttributes (Silicon_12VisAtt[k]);

		CsIcrystal1VisAtt[k] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
		CsIcrystal1VisAtt[k] -> SetForceWireframe(true);
		logicalDetector13[k] -> SetVisAttributes (CsIcrystal1VisAtt[k]);
	  //}
		Silicon_21VisAtt[k] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
		Silicon_21VisAtt[k]-> SetForceWireframe(true);
		logicalDetector21[k] -> SetVisAttributes (Silicon_21VisAtt[k]);

		Silicon_22VisAtt[k] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
		Silicon_22VisAtt[k] -> SetForceWireframe(true);
		logicalDetector22[k] -> SetVisAttributes (Silicon_22VisAtt[k]);

		CsIcrystal2VisAtt[k] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
		CsIcrystal2VisAtt[k] -> SetForceWireframe(true);
		logicalDetector23[k] -> SetVisAttributes (CsIcrystal2VisAtt[k]);
	//}
		Silicon_31VisAtt[k] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
		Silicon_31VisAtt[k]-> SetForceWireframe(true);
		logicalDetector31[k] -> SetVisAttributes (Silicon_31VisAtt[k]);

		Silicon_32VisAtt[k] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
		Silicon_32VisAtt[k] -> SetForceWireframe(true);
		logicalDetector32[k] -> SetVisAttributes (Silicon_32VisAtt[k]);

		CsIcrystal3VisAtt[k] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
		CsIcrystal3VisAtt[k] -> SetForceWireframe(true);
		logicalDetector33[k] -> SetVisAttributes (CsIcrystal3VisAtt[k]);
 //}		
		Silicon_41VisAtt[k] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
		Silicon_41VisAtt[k]-> SetForceWireframe(true);
		logicalDetector41[k] -> SetVisAttributes (Silicon_41VisAtt[k]);

		Silicon_42VisAtt[k] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
		Silicon_42VisAtt[k] -> SetForceWireframe(true);
		logicalDetector42[k] -> SetVisAttributes (Silicon_42VisAtt[k]);

		CsIcrystal4VisAtt[k] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
		CsIcrystal4VisAtt[k] -> SetForceWireframe(true);
		logicalDetector43[k] -> SetVisAttributes (CsIcrystal4VisAtt[k]);
  }
  return physWorld;
}
/*new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*degree, Angle1, 0*degree), G4ThreeVector(0*mm, 0*mm, 400*mm).rotateX(Angle1).rotateZ((360/CsIPhiN_1)*copyN*degree),"Si_PV", Si_LV_1, labPV, false, 111+copyN);*/
