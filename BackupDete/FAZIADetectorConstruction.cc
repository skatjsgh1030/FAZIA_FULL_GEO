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
#define Blktheta1 0.10728646 //0.10728646
#define Blktheta2 0.10062866 //0.10162866
#define Blktheta3 0.08820 //0.08920
#define Blktheta4 0.0812759 //0.0822759
#define Blkphi1 1.1902899
#define Blkphi2 1.3734008
#define Blkphi3 1.1071487
#define Blkphi4 1.3258177

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
  const int blkcor = 4;
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

  /*new G4PVPlacement(new G4RotationMatrix((360/CsIPhiN_1)*copyN*degree, Angle1, 0*degree), G4ThreeVector(0*mm, 0*mm, 400*mm).rotateX(Angle1).rotateZ((360/CsIPhiN_1)*copyN*degree),"Si_PV", Si_LV_1, labPV, false, 111+copyN);
	1000*(TMath::Sin(0.00142*(Blk+1)))*(TMath::Cos((pi/4.0)+Tel*(pi/2.0))*mm) , 1000*(TMath::Sin(0.00142*(Blk+1)))*(TMath::Sin((pi/4.0)+Tel*(pi/2.0)))*mm , 1000*(TMath::Cos(0.00142*(Blk+1)))*mm)*/
  //for(G4int Tel = 0; Tel < 4 ; Tel++){
  //}

  G4int silicon_11Num[blkcor][tel] = {{11,22,33,44},
	{111,122,133,144},
	{211,222,233,244},
	{311,322,333,344}};

  G4int silicon_12Num[blkcor][tel] = {{2011,2022,2033,2044},
	{2111,2122,2133,2144},
	{2211,2222,2233,2244},
	{2311,2322,2333,2344}};

  G4int csicrysial1Num[blkcor][tel] = {{3011,3022,3033,3044},
	{3111,3122,3133,3144},
	{3211,3322,3233,3244},
	{3311,3322,3333,3344}};

  G4int silicon_21Num[blkcor][tel] = {{12,23,34,41},
	{112,123,134,141},
	{212,223,234,241},
	{312,323,334,341}};

  G4int silicon_22Num[blkcor][tel] = {{2012,2023,2034,2041},
	{2112,2123,2134,2141},
	{2212,2223,2234,2241},
	{2312,2323,2334,2341}};

  G4int csicrysial2Num[blkcor][tel] = {{3012,3023,3034,3041},
	{3112,3123,3134,3141},
	{3212,3223,3234,3241},
	{3312,3323,3334,3341}};

  G4int silicon_31Num[blkcor][tel] = {{14,21,32,43},
	{114,121,132,143},
	{214,221,232,243},
	{314,321,332,343}};

  G4int silicon_32Num[blkcor][tel] = {{2014,2021,2032,2043},
	{2114,2121,2132,2143},
	{2214,2221,2232,2243},
	{2314,2321,2332,2343}};

  G4int csicrysial3Num[blkcor][tel] = {{3014,3021,3032,3043},
	{3114,3121,3132,3143},
	{3214,3221,3232,3243},
	{3314,3321,3332,3343}};

  G4int silicon_41Num[blkcor][tel] = {{13,24,31,42},
	{113,124,131,142},
	{213,224,231,242},
	{313,324,331,342}};

  G4int silicon_42Num[blkcor][tel] = {{2013,2024,2031,2042},
	{2113,2124,2131,2142},
	{2213,2224,2231,2242},
	{2313,2324,2331,2342}};

  G4int csicrysial4Num[blkcor][tel] = {{3013,3024,3031,3042},
	{3113,3124,3131,3142},
	{3213,3224,3231,3242},
	{3313,3324,3331,3342}};

  string silicon_11[blkcor][tel] ={ {"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"},
	{"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"},
	{"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"},
	{"S1Q1T1","S1Q2T2","S1Q3T3","S1Q4T4"}};

  string silicon_12[blkcor][tel] ={ {"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"},
	{"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"},
	{"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"},
	{"S2Q1T1","S2Q2T2","S2Q3T3","S2Q4T4"}};

  string csicrystal1[blkcor][tel] ={ {"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"},
	{"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"},
	{"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"},
	{"CsIQ1T1","CsIQ2T2","CsIQ3T3","CsIQ4T4"}};

  string silicon_21[blkcor][tel] ={ {"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"},
	{"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"},
	{"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"},
	{"S1Q1T2","S1Q2T3","S1Q3T4","S1Q4T1"}};

  string silicon_22[blkcor][tel] ={ {"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"},
	{"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"},
	{"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"},
	{"S2Q1T2","S2Q2T3","S2Q3T4","S2Q4T1"}};

  string csicrystal2[blkcor][tel] = {{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"},
	{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"},
	{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"},
	{"CsIQ1T2","CsIQ2T3","CsIQ3T4","CsIQ4T1"}};

  string silicon_31[blkcor][tel] ={ {"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"},
	{"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"},
	{"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"},
	{"S1Q1T4","S1Q2T1","S1Q3T2","S1Q4T3"}};

  string silicon_32[blkcor][tel] ={ {"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"},
	{"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"},
	{"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"},
	{"S2Q1T4","S2Q2T1","S2Q3T2","S2Q4T3"}};

  string csicrystal3[blkcor][tel] ={ {"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"},
	{"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"},
	{"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"},
	{"CsIQ1T4","CsIQ2T1","CsIQ3T2","CsIQ4T3"}};

  string silicon_41[blkcor][tel] ={ {"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"},
	{"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"},
	{"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"},
	{"S1Q1T3","S1Q2T4","S1Q3T1","S1Q4T2"}};

  string silicon_42[blkcor][tel] ={ {"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"},
	{"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"},
	{"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"},
	{"S2Q1T3","S2Q2T4","S2Q3T1","S2Q4T2"}};

  string csicrystal4[blkcor][tel] ={ {"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"},
	{"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"},
	{"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"},
	{"CsIQ1T3","CsIQ2T4","CsIQ3T1","CsIQ4T2"}};

  string TubeSurface[blkcor][tel] ={ {"CsI11","CsI12","CsI13","CsI14"},
	{"CsI11","CsI12","CsI13","CsI14"},
	{"CsI11","CsI12","CsI13","CsI14"},
	{"CsI11","CsI12","CsI13","CsI14"}};

  G4Box* Silicon_11[blkcor][tel];
  G4Box* Silicon_12[blkcor][tel];
  G4Trd* CsIcrystal1[blkcor][tel];

  G4Box* Silicon_21[blkcor][tel];
  G4Box* Silicon_22[blkcor][tel];
  G4Trd* CsIcrystal2[blkcor][tel];

  G4Box* Silicon_31[blkcor][tel];
  G4Box* Silicon_32[blkcor][tel];
  G4Trd* CsIcrystal3[blkcor][tel];

  G4Box* Silicon_41[blkcor][tel];
  G4Box* Silicon_42[blkcor][tel];
  G4Trd* CsIcrystal4[blkcor][tel];
  //G4Trd* CsIcrystal[1][qua][blkcor][tel];

  G4LogicalVolume* logicalDetector11[blkcor][tel];
  G4LogicalVolume* logicalDetector12[blkcor][tel];
  G4LogicalVolume* logicalDetector13[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector11[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector12[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector13[blkcor][tel];

  G4LogicalVolume* logicalDetector21[blkcor][tel];
  G4LogicalVolume* logicalDetector22[blkcor][tel];
  G4LogicalVolume* logicalDetector23[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector21[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector22[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector23[blkcor][tel];

  G4LogicalVolume* logicalDetector31[blkcor][tel];
  G4LogicalVolume* logicalDetector32[blkcor][tel];
  G4LogicalVolume* logicalDetector33[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector31[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector32[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector33[blkcor][tel];

  G4LogicalVolume* logicalDetector41[blkcor][tel];
  G4LogicalVolume* logicalDetector42[blkcor][tel];
  G4LogicalVolume* logicalDetector43[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector41[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector42[blkcor][tel];
  G4VPhysicalVolume* PhysicalDetector43[blkcor][tel];

  G4OpticalSurface* opTubeSurface[blkcor][tel];

  logicWorld -> SetVisAttributes (G4VisAttributes::GetInvisible());

  G4VisAttributes* Silicon_11VisAtt[4][4];
  G4VisAttributes* Silicon_12VisAtt[4][4];
  G4VisAttributes* CsIcrystal1VisAtt[4][4];
  G4VisAttributes* Silicon_21VisAtt[4][4];
  G4VisAttributes* Silicon_22VisAtt[4][4];
  G4VisAttributes* CsIcrystal2VisAtt[4][4];
  G4VisAttributes* Silicon_31VisAtt[4][4];
  G4VisAttributes* Silicon_32VisAtt[4][4];
  G4VisAttributes* CsIcrystal3VisAtt[4][4];
  G4VisAttributes* Silicon_41VisAtt[4][4];
  G4VisAttributes* Silicon_42VisAtt[4][4];
  G4VisAttributes* CsIcrystal4VisAtt[4][4];
  // for(G4int Blk = 0 ; Blk < 1 ; Blk++)
  //{
  //for(G4int Qua = 0 ; Qua < 1 ; Qua++)
  //{
  G4double AngleN = -(0.810)*CLHEP::degree;
  G4double AngleNN = -(0.80)*CLHEP::degree;

  G4RotationMatrix* Rot_1Si[blkcor][tel];
  G4RotationMatrix* Rot_2Si[blkcor][tel];
  G4RotationMatrix* Rot_3Si[blkcor][tel];
  G4RotationMatrix* Rot_4Si[blkcor][tel];

  //----------------------------------------------------------------------------------------------------------
  for(G4int Blkcor = 0 ; Blkcor < 4 ; Blkcor++)
  {
	for(G4int Tel = 0 ; Tel < 4 ; Tel++)
	{
	  Rot_1Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_2Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_3Si[Blkcor][Tel] = new G4RotationMatrix;
	  Rot_4Si[Blkcor][Tel] = new G4RotationMatrix;

	  if(Tel == 0){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 1){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 2){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*2);}
	  if(Tel == 3){
		Rot_1Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_1Si[Blkcor][Tel] -> rotateX(AngleNN*2);}

	  if(Tel == 0){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 1){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 2){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleNN*2);}
	  if(Tel == 3){
		Rot_2Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_2Si[Blkcor][Tel] -> rotateX(AngleN*1);}

	  if(Tel == 0){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleNN*2);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 1){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleNN*-2);}
	  if(Tel == 2){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleNN*-2);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleN*1);}
	  if(Tel == 3){
		Rot_3Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_3Si[Blkcor][Tel] -> rotateX(AngleNN*2);}

	  if(Tel == 0){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 1){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*-1);}
	  if(Tel == 2){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*-1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*1);}
	  if(Tel == 3){
		Rot_4Si[Blkcor][Tel] -> rotateY(AngleN*1);
		Rot_4Si[Blkcor][Tel] -> rotateX(AngleN*1);}

	  // u, v, w are the daughter axes, projected on the mother frame
	  //G4ThreeVector u = G4ThreeVector(1,0,0);
	  //G4ThreeVector v = G4ThreeVector(0,1,0);
	  //G4ThreeVector w = G4ThreeVector( 0,0,1);

	  //G4RotationMatrix rotm1  = G4RotationMatrix(u, v, w);
	  //G4cout << "\n --> phi = " << phi/deg << " deg;  direct rotation matrix : ";
	  //rotm1.print(G4cout);     
	  //G4ThreeVector position1 = og*w;

	  /*new G4PVPlacement(transform1,         //position, rotation        
		fTrdVolume,         //logical volume
		"Trd",              //name
		fWorldVolume,       //mother volume
		false,              //no boolean operation

		1);                 //copy number        3.14/4*/
	  //one block's geometry----------------------------------------------------------------------------------------------------------		
	  //G4RotationMatrix* rotm = new G4RotationMatrix(u,v,w);
	  {//1
		//Si1
		Silicon_11[Blkcor][Tel] = new G4Box(silicon_11[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector11[Blkcor][Tel] = new G4LogicalVolume(Silicon_11[Blkcor][Tel],Silicon_mat,silicon_11[Blkcor][Tel]);
		PhysicalDetector11[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta1))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector11[Blkcor][Tel],//theta == 0.0144
			silicon_11[Blkcor][Tel],
			logicWorld,
			false,
			silicon_11Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_12[Blkcor][Tel] = new G4Box(silicon_12[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector12[Blkcor][Tel] = new G4LogicalVolume(Silicon_12[Blkcor][Tel] ,Silicon_mat ,silicon_12[Blkcor][Tel]);

		PhysicalDetector12[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta1))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector12[Blkcor][Tel],
			silicon_12[Blkcor][Tel],
			logicWorld,
			false,
			silicon_12Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal1[Blkcor][Tel] = new G4Trd(csicrystal1[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector13[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal1[Blkcor][Tel] ,scintillator_mat ,csicrystal1[Blkcor][Tel]);

		PhysicalDetector13[Blkcor][Tel] = new G4PVPlacement(Rot_1Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta1))*TMath::Cos((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta1))*TMath::Sin((Blkphi1+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta1))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector13[Blkcor][Tel],
			csicrystal1[Blkcor][Tel],
			logicWorld,
			false,
			csicrysial1Num[Blkcor][Tel],
			true);
	  }
	  {//2
		Silicon_21[Blkcor][Tel] = new G4Box(silicon_21[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector21[Blkcor][Tel] = new G4LogicalVolume(Silicon_21[Blkcor][Tel],Silicon_mat,silicon_21[Blkcor][Tel]);
		PhysicalDetector21[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta2))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector21[Blkcor][Tel],//theta == 0.0316
			silicon_21[Blkcor][Tel],
			logicWorld,
			false,
			silicon_21Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_22[Blkcor][Tel] = new G4Box(silicon_22[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector22[Blkcor][Tel] = new G4LogicalVolume(Silicon_22[Blkcor][Tel] ,Silicon_mat ,silicon_22[Blkcor][Tel]);

		PhysicalDetector22[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta2))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector22[Blkcor][Tel],
			silicon_22[Blkcor][Tel],
			logicWorld,
			false,
			silicon_22Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal2[Blkcor][Tel] = new G4Trd(csicrystal2[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector23[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal2[Blkcor][Tel] ,scintillator_mat ,csicrystal2[Blkcor][Tel]);

		PhysicalDetector23[Blkcor][Tel] = new G4PVPlacement(Rot_2Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta2))*TMath::Cos(Blkphi2+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta2))*TMath::Sin(Blkphi2+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta2))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector23[Blkcor][Tel],
			csicrystal2[Blkcor][Tel],
			logicWorld,
			false,
			csicrysial2Num[Blkcor][Tel],
			true);
	  }
	  {//3
		Silicon_31[Blkcor][Tel] = new G4Box(silicon_31[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector31[Blkcor][Tel] = new G4LogicalVolume(Silicon_31[Blkcor][Tel],Silicon_mat,silicon_31[Blkcor][Tel]);
		PhysicalDetector31[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(pi*Blkcor/2))*mm, 1000*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(pi*Blkcor/2))*mm, 1000*(TMath::Cos(Blktheta3))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector31[Blkcor][Tel],//theta == 0.0147
			silicon_31[Blkcor][Tel],
			logicWorld,
			false,
			silicon_31Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_32[Blkcor][Tel] = new G4Box(silicon_32[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector32[Blkcor][Tel] = new G4LogicalVolume(Silicon_32[Blkcor][Tel] ,Silicon_mat ,silicon_32[Blkcor][Tel]);

		PhysicalDetector32[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(pi*Blkcor/2))*mm, 1000.625*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(pi*Blkcor/2))*mm, 1000.625*(TMath::Cos(Blktheta3))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector32[Blkcor][Tel],
			silicon_32[Blkcor][Tel],
			logicWorld,
			false,
			silicon_32Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal3[Blkcor][Tel] = new G4Trd(csicrystal3[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector33[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal3[Blkcor][Tel] ,scintillator_mat ,csicrystal3[Blkcor][Tel]);

		PhysicalDetector33[Blkcor][Tel] = new G4PVPlacement(Rot_3Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta3))*TMath::Cos(Blkphi3+(pi*Blkcor/2))*mm, 1051.058*(TMath::Sin(Blktheta3))*TMath::Sin(Blkphi3+(pi*Blkcor/2))*mm, 1051.058*(TMath::Cos(Blktheta3))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector33[Blkcor][Tel],
			csicrystal3[Blkcor][Tel],
			logicWorld,
			false,
			csicrysial3Num[Blkcor][Tel],
			true);
	  }
	  {//4
		Silicon_41[Blkcor][Tel] = new G4Box(silicon_41[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.3*mm)/2);
		logicalDetector41[Blkcor][Tel] = new G4LogicalVolume(Silicon_41[Blkcor][Tel],Silicon_mat,silicon_41[Blkcor][Tel]);
		PhysicalDetector41[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000*(TMath::Cos(Blktheta4))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector41[Blkcor][Tel],//theta == 0.0147
			silicon_41[Blkcor][Tel],
			logicWorld,
			false,
			silicon_41Num[Blkcor][Tel],
			true);

		//Si2		 
		Silicon_42[Blkcor][Tel] = new G4Box(silicon_42[Blkcor][Tel],(20*mm)/2,(20*mm)/2,(0.5*mm)/2);
		logicalDetector42[Blkcor][Tel] = new G4LogicalVolume(Silicon_42[Blkcor][Tel] ,Silicon_mat ,silicon_42[Blkcor][Tel]);

		PhysicalDetector42[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1000.625*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1000.625*(TMath::Cos(Blktheta4))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector42[Blkcor][Tel],
			silicon_42[Blkcor][Tel],
			logicWorld,
			false,
			silicon_42Num[Blkcor][Tel],
			true);

		//CsI crystal		  
		CsIcrystal4[Blkcor][Tel] = new G4Trd(csicrystal4[Blkcor][Tel],(20*mm)/2,(22.70*mm)/2,(20*mm)/2,(22.70*mm)/2,(100.01*mm)/2);
		logicalDetector43[Blkcor][Tel] = new G4LogicalVolume(CsIcrystal4[Blkcor][Tel] ,scintillator_mat ,csicrystal4[Blkcor][Tel]);

		PhysicalDetector43[Blkcor][Tel] = new G4PVPlacement(Rot_4Si[Blkcor][Tel],
			G4ThreeVector(1051.058*(TMath::Sin(Blktheta4))*TMath::Cos((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Sin(Blktheta4))*TMath::Sin((Blkphi4+(pi*Blkcor/2)/1.0))*mm, 1051.058*(TMath::Cos(Blktheta4))*mm).rotate((90*(Tel))*degree,G4ThreeVector(TMath::Sin(0.069886)*TMath::Cos(1.428893+(pi*Blkcor/2)),TMath::Sin(0.069886)*TMath::Sin(1.428893+(pi*Blkcor/2)),TMath::Cos(0.069886))),
			logicalDetector43[Blkcor][Tel],
			csicrystal4[Blkcor][Tel],
			logicWorld,
			false,
			csicrysial4Num[Blkcor][Tel],
			true);
	  }
	  //-------------------------------------------------------------------------------------------------------------------------------		
	  opTubeSurface[Blkcor][Tel] = new G4OpticalSurface(TubeSurface[Blkcor][Tel]);
	  opTubeSurface[Blkcor][Tel]->SetType(dielectric_metal);
	  opTubeSurface[Blkcor][Tel]->SetFinish(polished);
	  opTubeSurface[Blkcor][Tel]->SetModel(unified);

	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector13[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector23[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector33[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
	  new G4LogicalBorderSurface(TubeSurface[Blkcor][Tel],PhysicalDetector43[Blkcor][Tel],physWorld,opTubeSurface[Blkcor][Tel]);
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

	  opTubeSurface[Blkcor][Tel]->SetMaterialPropertiesTable(myST1);

	  Silicon_11VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_11VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector11[Blkcor][Tel] -> SetVisAttributes (Silicon_11VisAtt[Blkcor][Tel]);

	  Silicon_12VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_12VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector12[Blkcor][Tel] -> SetVisAttributes (Silicon_12VisAtt[Blkcor][Tel]);

	  CsIcrystal1VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal1VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector13[Blkcor][Tel] -> SetVisAttributes (CsIcrystal1VisAtt[Blkcor][Tel]);
	  //
	  Silicon_21VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_21VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector21[Blkcor][Tel] -> SetVisAttributes (Silicon_21VisAtt[Blkcor][Tel]);

	  Silicon_22VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_22VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector22[Blkcor][Tel] -> SetVisAttributes (Silicon_22VisAtt[Blkcor][Tel]);

	  CsIcrystal2VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal2VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector23[Blkcor][Tel] -> SetVisAttributes (CsIcrystal2VisAtt[Blkcor][Tel]);
	  //
	  Silicon_31VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_31VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector31[Blkcor][Tel] -> SetVisAttributes (Silicon_31VisAtt[Blkcor][Tel]);

	  Silicon_32VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_32VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector32[Blkcor][Tel] -> SetVisAttributes (Silicon_32VisAtt[Blkcor][Tel]);

	  CsIcrystal3VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal3VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector33[Blkcor][Tel] -> SetVisAttributes (CsIcrystal3VisAtt[Blkcor][Tel]);
	  //		
	  Silicon_41VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
	  Silicon_41VisAtt[Blkcor][Tel]-> SetForceWireframe(true);
	  logicalDetector41[Blkcor][Tel] -> SetVisAttributes (Silicon_41VisAtt[Blkcor][Tel]);

	  Silicon_42VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
	  Silicon_42VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector42[Blkcor][Tel] -> SetVisAttributes (Silicon_42VisAtt[Blkcor][Tel]);

	  CsIcrystal4VisAtt[Blkcor][Tel] = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
	  CsIcrystal4VisAtt[Blkcor][Tel] -> SetForceWireframe(true);
	  logicalDetector43[Blkcor][Tel] -> SetVisAttributes (CsIcrystal4VisAtt[Blkcor][Tel]);
	}
  }
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

  return physWorld;
}
