#include "DC.hh"
#include "detector.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4MaterialPropertiesTable.hh"

#include <cmath>

DC::DC(G4double density, G4double collimatorLength)
    : G4VUserDetectorConstruction(),
      fDensity(density),
      fCollimatorLength(collimatorLength),
      fGasMaterial(nullptr),
      fTeflon(nullptr),
      fHDPE(nullptr),
      fSilicon(nullptr),
      fAluminium(nullptr),
      fMylar(nullptr),
      fVacuum(nullptr),
      fLead(nullptr),
      fPolyethylene(nullptr),
      fBoratedPolyethylene(nullptr),
      fSipmSurf(nullptr),
      fMylarSurf(nullptr),
      fTeflonSurf(nullptr),
      fCf4SiSurface(nullptr),
      fWorldPhys(nullptr),
      fWorldLog(nullptr),
      fScoringVolume(nullptr),
      fCheckOverlaps(true)
{
}

DC::~DC()
{
    // Clean up optical surfaces
    if (fSipmSurf) delete fSipmSurf;
    if (fMylarSurf) delete fMylarSurf;
    if (fTeflonSurf) delete fTeflonSurf;
    if (fCf4SiSurface) delete fCf4SiSurface;
}

void DC::DefineMaterials()
{
    G4int ncomponents, natoms;
    
    // G4 materials 
    auto nistManager = G4NistManager::Instance();
    
    // Elements
    nistManager->FindOrBuildElement("C");
    nistManager->FindOrBuildElement("F");
    nistManager->FindOrBuildElement("Ar");
    nistManager->FindOrBuildElement("H");
    nistManager->FindOrBuildElement("O");
    nistManager->FindOrBuildElement("Si");
    nistManager->FindOrBuildElement("Al");
    nistManager->FindOrBuildElement("Pb");
    nistManager->FindOrBuildElement("B");
    
    auto C = G4Element::GetElement("C");
    auto F = G4Element::GetElement("F");
    auto H = G4Element::GetElement("H");
    auto O = G4Element::GetElement("O");
    auto Si = G4Element::GetElement("Si");
    auto Al = G4Element::GetElement("Al");
    auto Pb = G4Element::GetElement("Pb");
    auto Ar = G4Element::GetElement("Ar");
    auto B = G4Element::GetElement("B");
    
    // Vacuum
    G4double Vdens = 1.e-25 * g / cm3;
    G4double Vpres = 1.e-19 * pascal;
    G4double Vtemp = 0.1 * kelvin;
    fVacuum = new G4Material("Vacuum", 1, 1.01 * g / mole, Vdens, kStateGas, Vtemp, Vpres);
    
    // CF4
    CF4 = new G4Material("CF4", fDensity * mg / cm3, 2, kStateGas);
    CF4->AddElement(C, 1);
    CF4->AddElement(F, 4);

     const G4int iNbEntries = 300;
  std::vector<G4double>CF4PhotonMomentum  = {6.2*eV,6.138613861*eV,6.078431373*eV,6.019417476*eV,5.961538462*eV,5.904761905*eV,5.849056604*eV,5.794392523*eV,5.740740741*eV,5.688073394*eV,5.636363636*eV,5.585585586*eV,
						5.535714286*eV,5.486725664*eV,5.438596491*eV,5.391304348*eV,5.344827586*eV,5.299145299*eV,5.254237288*eV,5.210084034*eV,5.166666667*eV,5.123966942*eV,5.081967213*eV,5.040650407*eV,5*eV,4.96*eV,4.920634921*eV,
						4.881889764*eV,4.84375*eV,4.80620155*eV,4.769230769*eV,4.732824427*eV,4.696969697*eV,4.661654135*eV,4.626865672*eV,4.592592593*eV,4.558823529*eV,4.525547445*eV,4.492753623*eV,4.460431655*eV,
						4.428571429*eV,4.397163121*eV,4.366197183*eV,4.335664336*eV,4.305555556*eV,4.275862069*eV,4.246575342*eV,4.217687075*eV,4.189189189*eV,4.161073826*eV,4.133333333*eV,4.105960265*eV,4.078947368*eV,4.052287582*eV,
						4.025974026*eV,4*eV,3.974358974*eV,3.949044586*eV,3.924050633*eV,3.899371069*eV,3.875*eV,3.850931677*eV,3.827160494*eV,3.803680982*eV,3.780487805*eV,3.757575758*eV,3.734939759*eV,
						3.71257485*eV,3.69047619*eV,3.668639053*eV,3.647058824*eV,3.625730994*eV,3.604651163*eV,3.583815029*eV,3.563218391*eV,3.542857143*eV,3.522727273*eV,3.502824859*eV,3.483146067*eV,
						3.463687151*eV,3.444444444*eV,3.425414365*eV,3.406593407*eV,3.387978142*eV,3.369565217*eV,3.351351351*eV,3.333333333*eV,3.315508021*eV,3.29787234*eV,3.28042328*eV,3.263157895*eV,
						3.246073298*eV,3.229166667*eV,3.212435233*eV,3.195876289*eV,3.179487179*eV,3.163265306*eV,3.147208122*eV,3.131313131*eV,3.115577889*eV,3.1*eV,3.084577114*eV,3.069306931*eV,3.054187192*eV,
						3.039215686*eV,3.024390244*eV,3.009708738*eV,2.995169082*eV,2.980769231*eV,2.966507177*eV,2.952380952*eV,2.938388626*eV,2.924528302*eV,2.910798122*eV,2.897196262*eV,2.88372093*eV,
						2.87037037*eV,2.857142857*eV,2.844036697*eV,2.831050228*eV,2.818181818*eV,2.805429864*eV,2.792792793*eV,2.780269058*eV,2.767857143*eV,2.755555556*eV,2.743362832*eV,
						2.731277533*eV,2.719298246*eV,2.707423581*eV,2.695652174*eV,2.683982684*eV,2.672413793*eV,2.660944206*eV,2.64957265*eV,2.638297872*eV,2.627118644*eV,2.616033755*eV,2.605042017*eV,
						2.594142259*eV,2.583333333*eV,2.572614108*eV,2.561983471*eV,2.551440329*eV,2.540983607*eV,2.530612245*eV,2.520325203*eV,2.510121457*eV,2.5*eV,2.489959839*eV,2.48*eV,2.470119522*eV,2.46031746*eV,
						2.450592885*eV,2.440944882*eV,2.431372549*eV,2.421875*eV,2.412451362*eV,2.403100775*eV,2.393822394*eV,2.384615385*eV,2.375478927*eV,2.366412214*eV,2.357414449*eV,
						2.348484848*eV,2.339622642*eV,2.330827068*eV,2.322097378*eV,2.313432836*eV,2.304832714*eV,2.296296296*eV,2.287822878*eV,2.279411765*eV,2.271062271*eV,2.262773723*eV,2.254545455*eV,2.246376812*eV,2.238267148*eV,
						2.230215827*eV,2.222222222*eV,2.214285714*eV,2.206405694*eV,2.19858156*eV,2.190812721*eV,2.183098592*eV,2.175438596*eV,2.167832168*eV,2.160278746*eV,2.152777778*eV,2.14532872*eV,2.137931034*eV,2.130584192*eV,
						2.123287671*eV,2.116040956*eV,2.108843537*eV,2.101694915*eV,2.094594595*eV,2.087542088*eV,2.080536913*eV,2.073578595*eV,2.066666667*eV,2.059800664*eV,2.052980132*eV,2.04620462*eV,2.039473684*eV,2.032786885*eV,
						2.026143791*eV,2.019543974*eV,2.012987013*eV,2.006472492*eV,2*eV,1.993569132*eV,1.987179487*eV,1.980830671*eV,1.974522293*eV,1.968253968*eV,1.962025316*eV,1.955835962*eV,
						1.949685535*eV,1.943573668*eV,1.9375*eV,1.931464174*eV,1.925465839*eV,1.919504644*eV,1.913580247*eV,1.907692308*eV,1.901840491*eV,1.896024465*eV,1.890243902*eV,1.88449848*eV,1.878787879*eV,1.873111782*eV,
						1.86746988*eV,1.861861862*eV,1.856287425*eV,1.850746269*eV,1.845238095*eV,1.839762611*eV,1.834319527*eV,1.828908555*eV,1.823529412*eV,1.818181818*eV,
						1.812865497*eV,1.807580175*eV,1.802325581*eV,1.797101449*eV,1.791907514*eV,1.786743516*eV,1.781609195*eV,1.776504298*eV,1.771428571*eV,1.766381766*eV,1.761363636*eV,1.756373938*eV,1.751412429*eV,
						1.746478873*eV,1.741573034*eV,1.736694678*eV,1.731843575*eV,1.727019499*eV,1.722222222*eV,1.717451524*eV,1.712707182*eV,1.707988981*eV,1.703296703*eV,1.698630137*eV,1.693989071*eV,1.689373297*eV,1.684782609*eV,
						1.680216802*eV,1.675675676*eV,1.67115903*eV,1.666666667*eV,1.662198391*eV,1.657754011*eV,1.653333333*eV,1.64893617*eV,1.644562334*eV,1.64021164*eV,1.635883905*eV,1.631578947*eV,1.627296588*eV,1.623036649*eV,
						1.618798956*eV,1.614583333*eV,1.61038961*eV,1.606217617*eV,1.602067183*eV,1.597938144*eV,1.593830334*eV,1.58974359*eV,1.585677749*eV,1.581632653*eV,1.577608142*eV,1.573604061*eV,1.569620253*eV,
						1.565656566*eV,1.561712846*eV,1.557788945*eV,1.553884712*eV};
    // Sort the array in ascending order
  std::sort(CF4PhotonMomentum.begin(), CF4PhotonMomentum.end());

  std::vector<G4double>CF4Scintillation_Fast     = {0.0029,0.0029,0.0017,0.0024,0.0018,0.0011,0.0027,0.0009,0.0003,0.0019,0.0030,0.0024,0.0023,0.0036,0.0039,0.0056,
						0.0049,0.0061,0.0053,0.0052,0.0056,0.0064,0.0072,0.0064,0.0080,0.0071,0.0056,0.0069,0.0053,0.0070,0.0060,0.0057,0.0071,0.0066,0.0066,
						0.0055,0.0082,0.0076,0.0093,0.0089,0.0106,0.0109,0.0105,0.0102,0.0120,0.0121,0.0102,0.0097,0.0120,0.0126,0.0097,0.0103,0.0097,0.0084,
						0.0119,0.0112,0.0096,0.0171,0.0235,0.0078,0.0089,0.0071,0.0065,0.0074,0.0073,0.0074,0.0074,0.0080,0.0143,0.0522,0.0069,0.0076,0.0042,
						0.0059,0.0039,0.0053,0.0054,0.0185,0.0077,0.0599,0.0048,0.0034,0.0041,0.0041,0.0047,0.0059,0.0046,0.0065,0.0128,0.0037,0.0167,0.0053,
						0.0038,0.0042,0.0046,0.0032,0.0037,0.0073,0.0049,0.0067,0.0116,0.0054,0.0077,0.0111,0.0042,0.0043,0.0037,0.0046,0.0041,0.0028,0.0055,
						0.0031,0.0048,0.0057,0.0056,0.0035,0.0039,0.0068,0.0051,0.0037,0.0054,0.0048,0.0061,0.0033,0.0050,0.0052,0.0047,0.0014,0.0043,0.0041,
						0.0023,0.0062,0.0036,0.0038,0.0039,0.0043,0.0049,0.0049,0.0036,0.0048,0.0039,0.0023,0.0035,0.0025,0.0036,0.0010,0.0044,0.0013,0.0041,
						0.0021,0.0016,0.0046,0.0040,0.0034,0.0027,0.0026,0.0034,0.0004,0.0037,0.0004,0.0036,0.0029,0.0029,0.0036,0.0055,0.0034,0.0034,0.0025,
						0.0028,0.0055,0.0064,0.0037,0.0029,0.0047,0.0058,0.0040,0.0062,0.0055,0.0029,0.0067,0.0070,0.0080,0.0060,0.0094,0.0082,0.0072,0.0089,
						0.0117,0.0102,0.0134,0.0131,0.0131,0.0120,0.0135,0.0096,0.0107,0.0179,0.0210,0.0172,0.0165,0.0167,0.0176,0.0137,0.0196,0.0217,0.0175,
						0.0223,0.0192,0.0222,0.0188,0.0184,0.0183,0.0156,0.0098,0.0198,0.0268,0.0188,0.0236,0.0208,0.0171,0.0229,0.0228,0.0227,0.0204,0.0184,
						0.0190,0.0185,0.0145,0.0138,0.0122,0.0180,0.0132,0.0146,0.0087,0.0039,0.0147,0.0000,0.0000,0.0137,0.0084,0.0094,0.0114,0.0078,0.0100,
						0.0069,0.0055,0.0164,0.0113,0.0148,0.0053,0.0054,0.0065,0.0092,0.0000,0.0047,0.0000,0.0071,0.0000,0.0057,0.0063,0.0064,0.0050,0.0077,
						0.0034,0.0025,0.0000,0.0041,0.0025,0.0019,0.0042,0.0030,0.0000,0.0030,0.0000,0.0000,0.0000,0.0027,0.0000,0.0000,0.0000,0.0000,0.0006,
						0.0051,0.0083,0.0000,0.0000,0.0064,0.0003,0.0002,0.0074,0.0038,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};
   std::vector<G4double>CF4Scintillation_Slow   = {0.0029,0.0029,0.0017,0.0024,0.0018,0.0011,0.0027,0.0009,0.0003,0.0019,0.0030,0.0024,0.0023,0.0036,0.0039,0.0056,
						0.0049,0.0061,0.0053,0.0052,0.0056,0.0064,0.0072,0.0064,0.0080,0.0071,0.0056,0.0069,0.0053,0.0070,0.0060,0.0057,0.0071,0.0066,0.0066,
						0.0055,0.0082,0.0076,0.0093,0.0089,0.0106,0.0109,0.0105,0.0102,0.0120,0.0121,0.0102,0.0097,0.0120,0.0126,0.0097,0.0103,0.0097,0.0084,
						0.0119,0.0112,0.0096,0.0171,0.0235,0.0078,0.0089,0.0071,0.0065,0.0074,0.0073,0.0074,0.0074,0.0080,0.0143,0.0522,0.0069,0.0076,0.0042,
						0.0059,0.0039,0.0053,0.0054,0.0185,0.0077,0.0599,0.0048,0.0034,0.0041,0.0041,0.0047,0.0059,0.0046,0.0065,0.0128,0.0037,0.0167,0.0053,
						0.0038,0.0042,0.0046,0.0032,0.0037,0.0073,0.0049,0.0067,0.0116,0.0054,0.0077,0.0111,0.0042,0.0043,0.0037,0.0046,0.0041,0.0028,0.0055,
						0.0031,0.0048,0.0057,0.0056,0.0035,0.0039,0.0068,0.0051,0.0037,0.0054,0.0048,0.0061,0.0033,0.0050,0.0052,0.0047,0.0014,0.0043,0.0041,
						0.0023,0.0062,0.0036,0.0038,0.0039,0.0043,0.0049,0.0049,0.0036,0.0048,0.0039,0.0023,0.0035,0.0025,0.0036,0.0010,0.0044,0.0013,0.0041,
						0.0021,0.0016,0.0046,0.0040,0.0034,0.0027,0.0026,0.0034,0.0004,0.0037,0.0004,0.0036,0.0029,0.0029,0.0036,0.0055,0.0034,0.0034,0.0025,
						0.0028,0.0055,0.0064,0.0037,0.0029,0.0047,0.0058,0.0040,0.0062,0.0055,0.0029,0.0067,0.0070,0.0080,0.0060,0.0094,0.0082,0.0072,0.0089,
						0.0117,0.0102,0.0134,0.0131,0.0131,0.0120,0.0135,0.0096,0.0107,0.0179,0.0210,0.0172,0.0165,0.0167,0.0176,0.0137,0.0196,0.0217,0.0175,
						0.0223,0.0192,0.0222,0.0188,0.0184,0.0183,0.0156,0.0098,0.0198,0.0268,0.0188,0.0236,0.0208,0.0171,0.0229,0.0228,0.0227,0.0204,0.0184,
						0.0190,0.0185,0.0145,0.0138,0.0122,0.0180,0.0132,0.0146,0.0087,0.0039,0.0147,0.0000,0.0000,0.0137,0.0084,0.0094,0.0114,0.0078,0.0100,
						0.0069,0.0055,0.0164,0.0113,0.0148,0.0053,0.0054,0.0065,0.0092,0.0000,0.0047,0.0000,0.0071,0.0000,0.0057,0.0063,0.0064,0.0050,0.0077,
						0.0034,0.0025,0.0000,0.0041,0.0025,0.0019,0.0042,0.0030,0.0000,0.0030,0.0000,0.0000,0.0000,0.0027,0.0000,0.0000,0.0000,0.0000,0.0006,
						0.0051,0.0083,0.0000,0.0000,0.0064,0.0003,0.0002,0.0074,0.0038,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,0.0000};


  const G4int iNbEntries_1 = 3;
  G4double CF4PhotonMomentum_1[iNbEntries_1] = {200*eV,500*eV,798*eV};
  G4double CF4RefractiveIndex[iNbEntries_1]  = {1.004,1.004,1.004};
  G4double CF4AbsorbtionLength[iNbEntries_1] = {100.*cm, 100.*cm, 100.*cm};
  G4double CF4ScatteringLength[iNbEntries_1] = {30.*cm,  30.*cm,  30.*cm};
  G4MaterialPropertiesTable *CF4PropertiesTable = new G4MaterialPropertiesTable();
  CF4PropertiesTable->AddProperty("SCINTILLATIONCOMPONENT1", CF4PhotonMomentum, CF4Scintillation_Fast, iNbEntries);
  CF4PropertiesTable->AddProperty("SCINTILLATIONCOMPONENT2", CF4PhotonMomentum, CF4Scintillation_Slow, iNbEntries);
  CF4PropertiesTable->AddProperty("RINDEX", CF4PhotonMomentum_1, CF4RefractiveIndex, iNbEntries_1);
  CF4PropertiesTable->AddProperty("ABSLENGTH", CF4PhotonMomentum_1, CF4AbsorbtionLength, iNbEntries_1);
  CF4PropertiesTable->AddProperty("RAYLEIGH", CF4PhotonMomentum_1, CF4ScatteringLength, iNbEntries_1);
  CF4PropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 2500./keV,true);  // for electron recoil
  CF4PropertiesTable->AddConstProperty("RESOLUTIONSCALE", 1.0);
  CF4PropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 3.*ns,true);
  CF4PropertiesTable->AddConstProperty("SCINTILLATIONTIMECONSTANT2", 10.*ns,true);
  CF4PropertiesTable->AddConstProperty("YIELDRATIO", 1.0,true);
  CF4->SetMaterialPropertiesTable(CF4PropertiesTable);    
    // Ar:CF4 (90:10)
    G4Material* Ar_gas = new G4Material("Ar_gas", 1.782 * mg / cm3, 1, kStateGas);
    Ar_gas->AddElement(Ar, 1);
    
    fGasMaterial = new G4Material("Ar_CF4", 0.061 * mg / cm3, 2, kStateGas);
    fGasMaterial->AddMaterial(Ar_gas, 0.9);
    fGasMaterial->AddMaterial(CF4, 0.1);
    
    // Teflon
    fTeflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
    fTeflon->AddElement(C, 0.240183);
    fTeflon->AddElement(F, 0.759817);
    
    const G4int TNbEntries = 4;
    G4double pdTeflonPhotonMomentum[TNbEntries]  = {1.0*eV,6.91*eV, 6.98*eV, 7.05*eV};
    G4double pdTeflonRefractiveIndex[TNbEntries] = {1.34,1.34,1.34,1.34};
    G4double pdTeflonReflectivity[TNbEntries]    = {0.9,0.9,0.9,0.9};
    G4double pdTeflonAbsLength[TNbEntries]       = {0.2*mm,0.2*mm,0.2*mm,0.2*mm};
    G4double pdTeflonScatteringLength[TNbEntries] = {30*cm,30*cm,30*cm,30*cm};
    
    G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();
    pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, TNbEntries);
    pTeflonPropertiesTable->AddProperty("ABSLENGTH", pdTeflonPhotonMomentum, pdTeflonAbsLength, TNbEntries);
    pTeflonPropertiesTable->AddProperty("RAYLEIGH", pdTeflonPhotonMomentum, pdTeflonScatteringLength, TNbEntries);
    pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, TNbEntries);
    fTeflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);
    
    // HDPE (Polyethylene)
    fHDPE = new G4Material("HDPE", 0.93 * g/cm3, 2);
    fHDPE->AddElement(C, 2);
    fHDPE->AddElement(H, 4);
    fPolyethylene = fHDPE; // Keep reference for compatibility
    
    // Borated Polyethylene
    fBoratedPolyethylene = new G4Material("b_polyethylene", 0.94*g/cm3, 4, kStateSolid);
    fBoratedPolyethylene->AddElement(H, 11.6*perCent);
    fBoratedPolyethylene->AddElement(C, 61.2*perCent);
    fBoratedPolyethylene->AddElement(B, 5*perCent);
    fBoratedPolyethylene->AddElement(O, 22.2*perCent);
    
    // Silicon (SiPM)
    fSilicon = nistManager->FindOrBuildMaterial("G4_Si");
    
    // Aluminium
    fAluminium = nistManager->FindOrBuildMaterial("G4_Al");
    
    // Lead
    fLead = nistManager->FindOrBuildMaterial("G4_Pb");
    
    // Mylar
    fMylar = new G4Material("mylar", 1.39 * g/cm3, 3);
    fMylar->AddElement(C, 5);
    fMylar->AddElement(H, 4);
    fMylar->AddElement(O, 2);
    
    // Gas optical properties
    const G4int nScint = 18;
    G4double gasEnergy[nScint] = {
        1.55*eV, 1.61*eV, 1.70*eV, 1.78*eV, 1.88*eV,
        2.00*eV, 2.11*eV, 2.22*eV, 2.54*eV, 2.87*eV,
        3.42*eV, 3.89*eV, 4.35*eV, 4.72*eV, 4.88*eV,
        5.40*eV, 5.83*eV, 6.16*eV
    };
    
    G4double gasScintSp[nScint] = {
        0.01, 0.01, 0.04, 0.10, 0.20,
        0.26, 0.18, 0.04, 0.00, 0.02,
        0.04, 0.13, 0.30, 0.19, 0.25,
        0.03, 0.00, 0.00
    };
    
    G4double gasAbsLength[nScint];
    G4double gasScatLength[nScint];
    for (int i = 0; i < nScint; i++) {
        gasAbsLength[i] = 4.0 * m;
        gasScatLength[i] = 5.0 * cm;
    }
    
    auto gas_mpt = new G4MaterialPropertiesTable();
    gas_mpt->AddProperty("SCINTILLATIONCOMPONENT1", gasEnergy, gasScintSp, nScint);
    gas_mpt->AddProperty("RAYLEIGH", gasEnergy, gasScatLength, nScint);
    gas_mpt->AddProperty("ABSLENGTH", gasEnergy, gasAbsLength, nScint);
    gas_mpt->AddConstProperty("SCINTILLATIONYIELD", 5000/MeV);
    gas_mpt->AddConstProperty("SCINTILLATIONYIELD1", 1.0);
    gas_mpt->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 15*ns);
    gas_mpt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    fGasMaterial->SetMaterialPropertiesTable(gas_mpt);
}

G4VPhysicalVolume* DC::Construct()
{
    DefineMaterials();
    
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    
    // -------------------------- //
    // World volume - ArCF4 (90/10)
    // -------------------------- //
    G4double worldSize = 100 * cm;
    G4double worldDepth = 100 * cm;
    
    G4Box* worldBox = new G4Box("World", worldSize/2, worldSize/2, worldDepth/2);
    fWorldLog = new G4LogicalVolume(worldBox, CF4, "World");
    fWorldPhys = new G4PVPlacement(0, G4ThreeVector(), fWorldLog, "World", 0, false, 0, fCheckOverlaps);
    
    // Geometry parameters from the provided code
    G4double sipmSize = 1.0 * mm;
    G4double cellSize = sipmSize + 0.84 * mm;  // 1.84 mm
    G4double cellLength = fCollimatorLength;   // Use the passed parameter
    G4int nCells = 25;
    G4double collSize = nCells * cellSize;
    
    // Create collimator cell
    G4Box* cellBlock = new G4Box("cellBlock", cellLength/2, cellSize/2, cellSize/2);
    G4Box* cellHole = new G4Box("cellHole", cellLength/2, sipmSize/2, sipmSize/2);
    G4SubtractionSolid* cellSolid = new G4SubtractionSolid("cell", cellBlock, cellHole, 0, G4ThreeVector());
    
    G4LogicalVolume* cellLog = new G4LogicalVolume(cellSolid, fTeflon, "Cell");
    
    // Place collimator cells (4 sides)
    G4int copyNo = 0;
    for (G4int i = 0; i < 4; i++) {
        G4double angle = i * 90. * deg;
        auto cellRot = new G4RotationMatrix();
        cellRot->rotateZ(angle);
        
        for (G4int j = 0; j < nCells; j++) {
            G4double xPos, yPos;
            
            if (i == 0 || i == 2) { // Left and right sides
                xPos = (collSize/2 + cellLength/2) * std::cos(angle);
                yPos = ((-1.0 * nCells/2 + j + 0.5) * cellSize) + (collSize/2 + cellLength/2) * std::sin(angle);
            } else { // Top and bottom sides
                xPos = ((-1.0 * nCells/2 + j + 0.5) * cellSize) + (collSize/2 + cellLength/2) * std::cos(angle);
                yPos = (collSize/2 + cellLength/2) * std::sin(angle);
            }
            
            new G4PVPlacement(cellRot, G4ThreeVector(xPos, yPos, 0), cellLog, "Cell", fWorldLog, false, copyNo++, fCheckOverlaps);
        }
    }
    
    // Mylar electrodes (cathode and anode)
    G4double mylThickness = 0.012 * mm;
    G4double mylPos = cellSize;
    
    G4Box* mylBox = new G4Box("myl", collSize/2, collSize/2, mylThickness/2);
    G4LogicalVolume* mylLog = new G4LogicalVolume(mylBox, fAluminium, "Myl");
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, mylPos/2), mylLog, "MylA", fWorldLog, false, 0, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -mylPos/2), mylLog, "MylK", fWorldLog, false, 1, fCheckOverlaps);
    
    // HDPE converter (n,p converter)
    G4double convThickness = 0.01 * mm;
    G4double convPos = mylPos/2 + mylThickness/2 + convThickness/2;
    
    G4Box* convBox = new G4Box("conv", collSize/2, collSize/2, convThickness/2);
    G4LogicalVolume* convLog = new G4LogicalVolume(convBox, fHDPE, "Conv");
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, convPos), convLog, "Conv", fWorldLog, false, 0, fCheckOverlaps);
    
    // SiPMs - using your original volume names for compatibility
    G4double sipmWidth = 1.0 * mm;
    G4Box* sipmBox = new G4Box("sipm", sipmWidth/2, sipmSize/2, sipmSize/2);
    G4LogicalVolume* sipmLog = new G4LogicalVolume(sipmBox, fSilicon, "sensor_log1");
    
    // Create 4 different logical volumes with different names for your sensitive detector
    G4LogicalVolume* sipmLog1 = new G4LogicalVolume(sipmBox, fSilicon, "sensor_log1");
    G4LogicalVolume* sipmLog2 = new G4LogicalVolume(sipmBox, fSilicon, "sensor_log2");
    G4LogicalVolume* sipmLog3 = new G4LogicalVolume(sipmBox, fSilicon, "sensor_log3");
    G4LogicalVolume* sipmLog4 = new G4LogicalVolume(sipmBox, fSilicon, "sensor_log4");
    
    copyNo = 0;
    for (G4int i = 0; i < 4; i++) {
        G4double angle = i * 90. * deg;
        auto sipmRot = new G4RotationMatrix();
        sipmRot->rotateZ(angle);
        
        G4LogicalVolume* currentSipmLog;
        G4String volumeName;
        
        // Choose the right logical volume and name based on side
        if (i == 0) { // Left side
            currentSipmLog = sipmLog1;
            volumeName = "sensor_Vol1";
        } else if (i == 1) { // Top side
            currentSipmLog = sipmLog3;
            volumeName = "sensor_Vol3";
        } else if (i == 2) { // Right side
            currentSipmLog = sipmLog2;
            volumeName = "sensor_Vol2";
        } else { // Bottom side
            currentSipmLog = sipmLog4;
            volumeName = "sensor_Vol4";
        }
        
        for (G4int j = 0; j < nCells; j++) {
            G4double xPos, yPos;
            
            if (i == 0 || i == 2) { // Left or right
                xPos = (collSize/2 + cellLength + sipmWidth/2) * std::cos(angle);
                yPos = ((-1.0 * nCells/2 + j + 0.5) * cellSize) + (collSize/2 + cellLength/2) * std::sin(angle);
            } else { // Top or bottom
                xPos = ((-1.0 * nCells/2 + j + 0.5) * cellSize) + (collSize/2 + cellLength/2) * std::cos(angle);
                yPos = (collSize/2 + cellLength + sipmWidth/2) * std::sin(angle);
            }
            
            new G4PVPlacement(sipmRot, G4ThreeVector(xPos, yPos, 0), currentSipmLog, volumeName, fWorldLog, false, copyNo++, fCheckOverlaps);
        }
    }
    
    // Set the scoring volume (gas volume)
    fScoringVolume = fWorldLog;
    
    // Set visualization attributes
    cellLog->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 1.0, 0.0))); // Yellow
    mylLog->SetVisAttributes(G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.9))); // Grey
    convLog->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 1.0, 0.0))); // Green
    sipmLog1->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.0, 0.0))); // Red
    sipmLog2->SetVisAttributes(G4VisAttributes(G4Colour(0.0, 0.0, 1.0))); // Blue
    sipmLog3->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.5, 0.0))); // Orange
    sipmLog4->SetVisAttributes(G4VisAttributes(G4Colour(0.5, 0.0, 0.5))); // Purple
    
    fWorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    return fWorldPhys;
}

void DC::ConstructSDandField()
{
    // Set up sensitive detectors
    MySensitiveDetector* sensDet = new MySensitiveDetector("SensitiveDetector");
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    SDman->AddNewDetector(sensDet);
    
    // Set all SiPM logical volumes as sensitive
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
    
    for (auto lv : *lvStore) {
        G4String lvName = lv->GetName();
        if (lvName == "sensor_log1" || 
            lvName == "sensor_log2" || 
            lvName == "sensor_log3" || 
            lvName == "sensor_log4") {
            lv->SetSensitiveDetector(sensDet);
        }
    }
}
