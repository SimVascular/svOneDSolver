# include "cvOneDEnums.h"
# include "cvOneDGlobal.h"
# include "cvOneDMaterialManager.h"
# include "cvOneDMaterialOlufsen.h"
# include "cvOneDMaterialLinear.h"

// cvOneDMaterialManager* gMaterialManager = NULL;

// Constructor
cvOneDMaterialManager::cvOneDMaterialManager(){
  numMaterials = 0;
  for (int i=0;i<MAX_NUM_OF_MATERIAL_INSTANCES;i++){
    materials[i] = NULL;
  }
}

// Destructor
cvOneDMaterialManager::~cvOneDMaterialManager(){
}

int cvOneDMaterialManager::AddNewMaterial(MaterialType type, cvOneDMaterial* mat){
  int rtnID = numMaterials;
  types[rtnID] = type;
  materials[rtnID] = mat;
  numMaterials++;
  return rtnID;
}

int cvOneDMaterialManager::AddNewMaterialOlufsen(double density, double dynamicViscosity,
                                                 double profile_exponent, double pRef,
                                                 double *params){
  cvOneDMaterialOlufsen* olfmat = new cvOneDMaterialOlufsen();
  olfmat->SetDensity(density);
  olfmat->SetDynamicViscosity(dynamicViscosity);
  olfmat->SetProfileExponent(profile_exponent);
  olfmat->SetReferencePressure(pRef);
  olfmat->SetMaterialType(params,pRef);
  printf("new cvOneMaterialOlufsen called check pRef %f \n", olfmat->GetReferencePressure());
  return cvOneDGlobal::gMaterialManager->AddNewMaterial(MaterialType_MATERIAL_OLUFSEN,(cvOneDMaterial*)olfmat);
}

int cvOneDMaterialManager::AddNewMaterialLinear(double density, double dynamicViscosity,
                                                double profile_exponent, double pRef,
                                                double EHR){
  cvOneDMaterialLinear* linearmat = new cvOneDMaterialLinear();
  linearmat->SetDensity(density);
  linearmat->SetDynamicViscosity(dynamicViscosity);
  linearmat->SetProfileExponent(profile_exponent);
  linearmat->SetReferencePressure(pRef);
  linearmat->SetEHR(EHR,pRef);
  return cvOneDGlobal::gMaterialManager->AddNewMaterial(MaterialType_MATERIAL_LINEAR,(cvOneDMaterial*)linearmat);
}

// caller must deallocate material instance to avoid memory leak
cvOneDMaterial* cvOneDMaterialManager::GetNewInstance(int matID){
  if (types[matID] == MaterialType_MATERIAL_OLUFSEN) {
    cvOneDMaterialOlufsen* olfmat = new cvOneDMaterialOlufsen();
    printf("In GetNewInstance cvOneDMaterialOlufsen is called  matID=%i \n",matID);
    *olfmat = *((cvOneDMaterialOlufsen*)(materials[matID]));
    printf("In GetNewInstance cvOneDMaterialOlufsen* materials is called \n");
    return (cvOneDMaterial*)olfmat;
  }else if (types[matID] == MaterialType_MATERIAL_LINEAR) {
    cvOneDMaterialLinear* linearmat = new cvOneDMaterialLinear();
    *linearmat = *((cvOneDMaterialLinear*)(materials[matID]));
    return (cvOneDMaterial*)linearmat;
  }else{
    return NULL;
  }
}
