#ifndef CVONEDMATERIALMANAGER_H
#define CVONEDMATERIALMANAGER_H

# include "cvOneDEnums.h"
# include "cvOneDMaterial.h"

#define MAX_NUM_OF_MATERIAL_INSTANCES 256

class cvOneDMaterialManager{
  public:
    cvOneDMaterialManager();
    ~cvOneDMaterialManager();
 
    int AddNewMaterial(MaterialType type, cvOneDMaterial* mat);

    int AddNewMaterialOlufsen(double density, double dynamicViscosity,
                              double profile_exponent, double pRef, double *params);

    int AddNewMaterialLinear(double density, double dynamicViscosity,
                             double profile_exponent, double pRef, double EHR);

    // caller must deallocate material instance to avoid memory leak
    cvOneDMaterial* GetNewInstance(int matID);

  private:
    int numMaterials;
    MaterialType types[MAX_NUM_OF_MATERIAL_INSTANCES];
    cvOneDMaterial* materials[MAX_NUM_OF_MATERIAL_INSTANCES];
};

#endif // CVONEDMATERIALMANAGER_H
    
