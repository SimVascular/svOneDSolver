/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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
                              double profile_exponent, double pRef, 
                              double L_P, double P_ambient, double *params);

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
    
