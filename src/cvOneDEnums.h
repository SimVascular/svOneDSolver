#ifndef CVONEDENUMS_H
#define CVONEDENUMS_H

//
//  cvOneDEnums.h: Enumerated types and typedefs for me234c Project.
//  ~~~~~~~
//
//  NOTES:
//
//  In order to aviid name conflicts, both within the API, and with external
//  symbols, the namespace for enumerated types has been partitioned with 
//  structs (compare "Effective C++ by Scott Meyers).  For example, 
//  project-specific NO_ERROR is accessed as ErrorScope::NO_ERROR in order 
//  not to clash with other globally scoped symbols (e.g., in 
//  /usr/include/sys/cachectl.h).  While a little verbose, this is consistent 
//  and safe approach that also enforces strong typing.  For example:
//
//    ErrorType var = ErrorTypeScope::NO_ERROR;
//
//  The burden of fully qualifying the scope of the enumerated types can
//  be remedied in situations where name conflicts are not a problem by
//  introducing handy shortcuts.  For example, if application source file
//  does not include /usr/include/sys/cachectl.h or any other packages that
//  define a globbaly scoped NO_ERROR, one can define the following:
//
//    static ErrorType NO_ERROR = ErrorTypeScope::NO_ERROR;
//
//  or
//
//    #define NO_ERROR ErrorTypeScope::NO_ERROR
//
//  In either case, applications can then use:
//
//    ErrorType var = NO_ERROR;
//
//  Note, that the first approach is more inducive towards strong typing.
//
//

#include <iostream>

using namespace std;

// Error Handling Type
struct ErrorTypeScope{
  enum ErrorType{
    NO_ERROR      = 0,
    BAD_VALUE     = 1,
    BAD_ENUM      = 2,
    OUT_OF_MEMORY = 3,
    UNSUPPORTED   = 4,
    UNKNOWN       = 5,
    OTHER_ERROR   = 6
  };
};

ostream & operator << (ostream &, ErrorTypeScope::ErrorType &);
typedef ErrorTypeScope::ErrorType ErrorType;  

// Mesh Topology Type
struct MeshTypeScope{
  enum MeshType{
    UNIFORM       = 0,
    FROM_DENSE    = 1,
    TO_DENSE      = 2,
    USER_DEFINED  = 3
  };
};

ostream & operator << (ostream &, MeshTypeScope::MeshType &);
typedef MeshTypeScope::MeshType MeshType;

// Boundary Condition Type
struct BoundCondTypeScope{
  enum BoundCondType{
    PRESSURE        = 0,
    FLOW            = 2,
    RESISTANCE      = 3, // resistance with full pressure
    RESISTANCE_TIME = 4, // full pressure
    PRESSURE_WAVE   = 5, // pressure wave in time
    RCR             = 6, // RCR or Windkessel model, without assuming periodicity
    CORONARY        = 9, // Coronary BC, Jongmin Seo & Hyunjin Kim
    NOBOUND 
  };
};

ostream & operator << (ostream &, BoundCondTypeScope::BoundCondType &);
typedef BoundCondTypeScope::BoundCondType BoundCondType;

//Material property
struct MaterialTypeScope{
  enum MaterialType{
    MATERIAL_UNKNOWN = 0,
    MATERIAL_OLUFSEN = 1,
    MATERIAL_LINEAR  = 2,
  };
};

#define MaterialType_MATERIAL_UNKNOWN 0
#define MaterialType_MATERIAL_OLUFSEN 1
#define MaterialType_MATERIAL_LINEAR  2

typedef int MaterialType;
 
ostream & operator << (ostream&, MaterialType &);

// Losses    
struct MinorLossScope{
  enum MinorLoss{
    NONE                      = 0,
    STENOSIS                  = 1,
    BRANCH_THROUGH_DIVIDING   = 2,
    BRANCH_SIDE_DIVIDING      = 3,
    BRANCH_THROUGH_CONVERGING = 4,
    BRANCH_SIDE_CONVERGING    = 5,
    BIFURCATION_BRANCH        = 6,
    //DOUBLE_BRANCH_CONVERGE  = 7,
    //BEND_90     = ,
    //BEND_45     = ,
    //BEND_180   = ,
  };
};

ostream & operator << (ostream&, MinorLossScope::MinorLoss &);
typedef MinorLossScope::MinorLoss MinorLoss;

// We need a unique identifier for each math model type.
// However each math model type could have several implementations.
// The calling order of these math models is defined by the MthModelOrder

// Math Model Type

struct MthModelTypeScope{
  enum MthModelType{
    MthSingleSegment_Model_Type,
    MthBranch_Model_Type,
    MthStenosis_Model_Type,
    MthBad_Model_Type
  };
};
typedef MthModelTypeScope::MthModelType MthModelType;

//currently each model has only one implementation
struct MthModelIMPScope{
  enum MthModelIMP{
    MthSingleSegment_Model_IMP,
    MthBranch_Model_IMP,
    MthStenosis_Model_IMP
  };
};

typedef MthModelIMPScope::MthModelIMP MthModelIMP;

struct MthModelOrderScope{
  enum MthModelOrder{
    MthSingleSegment_Model_Order = 1000,
    MthBranch_Model_Order = 2000,
    MthStenosis_Model_Order = 3000
  };
};
typedef MthModelOrderScope::MthModelOrder MthModelOrder;

enum typeOfEquation {linear, nonlinear};

// Output Type
struct OutputTypeScope {
  enum OutputType {
    OUTPUT_TEXT = 0,
    OUTPUT_VTK  = 1,
    OUTPUT_BOTH = 2
  };
};


#endif // CVONEDENUMS_H
