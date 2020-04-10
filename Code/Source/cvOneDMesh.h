#ifndef CVONEDMESH_H
#define CVONEDMESH_H

//
//  Mesh.h:  Data structure to handle segment meshes
//

# include "cvOneDEnums.h"
# include "cvOneDVector.h"

// Have to define Node and Element structures
struct cvOneDMeshNode{
  long id; // Node ID

  double x; // Node Position;
  double y;
  double z;

  double Q; // Output Values;
  double A;       
};

struct cvOneDElement{
  long id; // Element ID
    
  long inNode; // inletNode;
  long outNode; // outletNode;
    
  double h; // Element size    
};

// Mesh simply maintains the Nodes and Elements;
struct cvOneDMesh{
  cvOneDVector<cvOneDMeshNode> NodesList;
  cvOneDVector<cvOneDElement>  ElementList;
};

#endif // CVONEDMESH_H

