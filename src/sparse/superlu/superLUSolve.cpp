#include "superLUSolve.h"

# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>
# include <time.h>

# include "st_to_cc.h"

using namespace std;

#define SuperLU_MT

#if defined( SuperLU )

# include "slu_ddefs.h"

// ROUTINE TO USE superLU ON A SINGLE-PROCESSOR 
int superLUSolve(cvOneDKentry* Kentries,  double *b,int numEntries, int NNZ,int nunknown,double u[]) {


  int i,j;
  double *M  = (double*)malloc(NNZ*sizeof(double));
  int    *IA = (int*)malloc(NNZ*sizeof(int));
  int    *JA = (int*)malloc(NNZ*sizeof(int));
  j = 0;

  clock_t tstart_transfer;
  tstart_transfer=clock();

  for (i = 0; i < numEntries; i++){
    if(Kentries[i].col!=-1){
      IA[j] = Kentries[i].row;
      JA[j] = Kentries[i].col;
      M[j]  = Kentries[i].value;
      j++;
    }
  }

  SuperMatrix A;
  double *acc;
  double *b2;
  SuperMatrix B;
  int *ccc;
  int *icc;
  int info;
  SuperMatrix L;
  int nrhs = 1;
  int ncc;
  superlu_options_t options;
  int *perm_c;
  int permc_spec;
  int *perm_r;
  SuperLUStat_t stat;
  SuperMatrix U;

  //
  //  convert  sparse triplet (ST) format, into compressed column (CC) format.
  //

  ncc=NNZ;

  //
  //  Create the CC indices.
  //

  icc = new int[ncc];
  ccc = new int[nunknown+1];
  st_to_cc_index (NNZ, IA, JA, ncc, nunknown, icc, ccc);

  //
  //  Create the CC values.
  //

  acc = st_to_cc_values ( NNZ, IA, JA, M, ncc, nunknown, icc, ccc );

  free(M);
  free(IA);
  free(JA);

  //
  //  Convert the compressed column (CC) matrix into a SuperMatrix A.
  //

  dCreate_CompCol_Matrix ( &A, nunknown, nunknown, ncc, acc, icc, ccc, SLU_NC, SLU_D, SLU_GE );

  //
  //  Create Super Right Hand Side.
  //

  dCreate_Dense_Matrix ( &B, nunknown, nrhs, b, nunknown, SLU_DN, SLU_D, SLU_GE );

  //
  //  Set space for the permutations.
  // 
  
  perm_r = new int[nunknown];
  perm_c = new int[nunknown];

  //
  //  Set the input options.
  //
  
  set_default_options ( &options );
  options.ColPerm = NATURAL;
  
  //
  //  Initialize the statistics variables.
  //
  
  StatInit ( &stat );
  
  //
  //  Solve the linear system.
  //
  
  dgssv( &options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info );
  cout << "  dgssv done info = " << info << "\n";

  for ( i = 0; i < nunknown; i++ ){
    u[i] =b[i];
  }
  
  //
  //  Free memory.
  //
  
  free ( b );
  free ( b2 );
  free ( perm_c );
  free ( perm_r );

  delete [] acc;
  delete [] ccc;
  delete [] icc;

  Destroy_SuperMatrix_Store ( &A );
  Destroy_SuperMatrix_Store ( &B );
  Destroy_SuperNode_Matrix ( &L );
  Destroy_CompCol_Matrix ( &U );
  StatFree ( &stat );

  //
  //  Terminate.
  //

  return 0;
}

#elif defined(SuperLU_MT)

# include "slu_mt_ddefs.h"

extern "C" {
  int_t sp_ienv(int_t ispec);
  int_t dPrint_Dense_Matrix(SuperMatrix *A);
}

// ROUTINE TO USE superLU ON A MULTIPLE PROCESSORS ON A SHARED MEMORY MACHINE
int superLUSolve(cvOneDKentry* Kentries,  double *b,int numEntries, int NNZ,int nunknown,double u[]) {

  int i,j;
  double *M  = (double*)malloc(NNZ*sizeof(double));
  int    *IA = (int*)malloc(NNZ*sizeof(int));
  int    *JA = (int*)malloc(NNZ*sizeof(int));
  j = 0;

  clock_t tstart_transfer;
  tstart_transfer=clock();
  for (i = 0; i < numEntries; i++) {
    if(Kentries[i].col!=-1){
      IA[j] = Kentries[i].row;
      JA[j] = Kentries[i].col;
      M[j]  = Kentries[i].value;
      j++;
    }
  }

  int ncc;
  double *acc;
  int *ccc;
  int *icc;
  SuperMatrix   A;
  // NCformat *Astore;
  double   *a;
  int_t      *asub, *xa;
  int_t      *perm_r; // row permutations from partial pivoting
  int_t      *perm_c; // column permutation vector
  SuperMatrix   L;       // factor L
  SCPformat *Lstore;
  SuperMatrix   U;       // factor U
  NCPformat *Ustore;
  SuperMatrix   B;
  int_t    nrhs, ldx, info;
  int_t    nprocs; // maximum number of processors to use.
  int_t    panel_size, relax, maxsup;
  int_t    permc_spec;
  trans_t  trans;
  nrhs   = 1;
  trans  = NOTRANS;
  nprocs = 1;
  
  panel_size        = sp_ienv(1);
  relax             = sp_ienv(2);
  maxsup            = sp_ienv(3);

  ncc = NNZ;
 
  //
  //  Create the CC indices.
  //
  
  icc = new int[ncc];
  ccc = new int[nunknown+1];
  st_to_cc_index ( NNZ, IA, JA, ncc, nunknown, icc, ccc );

  //
  //  Create the CC values.
  //
  
  acc = st_to_cc_values ( NNZ, IA, JA, M, ncc, nunknown, icc, ccc );
  free(M);
  free(IA);
  free(JA);

  dCreate_CompCol_Matrix ( &A, nunknown, nunknown, NNZ, acc, icc, ccc, SLU_NC, SLU_D, SLU_GE );
  dCreate_Dense_Matrix ( &B, nunknown, nrhs, b, nunknown, SLU_DN, SLU_D, SLU_GE );

  if (!(perm_r = intMalloc(nunknown))) SUPERLU_ABORT("Malloc fails for perm_r[].");
  if (!(perm_c = intMalloc(nunknown))) SUPERLU_ABORT("Malloc fails for perm_c[].");

  //  * Get column permutation vector perm_c[], according to permc_spec:
  //  *   permc_spec = 0: natural ordering
  //  *   permc_spec = 1: minimum degree ordering on structure of A'*A
  //  *   permc_spec = 2: minimum degree ordering on structure of A'+A
  //  *   permc_spec = 3: approximate minimum degree for unsymmetric matrices

  permc_spec = 1;
  get_perm_c(permc_spec, &A, perm_c);

  pdgssv(nprocs, &A, perm_c, perm_r, &L, &U, &B, &info);

  for ( i = 0; i < nunknown; i++ ){
    u[i] = b[i];
  }

  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_SCP(&L);
  Destroy_CompCol_NCP(&U);

  return 0;
}

#endif
