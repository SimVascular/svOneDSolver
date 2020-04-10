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

#include "cvOneDGlobal.h"

// GLOBAL FLAGS
bool cvOneDGlobal::isCreating = false;
bool cvOneDGlobal::isSolving = false;
int  cvOneDGlobal::outputType = 0; // Default Text Output
int  cvOneDGlobal::vtkOutputType = 0; // Default Multiple Files
int  cvOneDGlobal::CONSERVATION_FORM = 0;

// DEBUG MODE
bool cvOneDGlobal::debugMode = false;

// CURRENT MODEL INDEX
long cvOneDGlobal::currentModel = -1;

// VECTOR OF CREATED MODELS
vector<cvOneDModel*> cvOneDGlobal::gModelList;

// GLOBAL MATERIAL MANAGER OBJECT
cvOneDMaterialManager* cvOneDGlobal::gMaterialManager = NULL;

// GLOBAL Mth SEGMENT MODEL
cvOneDMthSegmentModel* cvOneDGlobal::gMthSegmentModel = NULL;

// GLOBAL SOLVER INSTANCE
cvOneDBFSolver* cvOneDGlobal::gBFSolver = NULL;

// GLOBAL VECTOR OF DATATABLES
vector<cvOneDDataTable*> cvOneDGlobal::gDataTables;

// GENERIC SOLVER INSTANCE
cvOneDLinearSolver* cvOneDGlobal::solver = NULL;
