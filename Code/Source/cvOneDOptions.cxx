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

#include "cvOneDOptions.h"

namespace cvOneD{

// =====================
// PERFORM DATA CHECKING
// =====================

namespace{

  template <class T>
  int checkForDuplicateEntry(vector<T> vec){
    for(int loopA=0;loopA<vec.size();loopA++){
      for(int loopB=loopA+1;loopB<vec.size();loopB++){
        if(vec[loopA] == vec[loopB]){
          return loopA;
        }
      }
    }
    return -1;
  }

  template <typename T>
  void checkForDuplicateValues(const vector<T>& vec, const string& type) {
      if (int dblIDX = checkForDuplicateEntry(vec); dblIDX >= 0) {
          throw cvException(string("ERROR: Duplicate " + type + ": " + vec[dblIDX] + "\n").c_str());
      }
  }

  void checkForDuplicateValues(const cvLongVec& vec, const string& type) {
      if (int dblIDX = checkForDuplicateEntry(vec); dblIDX >= 0) {
          throw cvException(string("ERROR: Duplicate " + type + ": " + to_string(vec[dblIDX]) + "\n").c_str());
      }
  }


  int checkForPositiveVal(cvDoubleVec vec){
    for(int loopA=0;loopA<vec.size();loopA++){
      if(vec[loopA] < 0.0){
        return loopA;
      }
    }
    return -1;
  }

  template <typename T>
  void checkForNegativeValues(const vector<T>& vec, const vector<string>& names, const string& type) {
      if (int chkIDX = checkForPositiveVal(vec); chkIDX >= 0) {
          throw cvException(string("ERROR: Negative " + type + " in segment: " + names[chkIDX] + "\n").c_str());
      }
  }

  void checkForDuplicates(options const& opts){
    checkForDuplicateValues(opts.nodeName, "Node Name");
    checkForDuplicateValues(opts.jointName, "Joint Name");
    checkForDuplicateValues(opts.jointInletListNames, "JointInlet Name");
    checkForDuplicateValues(opts.jointOutletListNames, "JointOutlet Name");
    checkForDuplicateValues(opts.materialName, "Material Name");
    checkForDuplicateValues(opts.dataTableName, "Data Table Name");
    checkForDuplicateValues(opts.segmentName, "Segment Name");
    checkForDuplicateValues(opts.segmentID, "Segment ID");
  }

  void checkForNegatives(options const& opts){
    // Check for negative areas
    checkForNegativeValues(opts.segmentInInletArea, opts.segmentName, "Inlet area");
    checkForNegativeValues(opts.segmentInOutletArea, opts.segmentName, "Outlet area");
  }

  void checkSegmentLengthConsistency(options const& opts){
    bool inconsistencyFound = false;
    // Get total number of segments
    int totSegs = opts.segmentLength.size();
    double segLength = 0.0;
    int inNode = 0;
    int outNode = 0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double nodeDist = 0.0;
    for(int loopA=0;loopA<totSegs;loopA++){
      // Get Current Segment Length
      segLength = opts.segmentLength[loopA];
      // Get end nodes
      inNode = opts.segmentInNode[loopA];
      outNode = opts.segmentOutNode[loopA];
      // Get node spatial distance
      dx = opts.nodeXcoord[outNode] - opts.nodeXcoord[inNode];
      dy = opts.nodeYcoord[outNode] - opts.nodeYcoord[inNode];
      dz = opts.nodeZcoord[outNode] - opts.nodeZcoord[inNode];

      nodeDist = sqrt(dx*dx + dy*dy + dz*dz);
    }
    if(inconsistencyFound){
      printf("WARNING: Inconsistency detected between segment length and end node distance.\n");
      printf("Changing the segment lengths.\n");
    }
  }

  template <class T>
  bool checkContains(T value, vector<T> vec){
    for(int loopA=0;loopA<vec.size();loopA++){
      if(vec[loopA] == value){
        return true;
      }
    }
    return false;
  }

  void checkSegmentHasNodes(options const& opts){
    int inNode = 0;
    int outNode = 0;
    for(int loopA=0;loopA<opts.segmentName.size();loopA++){
      // Get end nodes
      inNode = opts.segmentInNode[loopA];
      outNode = opts.segmentOutNode[loopA];
      // Check
      if((inNode < 0)||(inNode >= opts.nodeName.size())){
        throw cvException(string("ERROR: Missing Node in Segment: " + opts.segmentName[loopA] + "\n").c_str());
      }
      if((outNode < 0)||(outNode >= opts.nodeName.size())){
        throw cvException(string("ERROR: Missing Node in Segment: " + opts.segmentName[loopA] + "\n").c_str());
      }
    }
  }

  void checkSegmentHasMaterials(options const& opts) {
    for (int segmentIndex = 0;
        segmentIndex < opts.segmentName.size();
        ++segmentIndex) {
      const auto& materialName = opts.segmentMatName[segmentIndex];

      if (!checkContains(materialName, opts.materialName)) {
        const auto errMsg =
            "ERROR: Segment '" + opts.segmentName[segmentIndex] +
            "' references material '" + materialName +
            "', but that material is not defined in the 1D input file.";
        throw cvException(errMsg.c_str());
      }
    }
  }

  void checkJointHasNodes(options const& opts){
    string currNodeName;
    for(int loopA=0;loopA<opts.jointName.size();loopA++){
      // Get end nodes
      currNodeName = opts.jointNode[loopA];
      // Check If
      if(!checkContains(currNodeName,opts.nodeName)){
        throw cvException(string("ERROR: Missing Node " + currNodeName + " in Joint: " + opts.jointName[loopA] + "\n").c_str());
      }
    }
  }

} // namespace

void validateOptions(options const& opts){

  // Check if there are duplicate options
  checkForDuplicates(opts);

  // Check for negative values when there shouldn't be
  checkForNegatives(opts);

  // Check the consistency between node coords and segment lengths
  checkSegmentLengthConsistency(opts);

  // Check if the segments refer to a node that is not there
  checkSegmentHasNodes(opts);

  // Check if the segments refer to a material that is not defined
  checkSegmentHasMaterials(opts);

  // Check if the joints refer to a node that is not there
  checkJointHasNodes(opts);

}

} // namespace cvOneD