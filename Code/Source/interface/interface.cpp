// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "interface.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <optional>
#include <algorithm>

#include "cvOneDGlobal.h"
#include "cvOneDModelManager.h"
#include "cvOneDOptions.h"
#include "cvOneDOptionsJsonParser.h"
#include "cvOneDOptionsJsonSerializer.h"
#include "cvOneDOptionsLegacySerializer.h"
#include "cvOneDUtility.h"

#ifndef NDEBUG
#define SVONED_DEBUG_LOG(message) \
  (std::cout << message << '\n')
#else
#define SVONED_DEBUG_LOG(message) ((void)0)
#endif

//////////////////////////////////////////////////////////
//      Static member data initialization               //
//////////////////////////////////////////////////////////

int OneDSolverInterface::problem_id_count_ = 0;
std::map<int, OneDSolverInterface*> OneDSolverInterface::interface_list_;

//-----------------------
// OneDSolverInterface
//-----------------------
/**
 * @brief Constructor for OneDSolverInterface.
 *
 * @param input_file_name The 1D input file name
 */
OneDSolverInterface::OneDSolverInterface(const std::string& input_file_name)
    : input_file_name_(input_file_name) {
  problem_id_ = problem_id_count_++;
  OneDSolverInterface::interface_list_[problem_id_] = this;
}

/**
 * @brief Destructor for OneDSolverInterface.
 */
OneDSolverInterface::~OneDSolverInterface() {}

namespace {

void createModelNodes(
    const cvOneD::options& opts,
    cvOneDModelManager* modelManager) {
    SVONED_DEBUG_LOG("Creating Nodes ... ");

    const int totalNodes = static_cast<int>(opts.nodeName.size());

    for (int nodeIndex = 0; nodeIndex < totalNodes; ++nodeIndex) {
    modelManager->CreateNode(
        opts.nodeName[nodeIndex].c_str(),
        opts.nodeXcoord[nodeIndex],
        opts.nodeYcoord[nodeIndex],
        opts.nodeZcoord[nodeIndex]);
    }
}

void createModelJoints(
    const cvOneD::options& opts,
    cvOneDModelManager* modelManager) {
  SVONED_DEBUG_LOG("Creating Joints ... ");

  const int totalJoints = static_cast<int>(opts.jointName.size());

  for (int jointIndex = 0; jointIndex < totalJoints; ++jointIndex) {
    const auto& inletName = opts.jointInletName[jointIndex];
    const auto& outletName = opts.jointOutletName[jointIndex];

    const int inletID =
        getListIDWithStringKey(inletName, opts.jointInletListNames);
    if (inletID < 0) {
        const auto errMsg =
            "ERROR: Cannot Find JOINTINLET for key " + inletName;
        throw cvException(errMsg.c_str());
    }

    const int outletID =
        getListIDWithStringKey(outletName, opts.jointOutletListNames);
    if (outletID < 0) {
        const auto errMsg =
            "ERROR: Cannot Find JOINTOUTLET for key " + outletName;
        throw cvException(errMsg.c_str());
    }

    const int totalInlets = opts.jointInletListNumber[inletID];
    const int totalOutlets = opts.jointOutletListNumber[outletID];

    std::vector<int> inletSegments(opts.jointInletList[inletID].begin(), opts.jointInletList[inletID].end());

    std::vector<int> outletSegments;
    outletSegments.reserve(totalOutlets);
    for (int i = 0; i < totalOutlets; ++i) {
      outletSegments.push_back(opts.jointOutletList[outletID][i]);
    }

    const auto& jointName = opts.jointName.at(jointIndex);
    const std::size_t nodeIndex = findJointNodeIndexOrThrow(
        opts.jointNode.at(jointIndex),
        opts.nodeName,
        jointName);
    modelManager->CreateJoint(
        jointName.c_str(),
        opts.nodeXcoord[nodeIndex],
        opts.nodeYcoord[nodeIndex],
        opts.nodeZcoord[nodeIndex],
        totalInlets,
        totalOutlets,
        inletSegments.empty() ? nullptr : inletSegments.data(),
        outletSegments.empty() ? nullptr : outletSegments.data());
    }
}

void createModelMaterials(
    const cvOneD::options& opts,
    cvOneDModelManager* modelManager) {
  SVONED_DEBUG_LOG("Creating Materials ... ");

  const int totalMaterials =
      static_cast<int>(opts.materialName.size());

  for (int materialIndex = 0;
       materialIndex < totalMaterials;
       ++materialIndex) {
    const bool isOlufsen =
        upper_string(opts.materialType[materialIndex]) == "OLUFSEN";

    const std::string materialType =
        isOlufsen ? "MATERIAL_OLUFSEN" : "MATERIAL_LINEAR";
    const int numberOfParameters = isOlufsen ? 3 : 1;

    double parameters[3] = {
        opts.materialParam1[materialIndex],
        opts.materialParam2[materialIndex],
        opts.materialParam3[materialIndex]
    };

    int materialID = 0;
    const int errorCode = modelManager->CreateMaterial(
        const_cast<char*>(opts.materialName[materialIndex].c_str()),
        const_cast<char*>(materialType.c_str()),
        opts.materialDensity[materialIndex],
        opts.materialViscosity[materialIndex],
        opts.materialExponent[materialIndex],
        opts.materialPRef[materialIndex],
        numberOfParameters,
        parameters,
        &materialID);

    if (errorCode == CV_ERROR) {
      const std::string errMsg =
          "ERROR: Error Creating MATERIAL " +
          std::to_string(materialIndex) + "\n";
      throw cvException(errMsg.c_str());
    }
  }
}

void createModelDataTables(
    const cvOneD::options& opts,
    cvOneDModelManager* modelManager) {
  SVONED_DEBUG_LOG("Creating Data Tables ... ");

  const int totalDataTables =
      static_cast<int>(opts.dataTableName.size());

  for (int tableIndex = 0;
       tableIndex < totalDataTables;
       ++tableIndex) {
    const int errorCode = modelManager->CreateDataTable(
        const_cast<char*>(opts.dataTableName[tableIndex].c_str()),
        const_cast<char*>(opts.dataTableType[tableIndex].c_str()),
        opts.dataTableVals[tableIndex]);

    if (errorCode == CV_ERROR) {
      const std::string errMsg =
          "ERROR: Error Creating DATATABLE " +
          std::to_string(tableIndex) + "\n";
      throw cvException(errMsg.c_str());
    }
  }
}

void createModelSegments(
    const cvOneD::options& opts,
    cvOneDModelManager* modelManager) {
  SVONED_DEBUG_LOG("Creating Segments ... ");

  const int totalSegments =
      static_cast<int>(opts.segmentName.size());

  for (int segmentIndex = 0;
       segmentIndex < totalSegments;
       ++segmentIndex) {
    const auto& materialName =
        opts.segmentMatName[segmentIndex];
    const int materialID =
        getListIDWithStringKey(materialName, opts.materialName);

    if (materialID < 0) {
    const auto errMsg =
        "ERROR: Segment '" + opts.segmentName[segmentIndex] +
        "' references material '" + materialName +
        "', but that material is not defined in the 1D input file.";
    throw cvException(errMsg.c_str());
    }

    std::vector<double> curveTimes;
    std::vector<double> curveValues;

    const auto& curveName =
        opts.segmentDataTableName[segmentIndex];

    if (upper_string(curveName) != "NONE") {
      const int dataTableID =
          getDataTableID(curveName);
      const int curveSize =
          cvOneDGlobal::gDataTables[dataTableID]->getSize();

      curveTimes.reserve(curveSize);
      curveValues.reserve(curveSize);

      for (int curveIndex = 0;
           curveIndex < curveSize;
           ++curveIndex) {
        curveTimes.push_back(
            cvOneDGlobal::gDataTables[dataTableID]->getTime(curveIndex));
        curveValues.push_back(
            cvOneDGlobal::gDataTables[dataTableID]->getValues(curveIndex));
      }
    } else {
      curveTimes.push_back(0.0);
      curveValues.push_back(0.0);
    }

    modelManager->CreateSegment(
        const_cast<char*>(opts.segmentName[segmentIndex].c_str()),
        static_cast<long>(opts.segmentID[segmentIndex]),
        opts.segmentLength[segmentIndex],
        static_cast<long>(opts.segmentTotEls[segmentIndex]),
        static_cast<long>(opts.segmentInNode[segmentIndex]),
        static_cast<long>(opts.segmentOutNode[segmentIndex]),
        opts.segmentInInletArea[segmentIndex],
        opts.segmentInOutletArea[segmentIndex],
        opts.segmentInFlow[segmentIndex],
        materialID,
        const_cast<char*>(opts.segmentLossType[segmentIndex].c_str()),
        opts.segmentBranchAngle[segmentIndex],
        opts.segmentUpstreamSegment[segmentIndex],
        opts.segmentBranchSegment[segmentIndex],
        const_cast<char*>(opts.segmentBoundType[segmentIndex].c_str()),
        curveValues.data(),
        curveTimes.data(),
        static_cast<int>(curveTimes.size()));
  }
}

struct InletCurveData {
  std::vector<double> times;
  std::vector<double> values;
};

InletCurveData createInletCurveData(
    const cvOneD::options& opts,
    const std::string& couplingType) {
  InletCurveData curveData;
  const std::string& curveName = opts.inletDataTableName;

  if (upper_string(curveName) == "NONE") {
    curveData.times.push_back(0.0);
    curveData.values.push_back(0.0);
    return curveData;
  }

  const int dataTableID =
      getDataTableID(curveName);
  const int curveSize =
      cvOneDGlobal::gDataTables[dataTableID]->getSize();

  curveData.times.resize(curveSize);
  curveData.values.resize(curveSize);

  if (couplingType == "DIR") {
    SVONED_DEBUG_LOG(
        "Dirichlet coupling. Get inlet flow data");

    for (int curveIndex = 0;
         curveIndex < curveSize;
         ++curveIndex) {
      curveData.times[curveIndex] =
          cvOneDGlobal::gDataTables[dataTableID]->getTime(curveIndex);
      curveData.values[curveIndex] =
          cvOneDGlobal::gDataTables[dataTableID]->getValues(curveIndex);
    }
  }

  return curveData;
}

BoundCondTypeScope::BoundCondType getBoundaryConditionType(
    const std::string& boundaryType) {
  const std::string normalizedType = upper_string(boundaryType);

  SVONED_DEBUG_LOG(
      "Inlet Condition Type: " << normalizedType);

  if (normalizedType == "NOBOUND") {
    return BoundCondTypeScope::NOBOUND;
  }
  if (normalizedType == "PRESSURE") {
    return BoundCondTypeScope::PRESSURE;
  }
  if (normalizedType == "PRESSURE_WAVE") {
    return BoundCondTypeScope::PRESSURE_WAVE;
  }
  if (normalizedType == "FLOW") {
    return BoundCondTypeScope::FLOW;
  }
  if (normalizedType == "RESISTANCE") {
    return BoundCondTypeScope::RESISTANCE;
  }
  if (normalizedType == "RESISTANCE_TIME") {
    return BoundCondTypeScope::RESISTANCE_TIME;
  }
  if (normalizedType == "RCR") {
    return BoundCondTypeScope::RCR;
  }
  if (normalizedType == "CORONARY") {
    return BoundCondTypeScope::CORONARY;
  }
  if (normalizedType == "COUPLED") {
    return BoundCondTypeScope::COUPLED;
  }

  const std::string errMsg =
      "ERROR: Invalid inlet boundary condition type: " +
      boundaryType;
  throw cvException(errMsg.c_str());
}

void initializeSolver(
    const cvOneD::options& opts,
    char* couplingTypes,
    int& systemSize) {
  const std::string couplingType(couplingTypes);
  InletCurveData inletCurveData =
      createInletCurveData(opts, couplingType);

  cvOneDGlobal::isCreating = false;
  const auto boundaryType =
      getBoundaryConditionType(opts.boundaryType);

  cvOneDMthSegmentModel::STABILIZATION = opts.useStab;
  cvOneDGlobal::CONSERVATION_FORM = opts.useIV;
  cvOneDBFSolver::ASCII = 1;

  cvOneDBFSolver::SetModelPtr(
      cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]);

  cvOneDBFSolver::SetDeltaTime(opts.timeStep);
  cvOneDBFSolver::SetStepSize(opts.stepSize);
  cvOneDBFSolver::SetMaxStep(opts.maxStep);
  cvOneDBFSolver::SetQuadPoints(opts.quadPoints);
  cvOneDBFSolver::SetInletBCType(boundaryType);

  if (couplingType == "DIR") {
    cvOneDBFSolver::DefineInletFlow(
        inletCurveData.times.data(),
        inletCurveData.values.data(),
        static_cast<int>(inletCurveData.times.size()));
  }

  cvOneDBFSolver::SetConvergenceCriteria(
      opts.convergenceTolerance);

  cvOneDGlobal::isSolving = true;

  cvOneDBFSolver::Solve_initi(systemSize, couplingTypes);

  SVONED_DEBUG_LOG(
      "[initialize_1d] Model initialized, "
      "preparing for time-stepping...");

  cvOneDBFSolver::InitializeAllEquations();

  SVONED_DEBUG_LOG(
      "[initialize_1d] Equations initialized successfully");
}

}  // namespace

OneDSolverInterface* OneDSolverInterface::get_interface(
    int problem_id) {
  const auto iter = interface_list_.find(problem_id);

  if (iter == interface_list_.end()) {
    const std::string errMsg =
        "ERROR: 1D interface problem_id " +
        std::to_string(problem_id) +
        " was not found.";
    throw cvException(errMsg.c_str());
  }

  return iter->second;
}

//////////////////////////////////////////////////////////
//          Callable C interface functions              //
//////////////////////////////////////////////////////////

extern "C" void initialize_1d(const char* input_file, int& problem_id,  int& systemSize,
                              char* coupling_types);

extern "C" void set_external_step_size_1d(int problem_id,
                                          double external_step_size);

extern "C" void return_1d_solution(int problem_id, double* solution_1d, int solution_size);

extern "C" void update_1d_solution(int problem_id, const double* previous_solution_data, int solution_size);

extern "C" void run_1d_simulation_step_1d(int problem_id, double current_time, int save_time, char* coupling_types, double* params,
                                          double* solution_vector, double& cplBCvalue, char* last_flag, int& error_code);
extern "C" void extract_coupled_dof(int problem_id, int& coupled_dof, char* coupling_types);

/**
 * @brief Initialize the 1D solver interface for 3D-1D coupling.
 *
 * @param input_file The name of the 1D solver input file (.in format)
 * @param problem_id The returned ID used to identify the 1D problem
 * @param coupling_types Type of coupling for each surface ("DIR" or "NEU")
 */
void initialize_1d(const char* input_file, int& problem_id, int& systemSize,
                   char* coupling_types) {
  SVONED_DEBUG_LOG("========== svOneD initialize ==========");
  try {
    std::string input_file_str(input_file);
    SVONED_DEBUG_LOG("[initialize] input_file: " << input_file_str);
    // Create interface object
    auto interface = new OneDSolverInterface(input_file_str);
    problem_id = interface->problem_id_;
    SVONED_DEBUG_LOG("[initialize] problem_id: " << (int)problem_id);
    // Parse input file and read 1D configuration
    char* argv[] = { (char*)"svOneDSolver", const_cast<char*>(input_file) };
    auto const simulationOptions = parseArgsAndHandleOptions(2, argv);

    if (simulationOptions) {
        SVONED_DEBUG_LOG("[initialize] Simulation options loaded successfully");
        SVONED_DEBUG_LOG("[initialize] Model name: " << simulationOptions->modelName);
        SVONED_DEBUG_LOG("[initialize] Time step: " << simulationOptions->timeStep);
        SVONED_DEBUG_LOG("[initialize] Step size: " << simulationOptions->stepSize);
        SVONED_DEBUG_LOG("[initialize] Max steps: " << simulationOptions->maxStep);
        // Store options in interface for later use
        interface->model_name_ = simulationOptions->modelName;
        interface->time_step_size_ = simulationOptions->stepSize;
        interface->max_step_ = simulationOptions->maxStep;

        // Store coupling options in interface
        interface->coupling_status_ = simulationOptions->couplingStatus;
        interface->coupling_type_ = simulationOptions->couplingType;
        interface->coupling_substeps_ = simulationOptions->couplingSubsteps;

        SVONED_DEBUG_LOG("[initialize] Coupling status: " << interface->coupling_status_);
        SVONED_DEBUG_LOG("[initialize] Coupling type: " << interface->coupling_type_);
        SVONED_DEBUG_LOG("[initialize] Coupling substeps: " << static_cast<int>(interface->coupling_substeps_));
        //// use functions inside runOneDSolver
        // Model checking
        const cvOneD::options& opts = *simulationOptions;
        cvOneD::validateOptions(opts);
        SVONED_DEBUG_LOG("[initialize] validationOptions completed");
        // Not sure what this function do right now
        setOutputGlobals(opts);
        SVONED_DEBUG_LOG("[initialize] setOutputGlobals completed");
        // Use functions inside createAndRunModel
        // only use functions for creating model
        SVONED_DEBUG_LOG("[initialize] creating model: " << opts.modelName);
        // Create model manager
        cvOneDModelManager* oned = new cvOneDModelManager((char*)opts.modelName.c_str());

        // Create nodes
        createModelNodes(opts, oned);
        // Create joints
        createModelJoints(opts, oned);
        // Create materials
        createModelMaterials(opts, oned);
        // Create datatables
        createModelDataTables(opts, oned);
        // Create segments
        createModelSegments(opts, oned);
        // Initialize solver
        initializeSolver(opts, coupling_types, systemSize);

        // Time loop initialization
        interface->time_step_ = 0;

        SVONED_DEBUG_LOG("[initialize_1d] 1D model is ready for time-stepping");
    }else {
        SVONED_DEBUG_LOG("[initialize_1d] WARNING: No simulation options found");
    }
    SVONED_DEBUG_LOG("[initialize_1d] 1D model initialized successfully");
    } catch (const std::exception& e) {
        problem_id = -1;

        const std::string errMsg =
            "ERROR: Failed to initialize the 1D solver: " +
            std::string(e.what());

        throw cvException(errMsg.c_str());
    }
}

/**
 * @brief Set the external (3D solver) time step size for the 1D interface.
 *
 * @param problem_id The ID used to identify the 1D problem.
 * @param external_step_size The time step size of the external program.
 */
void set_external_step_size_1d(int problem_id, double external_step_size) {
  auto interface = OneDSolverInterface::get_interface(problem_id);

  double oned_step_size = external_step_size/ (double(interface->coupling_substeps_));
  // interface->coupling_substeps_ : number of sub steps that 1D solver will take within one external step.
  // For example, if coupling_substeps_ = 2, then 1D solver will take 2 steps within one external step
  // so each 1D step size will be external_step_size/2.

  interface->external_step_size_ = oned_step_size;//now 3D and 1D has same step size

  cvOneDBFSolver::SetDeltaTime(oned_step_size);// set deltaTime in solver as oned_step_size


  SVONED_DEBUG_LOG("[set_external_step_size_1d] Step size: " << oned_step_size << " s");
}

/**
 * @brief copy current 1D solution from 1D solver to 3D solver
 *
 * @param problem_id The ID used to identify the 1D problem.
 * @param solution_1d 1D solution of 3D solver
 */
void return_1d_solution(int problem_id, double* solution_1d, int solution_size){
    auto interface = OneDSolverInterface::get_interface(problem_id);

    try {
        SVONED_DEBUG_LOG(
            "[return_1d_solution] Extracting solution for problem_id: "
            << static_cast<int>(problem_id));

        cvOneDBFSolver::GetCurrentSolution(solution_1d, solution_size);

        SVONED_DEBUG_LOG("[return_1d_solution] Solution extracted successfully");
    } catch (const std::exception& e) {
    const std::string errMsg =
        "ERROR: Failed to return the 1D solution: " +
        std::string(e.what());
    throw cvException(errMsg.c_str());
    }
}

/**
 * @brief Reset the 1D solver solution vectors for coupling with 3D solver.
 *
 * This function initializes previousSolution and currentSolution with values from
 * the previous time step, which is required when the 1D solver is called multiple times
 * within a 3D Newton loop.
 *
 * @param problem_id The ID of the 1D problem
 * @param previous_solution_data Array containing the solution from previous time step
 * @param solution_size The size of the solution array
 */
void update_1d_solution(int problem_id, const double* previous_solution_data, int solution_size) {
    auto interface = OneDSolverInterface::get_interface(problem_id);

    try {
        SVONED_DEBUG_LOG("[update_1d_solution] ========================================");
        // call function from cvOneDBFSolver
        cvOneDBFSolver::InitializeSolutionFromVector(previous_solution_data, solution_size);

        SVONED_DEBUG_LOG("[update_1d_solution] Solution vectors reset successfully");
    } catch (const std::exception& e) {
        const std::string errMsg =
            "ERROR: Failed to update the 1D solution: " +
            std::string(e.what());
        throw cvException(errMsg.c_str());
    }
}


/**
 * @brief Run one time step of the 1D simulation.
 *
 * Performs a single time step of 1D blood flow simulation by calling the
 * underlying solver's SolveSingleTimeStep function. Returns the solution
 * vector (pressure and flow) after the time step.
 *
 * @param problem_id The ID used to identify the 1D problem.
 * @param current_time Current simulation time (seconds)
 * @param save_time Save 1D results every save_time 3D time steps
 * @param coupling_types Type of coupling for each surface ("DIR" or "NEU")
 * @param params bc data from 3D solver for interpolation [2, t1, t2, value1, value2]
 * @param solution_vector Output vector containing pressure and flow values
 *                       Size: 2*system_size (pressure[0..system_size-1], flow[system_size..2*system_size-1])
 * @param error_code Output error code (0 = success, <0 = error)
 */
void run_1d_simulation_step_1d(int problem_id, double current_time, int save_time, char* coupling_types, double* params,
                                double* solution_vector, double& cplBCvalue, char* last_flag,int& error_code) {
  auto interface = OneDSolverInterface::get_interface(problem_id);

  error_code = 0;
  auto oned_substeps = interface->coupling_substeps_;

  try {
    SVONED_DEBUG_LOG(
    "[run_1d_simulation_step_1d] Time step "
    << static_cast<int>(interface->time_step_)
    << ", current_time = " << current_time << " s");


    // Update current time in solver
    cvOneDBFSolver::currentTime = current_time;

    // now params[0] is always 2. not used here
    cvOneDFEAVector* final_solution = nullptr;
    double t1 = params[1];
    double t2 = params[2];
    double val1 = params[3];
    double val2 = params[4];
    double alpha = 0.0;
    double interpolated_value = 0.0;

    for (int i = 0; i < oned_substeps; i++) {
        // Evaluate the 3D-coupled BC value at the end of this 1D substep.
        // If dt_1D == dt_3D, eval_time becomes t_new and params[4] is used.
        const double eval_time = current_time + interface->external_step_size_;

        if (eval_time <= t1) {
            interpolated_value = val1;
        } else if (eval_time >= t2) {
            interpolated_value = val2;
        } else {
            alpha = (eval_time - t1) / (t2 - t1);
            interpolated_value = val1 + alpha * (val2 - val1);
        }

        final_solution = cvOneDBFSolver::SolveSingleTimeStep(eval_time, interpolated_value);

        current_time = eval_time;
    }

    // Check if the solver returned a valid solution
    if (final_solution == nullptr) {
    throw cvException(
        "ERROR: SolveSingleTimeStep returned a null solution.");
    }

    // Step 1: Copy final_solution directly to solution_vector [area1][flow1][area2][flow2]...
    int solution_size = final_solution->GetDimension();
    double* solution_data_tmp = final_solution->GetEntries();
    for(int i = 0; i < solution_size; i++) {
        solution_vector[i] = solution_data_tmp[i];
    }

    // Step 2: Convert to flow and pressure format separately (for internal use or debugging)
    // Create a separate vector for converted solution
    // Make solution vector to transfer to 3D solver
    // converted_solution format: [flow1][pressure1][flow2][pressure2]... for each nodes
    // current final_solution: [area1][flow1][area2][flow2]... for each nodes
    std::vector<double> converted_solution(solution_size);

    cvOneDBFSolver::ConvertSolutionToFlowPressure(
        final_solution, converted_solution.data());

    cvOneDBFSolver::extractCplBC(
        converted_solution.data(),
        cplBCvalue,
        coupling_types);

    // last_flag is "L" when it is at the last iteration
    // last_flag is "D: when it is at the middel of 3D newton iteration
    if (std::string(last_flag) == "L"){
        // Update interface internal states
        interface->time_step_++; // this is for 1D internal use only. substep is not inlcuded.
        // this is same time_step with 3D solver

        // Write the accepted 1D solution at the requested output interval.
        // This is executed only after the final 3D Newton iteration.
        if (interface->time_step_ % save_time == 0) {
            SVONED_DEBUG_LOG("generate vtk file at time step: "<< static_cast<int>(interface->time_step_));
            cvOneDBFSolver::postprocess_VTK_XML3D_SingleTimeStep(interface->time_step_, final_solution);
            SVONED_DEBUG_LOG("generate txt files at time step: "<< static_cast<int>(interface->time_step_));
            cvOneDBFSolver::postprocess_Text_SingleTimeStep(interface->time_step_, final_solution);
        }

    }

    SVONED_DEBUG_LOG("[run_1d_simulation_step_1d] Time step completed");
    } catch (const std::exception& e) {
    error_code = -1;

    const std::string errMsg =
        "ERROR: Failed to run a 1D simulation step: " +
        std::string(e.what());
    throw cvException(errMsg.c_str());
    }
}

void extract_coupled_dof(int problem_id, int& coupled_dof, char* coupling_types){
    auto interface = OneDSolverInterface::get_interface(problem_id);

    try {
        SVONED_DEBUG_LOG("[extract_coupled_dof] ========================================");
        // call function from cvOneDBFSolver
        cvOneDBFSolver::extractCplDOF(coupled_dof, coupling_types);

        SVONED_DEBUG_LOG("[extract_coupled_dof] Coupled DOF extracted successfully");
    } catch (const std::exception& e) {
    const std::string errMsg =
        "ERROR: Failed to extract the coupled 1D DOF: " +
        std::string(e.what());
    throw cvException(errMsg.c_str());
    }

}