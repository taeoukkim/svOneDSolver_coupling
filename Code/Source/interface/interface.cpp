// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the
// University of California, and others. SPDX-License-Identifier: BSD-3-Clause

#include "interface.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <optional>
#include <algorithm>

#include "cvOneDGlobal.h"
#include "cvOneDModelManager.h"
#include "cvOneDOptions.h"
#include "cvOneDOptionsJsonParser.h"
#include "cvOneDOptionsJsonSerializer.h"
#include "cvOneDOptionsLegacySerializer.h"

// ===== Resolve operator<< ambiguity =====
using std::cout;
using std::cerr;
using std::endl;
using std::string;

//////////////////////////////////////////////////////////
//      Helper functions from main.cxx                  //
//////////////////////////////////////////////////////////

std::string removeQuotesIfPresent(const std::string& str) {
    std::string result = str;
    if (!result.empty() && result.front() == '"' && result.back() == '"') {
        result = result.substr(1, result.length() - 2); 
    }
    return result;
}

struct ArgOptions {
  std::optional<std::string> jsonInput = std::nullopt;
  std::optional<std::string> legacyConversionInput = std::nullopt;
  std::optional<std::string> jsonConversionOutput = std::nullopt;
};


cvOneD::options readLegacyOptions(std::string const& inputFile) {
  cvOneD::options opts{};
  cvOneD::readOptionsLegacyFormat(inputFile, &opts);
  return opts;
}

std::optional<cvOneD::options> parseArgsAndHandleOptions(int argc, char** argv) {
  // Legacy behavior: single input file
  if(argc == 2) {
    string inputFile{removeQuotesIfPresent(argv[1])};
    return readLegacyOptions(inputFile);
  }
  
  return std::nullopt;
}

void setOutputGlobals(const cvOneD::options& opts){  

  if(upper_string(opts.outputType) == "TEXT"){
    cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_TEXT;
  }else if(upper_string(opts.outputType) == "VTK"){
    cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_VTK;
  }else if(upper_string(opts.outputType) == "BOTH"){
    cvOneDGlobal::outputType = OutputTypeScope::OUTPUT_BOTH;
  }else{
    throw cvException("ERROR: Invalid OUTPUT Type.\n");
  }

  if(opts.vtkOutputType){
    cvOneDGlobal::vtkOutputType = *opts.vtkOutputType;
    if(cvOneDGlobal::vtkOutputType > 1){
      throw cvException("ERROR: Invalid OUTPUT VTK Type.\n");
    }
  }
  
}

size_t findJointNodeIndexOrThrow(const auto& jointNodeName, const auto& nodeNames, const auto& jointName){
  // Return the index of "nodeNames" that corresponds to the "jointNodeName"
  // or throw a context-specific error.

  auto const iter = std::find(nodeNames.begin(), nodeNames.end(), jointNodeName);
  if (iter == nodeNames.end()) {
    std::string const errMsg = "ERROR: The node '" + jointNodeName + "' required by joint '" 
      + jointName + "' was not found in the list of nodes.";
    throw cvException(errMsg.c_str());
  }
  return std::distance(nodeNames.begin(), iter); 
}

int getDataTableIDFromStringKey(string key){
  bool found = false;
  int count = 0;
  while((!found)&&(count<cvOneDGlobal::gDataTables.size())){
    found = upper_string(key) == upper_string(cvOneDGlobal::gDataTables[count]->getName());
    // Update Counter
    if(!found){
      count++;
    }
  }
  if(!found){
    throw cvException(string("ERROR: Cannot find data table entry from string key: " + key + ".\n").c_str());
    return -1;
  }else{
    return count;
  }
}

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

//////////////////////////////////////////////////////////
//          Callable C interface functions              //
//////////////////////////////////////////////////////////

extern "C" void initialize_1d(const char* input_file, int& problem_id,  int& systemSize,
                              char* coupling_types);

extern "C" void set_external_step_size_1d(int problem_id,
                                          double external_step_size);

extern "C" void update_coupled_bc_1d(int problem_id, int num_surfaces,
                                     double* input_values,
                                     const char* coupling_type);

extern "C" void get_coupled_solution_1d(int problem_id, int num_surfaces,
                                        double* flows_out,
                                        double* pressures_out);

extern "C" void run_1d_simulation_step_1d(int problem_id, double current_time,
                                          double* solution_vector, int& error_code);

extern "C" void get_resistance_matrix_1d(int problem_id, int num_surfaces,
                                         double** resistance_matrix);

extern "C" void cleanup_1d(int problem_id);

/**
 * @brief Initialize the 1D solver interface for 3D-1D coupling.
 *
 * @param input_file The name of the 1D solver input file (.in format)
 * @param problem_id The returned ID used to identify the 1D problem
 * @param coupling_types Type of coupling for each surface ("DIR" or "NEU")
 */
void initialize_1d(const char* input_file, int& problem_id, int& systemSize,
                   char* coupling_types) {
  cout << "========== svOneD initialize ==========" << endl;
  
  try {
    std::string input_file_str(input_file);
    cout << "[initialize] input_file: " << input_file_str << endl;

    // Create interface object
    auto interface = new OneDSolverInterface(input_file_str);
    problem_id = interface->problem_id_;
    cout << "[initialize] problem_id: " << (int)problem_id << endl;

    // Parse input file and read 1D configuration
    char* argv[] = { (char*)"svOneDSolver", const_cast<char*>(input_file) };
    auto const simulationOptions = parseArgsAndHandleOptions(2, argv);

    if (simulationOptions) {
        cout << "[initialize] Simulation options loaded successfully" << endl;
        cout << "[initialize] Model name: " << simulationOptions->modelName << endl;
        cout << "[initialize] Time step: " << simulationOptions->timeStep << endl;
        cout << "[initialize] Step size: " << simulationOptions->stepSize << endl;
        cout << "[initialize] Max steps: " << simulationOptions->maxStep << endl;
        
        // Store options in interface for later use
        interface->model_name_ = simulationOptions->modelName;
        interface->time_step_size_ = simulationOptions->stepSize;
        interface->max_step_ = simulationOptions->maxStep;

        //// use functions inside runOneDSolver
        // Model checking
        const cvOneD::options& opts = *simulationOptions;
        cvOneD::validateOptions(opts);
        cout << "[initialize] validationOptions completed" << endl;

        // Not sure what this function do right now
        setOutputGlobals(opts);
        cout << "[initialize] setOutputGlobals completed" << endl;

        // Use functions inside createAndRunModel
        // only use functions for creating model
        cout << "[initialize] creating model: " << opts.modelName << endl;
        // Create model manager
        cvOneDModelManager* oned = new cvOneDModelManager((char*)opts.modelName.c_str());

        // Create nodes
        printf("Creating Nodes ... \n");
        int totNodes = opts.nodeName.size();
        int nodeError = CV_OK;
        for(int loopA = 0; loopA < totNodes; loopA++) {
            // Finally Create Joint
            nodeError = oned->CreateNode((char*)opts.nodeName[loopA].c_str(),
                                        opts.nodeXcoord[loopA], opts.nodeYcoord[loopA], opts.nodeZcoord[loopA]);
            if(nodeError == CV_ERROR) {
            throw cvException(string("ERROR: Error Creating NODE " + to_string(loopA) + "\n").c_str());
            }
        }

        // Create joints
        printf("Creating Joints ... \n");
        int totJoints = opts.jointName.size();
        int jointError = CV_OK;
        int* asInlets = nullptr;
        int* asOutlets = nullptr;
        string currInletName;
        string currOutletName;
        int jointInletID = 0;
        int jointOutletID = 0;
        int totJointInlets = 0;
        int totJointOutlets = 0;
        for(int loopA = 0; loopA < totJoints; loopA++) {
            // GET NAMES FOR INLET AND OUTLET
            currInletName = opts.jointInletName[loopA];
            currOutletName = opts.jointOutletName[loopA];
            // FIND JOINTINLET INDEX
            jointInletID = getListIDWithStringKey(currInletName, opts.jointInletListNames);
            if(jointInletID < 0) {
            throw cvException(string("ERROR: Cannot Find JOINTINLET for key " + currInletName).c_str());
            }
            totJointInlets = opts.jointInletListNumber[jointInletID];
            // FIND JOINTOUTLET INDEX
            jointOutletID = getListIDWithStringKey(currOutletName, opts.jointOutletListNames);
            if(jointInletID < 0) {
            throw cvException(string("ERROR: Cannot Find JOINTOUTLET for key " + currOutletName).c_str());
            }
            // GET TOTALS
            totJointInlets = opts.jointInletListNumber[jointInletID];
            totJointOutlets = opts.jointOutletListNumber[jointOutletID];
            // ALLOCATE INLETS AND OUTLET LIST
            asInlets = nullptr;
            asOutlets = nullptr;
            if(totJointInlets > 0) {
            asInlets = new int[totJointInlets];
            for(int loopB = 0; loopB < totJointInlets; loopB++) {
                asInlets[loopB] = opts.jointInletList[jointInletID][loopB];
            }
            }
            if(totJointOutlets > 0) {
            asOutlets = new int[totJointOutlets];
            for(int loopB = 0; loopB < totJointOutlets; loopB++) {
                asOutlets[loopB] = opts.jointOutletList[jointOutletID][loopB];
            }
            }

            // Find the index of the indicated node.
            auto const jointName = opts.jointName.at(loopA);
            auto const nodeIndex = findJointNodeIndexOrThrow( 
            opts.jointNode.at(loopA), opts.nodeName, jointName);

            // Finally Create Joint
            jointError = oned->CreateJoint(jointName.c_str(),
                                        opts.nodeXcoord[nodeIndex], opts.nodeYcoord[nodeIndex], opts.nodeZcoord[nodeIndex],
                                        totJointInlets, totJointOutlets, asInlets, asOutlets);
            if(jointError == CV_ERROR) {
            throw cvException(string("ERROR: Error Creating JOINT " + to_string(loopA) + "\n").c_str());
            }
            // Deallocate
            delete[] asInlets;
            delete[] asOutlets;
            asInlets = nullptr;
            asOutlets = nullptr;
        }

        // Create materials
        printf("Creating Materials ... \n");
        int totMaterials = opts.materialName.size();
        int matError = CV_OK;
        double doubleParams[3];
        int matID = 0;
        string currMatType = "MATERIAL_OLUFSEN";
        int numParams = 0;
        for(int loopA = 0; loopA < totMaterials; loopA++) {
            if(upper_string(opts.materialType[loopA]) == "OLUFSEN") {
            currMatType = "MATERIAL_OLUFSEN";
            numParams = 3;
            } else {
            currMatType = "MATERIAL_LINEAR";
            numParams = 1;
            }
            doubleParams[0] = opts.materialParam1[loopA];
            doubleParams[1] = opts.materialParam2[loopA];
            doubleParams[2] = opts.materialParam3[loopA];
            // CREATE MATERIAL
            matError = oned->CreateMaterial((char*)opts.materialName[loopA].c_str(),
                                            (char*)currMatType.c_str(),
                                            opts.materialDensity[loopA],
                                            opts.materialViscosity[loopA],
                                            opts.materialExponent[loopA],
                                            opts.materialPRef[loopA],
                                            numParams, doubleParams,
                                            &matID);
            if(matError == CV_ERROR) {
            throw cvException(string("ERROR: Error Creating MATERIAL " + to_string(loopA) + "\n").c_str());
            }
        }

        // Create datatables
        printf("Creating Data Tables ... \n");
        int totCurves = opts.dataTableName.size();
        int curveError = CV_OK;
        for(int loopA = 0; loopA < totCurves; loopA++) {
            curveError = oned->CreateDataTable((char*)opts.dataTableName[loopA].c_str(),(char*)opts.dataTableType[loopA].c_str(), opts.dataTableVals[loopA]);
            if(curveError == CV_ERROR) {
            throw cvException(string("ERROR: Error Creating DATATABLE " + to_string(loopA) + "\n").c_str());
            }
        }

        // Create segments
        printf("Creating Segments ... \n");
        int segmentError = CV_OK;
        int totalSegments = opts.segmentName.size();
        int curveTotals = 0;
        double* curveTime = nullptr;
        double* curveValue = nullptr;
        string matName;
        string curveName;
        int currMatID = 0;
        int dtIDX = 0;
        for(int loopA = 0; loopA < totalSegments; loopA++) {

            // GET MATERIAL
            matName = opts.segmentMatName[loopA];
            currMatID = getListIDWithStringKey(matName, opts.materialName);
            if(currMatID < 0) {
            throw cvException(string("ERROR: Cannot Find Material for key " + matName).c_str());
            }

            // GET CURVE DATA
            curveName = opts.segmentDataTableName[loopA];

            if(upper_string(curveName) != "NONE") {
            dtIDX = getDataTableIDFromStringKey(curveName);
            curveTotals = cvOneDGlobal::gDataTables[dtIDX]->getSize();
            curveTime = new double[curveTotals];
            curveValue = new double[curveTotals];
            for(int loopA = 0; loopA < curveTotals; loopA++) {
                curveTime[loopA] = cvOneDGlobal::gDataTables[dtIDX]->getTime(loopA);
                curveValue[loopA] = cvOneDGlobal::gDataTables[dtIDX]->getValues(loopA);
            }
            } else {
            curveTotals = 1;
            curveTime = new double[curveTotals];
            curveValue = new double[curveTotals];
            curveTime[0] = 0.0;
            curveValue[0] = 0.0;
            }
            segmentError = oned->CreateSegment((char*)opts.segmentName[loopA].c_str(),
                                            (long)opts.segmentID[loopA],
                                            opts.segmentLength[loopA],
                                            (long)opts.segmentTotEls[loopA],
                                            (long)opts.segmentInNode[loopA],
                                            (long)opts.segmentOutNode[loopA],
                                            opts.segmentInInletArea[loopA],
                                            opts.segmentInOutletArea[loopA],
                                            opts.segmentInFlow[loopA],
                                            currMatID,
                                            (char*)opts.segmentLossType[loopA].c_str(),
                                            opts.segmentBranchAngle[loopA],
                                            opts.segmentUpstreamSegment[loopA],
                                            opts.segmentBranchSegment[loopA],
                                            (char*)opts.segmentBoundType[loopA].c_str(),
                                            curveValue,
                                            curveTime,
                                            curveTotals);
            if(segmentError == CV_ERROR) {
            throw cvException(string("ERROR: Error Creating SEGMENT " + to_string(loopA) + "\n").c_str());
            }
            // Deallocate
            delete[] curveTime;
            curveTime = nullptr;
            delete[] curveValue;
            curveValue = nullptr;
        }// until this part of the code is the part of the code of createAndUnModel before SOLVE MODEL

        ////////////////////////////////////////////////////////////////////////////////////////////
        // 이건 1D 모델에서 inlet boundary condition을 설정하는 부분
        // 커플링이 Dirichlet인 경우 여전히 필요
        // 커플링이 Neumann인 경우에는 flow 값을 각 시간에서 풀때 3D 솔버에서 받아와서 업데이트 해주는 방식으로 구현할 예정
        // 일딴 지금 계발에선 항상 커플링이 Dirichlet인 경우로 가정
        string inletCurveName = opts.inletDataTableName;
        int inletCurveIDX = getDataTableIDFromStringKey(inletCurveName);
        int inletCurveTotals = cvOneDGlobal::gDataTables[inletCurveIDX]->getSize();
        double* inletCurveTime = new double[inletCurveTotals];
        double* inletCurveValue = new double[inletCurveTotals];
        if(std::string(coupling_types) == "DIR"){
            // Dirichlet coupling인 경우 inlet boundary condition으로 사용할 데이터를 미리 저장해 놓는 과정
            printf("Dirichlet coupling. Get inlet flow data ... \n");
            for(int loopB = 0; loopB < inletCurveTotals; loopB++) {
                inletCurveTime[loopB] = cvOneDGlobal::gDataTables[inletCurveIDX]->getTime(loopB);
                inletCurveValue[loopB] = cvOneDGlobal::gDataTables[inletCurveIDX]->getValues(loopB);
            }
        }
        // coupling code에서 반드시 필요한 부분임. if문만 나중에 다듬기
        ////////////////////////////////////////////////////////////////////////////////////////////

        // 여기부터는 SolveModel중 초기화에 해당되는 부분을 추가
        int solveError = CV_OK;
        // inlet boundary type이 pts.boundaryType에 저장되어있는데 (.in 파일에서 SOLVEROPTIONS에서 읽은것)
        // 아래 코드는 SolveModel에서 그대로 가지고 온건데 왜 여러 다른 boundary type이 있는지 아직 잘 모르겠음.
        // TODO: 나중에 Neumann coupling의 경우 이 inlet BC이 coupled 같은게 될것이므로 수정 필요
        // 일딴은 그냥 넣어보자
        BoundCondTypeScope::BoundCondType boundT;

        // set the creation flag to off.
        cvOneDGlobal::isCreating = false;
        char* boundType_tmp = (char*)opts.boundaryType.c_str();
        // convert char string to boundary condition type
        if(!strcmp( boundType_tmp, "NOBOUND")){
            boundT = BoundCondTypeScope::NOBOUND;
            printf("Inlet Condition Type: NOBOUND\n");
        }else if(!strcmp( boundType_tmp, "PRESSURE")){
            boundT = BoundCondTypeScope::PRESSURE;
            printf("Inlet Condition Type: PRESSURE\n");
        }else if(!strcmp( boundType_tmp, "PRESSURE_WAVE")){
            boundT = BoundCondTypeScope::PRESSURE_WAVE;
            printf("Inlet Condition Type: PRESSURE_WAVE\n");
        }else if(!strcmp( boundType_tmp, "FLOW")){
            boundT = BoundCondTypeScope::FLOW;
            printf("Inlet Condition Type: FLOW\n");
        }else if(!strcmp( boundType_tmp, "RESISTANCE")){
            boundT = BoundCondTypeScope::RESISTANCE;
            printf("Inlet Condition Type: RESISTANCE\n");
        }else if(!strcmp( boundType_tmp, "RESISTANCE_TIME")){
            boundT = BoundCondTypeScope::RESISTANCE_TIME;
            printf("Inlet Condition Type: RESISTANCE_TIME\n");
        }else if(!strcmp( boundType_tmp, "RCR")){
            boundT = BoundCondTypeScope::RCR;
            printf("Inlet Condition Type: RCR\n");
        }else if(!strcmp( boundType_tmp, "CORONARY")){
            boundT = BoundCondTypeScope::CORONARY;
            printf("Inlet Condition Type: CORONARY\n");
        }else{
            solveError = CV_ERROR;
        }

        // Set Solver Options
        cvOneDMthSegmentModel::STABILIZATION = opts.useStab; // 1=stabilization, 0=none
        cvOneDGlobal::CONSERVATION_FORM = opts.useIV;
        cvOneDBFSolver::ASCII = 1;

        // enroll model inside solver of cvOneDBFSolver file
        cvOneDBFSolver::SetModelPtr(cvOneDGlobal::gModelList[cvOneDGlobal::currentModel]);

        // We need to get these from the solver
        cvOneDBFSolver::SetDeltaTime(opts.timeStep);
        cvOneDBFSolver::SetStepSize(opts.stepSize);
        cvOneDBFSolver::SetMaxStep(opts.maxStep);
        cvOneDBFSolver::SetQuadPoints(opts.quadPoints);
        cvOneDBFSolver::SetInletBCType(boundT);
        if(std::string(coupling_types) == "DIR"){
            cvOneDBFSolver::DefineInletFlow(inletCurveTime, inletCurveValue, inletCurveTotals); // this need to be chaged it only happen for Dir
        }
        cvOneDBFSolver::SetConvergenceCriteria(opts.convergenceTolerance);

        cvOneDGlobal::isSolving = true;
        // 여기까지가 SolveModel에서 Solve() 함수 이전까지의 부분

        // 여기부터는 Solve() 함수 부분 중에서도 GenerateSolution 이전부분
        cvOneDBFSolver::Solve_initi(systemSize); // TODO: 지금 이 안에도 커플링위해서 바꿔야하는 함수들 많음
        // output: systemSize, which is the total number of unknowns in the system. #NODE x 2 (flow&area)
        // cout << "system size: " << static_cast<int>(systemSize) << endl; // total number of unknows in the system. #NODE x 2 (flow&area) 


        // TODO: 이제 generateSolution 에 들어왔음.
        // time loop 시작 전에 필요한 초기화 작업들
        // 여기부터는 Solve() 함수 부분 중에서도 GenerateSolution 이전부분
        cout << "[initialize_1d] Model initialized, preparing for time-stepping..." << endl;

        //EquationInitialize() 호출
        try {
            cvOneDBFSolver::InitializeAllEquations();
            cout << "[initialize_1d] Equations initialized successfully" << endl;
        } catch (const std::exception& e) {
            cerr << "[initialize_1d] Error initializing equations: " << e.what() << endl;
            throw;
        }

        // Time loop initialization
        interface->current_time_ = 0.0;
        interface->time_step_ = 0;

        cout << "[initialize_1d] 1D model is ready for time-stepping" << endl;
      



        cout << "TEST~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;








    }else {
        cout << "[initialize_1d] WARNING: No simulation options found" << endl;
    }
    cout << "[initialize_1d] 1D model initialized successfully" << endl;
    
  } catch (const std::exception& e) {
    cerr << "Error in initialize_1d: " << e.what() << endl;
    problem_id = -1;
  }
}

/**
 * @brief Set the external (3D solver) time step size for the 1D interface.
 *
 * @param problem_id The ID used to identify the 1D problem.
 * @param external_step_size The time step size of the external program.
 */
void set_external_step_size_1d(int problem_id, double external_step_size) {
  auto it = OneDSolverInterface::interface_list_.find(problem_id);
  if (it == OneDSolverInterface::interface_list_.end()) {
    cerr << "Error: problem_id " << (int)problem_id << " not found" << endl;
    return;
  }
  
  auto interface = it->second;
  interface->external_step_size_ = external_step_size;
  cout << "[set_external_step_size_1d] Step size: " << external_step_size << " s" << endl;
}

/**
 * @brief Update boundary conditions for coupled surfaces.
 *
 * @param problem_id The ID used to identify the 1D problem.
 * @param num_surfaces Number of coupled surfaces
 * @param input_values Array of input values (pressures or flows)
 * @param coupling_type Type of coupling ("DIR" for pressure or "NEU" for flow)
 */
void update_coupled_bc_1d(int problem_id, int num_surfaces,
                          double* input_values,
                          const char* coupling_type) {
  auto it = OneDSolverInterface::interface_list_.find(problem_id);
  if (it == OneDSolverInterface::interface_list_.end()) {
    cerr << "Error: problem_id " << (int)problem_id << " not found" << endl;
    return;
  }

  auto interface = it->second;
  std::string coup_type(coupling_type);

  cout << "[update_coupled_bc_1d] Coupling type: " << coup_type << endl;

  for (int i = 0; i < num_surfaces; i++) {
    if (coup_type == "DIR") {
      interface->current_pressures_[i] = input_values[i];
      cout << "  Surface " << (int)i << " (DIR): P = " << input_values[i] << " mmHg" << endl;
    } else if (coup_type == "NEU") {
      interface->current_flows_[i] = input_values[i];
      cout << "  Surface " << (int)i << " (NEU): Q = " << input_values[i] << " mL/s" << endl;
    }
  }
}

/**
 * @brief Get the current solution (flows and pressures) at coupled surfaces.
 *
 * @param problem_id The ID used to identify the 1D problem.
 * @param num_surfaces Number of coupled surfaces
 * @param flows_out Output array for flows at coupled surfaces
 * @param pressures_out Output array for pressures at coupled surfaces
 */
void get_coupled_solution_1d(int problem_id, int num_surfaces,
                             double* flows_out,
                             double* pressures_out) {
  auto it = OneDSolverInterface::interface_list_.find(problem_id);
  if (it == OneDSolverInterface::interface_list_.end()) {
    cerr << "Error: problem_id " << (int)problem_id << " not found" << endl;
    return;
  }

  auto interface = it->second;

  for (int i = 0; i < num_surfaces; i++) {
    flows_out[i] = interface->current_flows_[i];
    pressures_out[i] = interface->current_pressures_[i];
    
    cout << "[get_coupled_solution_1d] Surface " << (int)i 
         << ": Q = " << flows_out[i] 
         << " mL/s, P = " << pressures_out[i] << " mmHg" << endl;
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
 * @param solution_vector Output vector containing pressure and flow values
 *                       Size: 2*system_size (pressure[0..system_size-1], flow[system_size..2*system_size-1])
 * @param error_code Output error code (0 = success, <0 = error)
 */
void run_1d_simulation_step_1d(int problem_id, double current_time, 
                                double* solution_vector, int& error_code) {
  auto it = OneDSolverInterface::interface_list_.find(problem_id); //problem_id에 해당되는 interface 객체를 찾음
  if (it == OneDSolverInterface::interface_list_.end()) {
    cerr << "[run_1d_simulation_step_1d] Error: problem_id " 
         << static_cast<int>(problem_id) << " not found" << endl;
    error_code = -1;
    return;
  }

  auto interface = it->second; // map에서 problem_id에 해당하는 interface 객체를 가져옴(second: value)
  error_code = 0;

  try {
    cout << "[run_1d_simulation_step_1d] ========================================" << endl;
    cout << "[run_1d_simulation_step_1d] Time step " 
         << static_cast<int>(interface->time_step_)
         << ", current_time = " << current_time << " s" << endl;

    // Store previous solution for coupling calculations
    interface->previous_flows_ = interface->current_flows_;
    interface->previous_pressures_ = interface->current_pressures_;

    // Update current time in solver
    cvOneDBFSolver::currentTime = current_time;
    
    // Call the underlying 1D solver to compute one time step
    // SolveSingleTimeStep returns a pointer to the solution vector
    cvOneDFEAVector* solution_ptr = cvOneDBFSolver::SolveSingleTimeStep(current_time);
    
    if (solution_ptr == nullptr) {
      throw std::runtime_error("SolveSingleTimeStep returned null pointer");
    }
    
    // Extract solution from solver
    // The solution is stored in TotalSolution matrix:
    // TotalSolution[i][0] = pressure at node i
    // TotalSolution[i][1] = flow at node i
    cvOneDModel* model = cvOneDBFSolver::GetModelPtr();
    if (!model) {
      throw std::runtime_error("Model pointer is null in cvOneDBFSolver");
    }

    int num_nodes = model->getNumberOfNodes();
    
    // Fill solution vector: [pressure_0, pressure_1, ..., flow_0, flow_1, ...]
    for (int i = 0; i < num_nodes; i++) {
      // Pressure values (first half)
      solution_vector[i] = cvOneDBFSolver::GetSolution(i, 0);
      
      // Flow values (second half)
      solution_vector[num_nodes + i] = cvOneDBFSolver::GetSolution(i, 1);
    }

    // Update interface internal states
    interface->time_step_++;
    interface->current_time_ = current_time;

    // // Update coupled surface flows and pressures from solution
    // // (This depends on your coupled segment configuration)
    // for (size_t i = 0; i < interface->coupled_segment_ids_.size(); i++) {
    //   int seg_id = interface->coupled_segment_ids_[i];
    //   // Extract flow and pressure from coupled segment
    //   // This requires mapping segment ID to node indices
    //   // TODO: Implement segment-to-node mapping based on your model structure
    // }

    cout << "[run_1d_simulation_step_1d] Time step completed" << endl;
    cout << "[run_1d_simulation_step_1d] Returned " << num_nodes * 2 
         << " solution values" << endl;
    cout << "[run_1d_simulation_step_1d] ========================================" << endl;

  } catch (const std::exception& e) {
    cerr << "[run_1d_simulation_step_1d] Error: " << e.what() << endl;
    error_code = -1;
  }
}



/**
 * @brief Get the resistance matrix (sensitivity) dP/dQ for coupling.
 *
 * @param problem_id The ID used to identify the 1D problem.
 * @param num_surfaces Number of coupled surfaces
 * @param resistance_matrix Output resistance matrix (num_surfaces x num_surfaces)
 */
void get_resistance_matrix_1d(int problem_id, int num_surfaces,
                              double** resistance_matrix) {
  auto it = OneDSolverInterface::interface_list_.find(problem_id);
  if (it == OneDSolverInterface::interface_list_.end()) {
    cerr << "Error: problem_id " << (int)problem_id << " not found" << endl;
    return;
  }

  auto interface = it->second;

  const double perturbation = 1.0e-7;

  // Compute finite difference sensitivity: dP/dQ
  for (int j = 0; j < num_surfaces; j++) {
    auto saved_flows = interface->current_flows_;
    interface->current_flows_[j] += perturbation;

    for (int i = 0; i < num_surfaces; i++) {
      double dP = interface->current_pressures_[i] - 
                  interface->previous_pressures_[i];
      double dQ = perturbation;
      resistance_matrix[i][j] = (dQ != 0.0) ? dP / dQ : 0.0;
    }

    interface->current_flows_ = saved_flows;
  }

  cout << "[get_resistance_matrix_1d] Resistance matrix computed" << endl;
}

/**
 * @brief Cleanup and destroy the 1D interface.
 *
 * @param problem_id The ID used to identify the 1D problem.
 */
void cleanup_1d(int problem_id) {
  auto it = OneDSolverInterface::interface_list_.find(problem_id);
  if (it == OneDSolverInterface::interface_list_.end()) {
    cerr << "Error: problem_id " << (int)problem_id << " not found" << endl;
    return;
  }

  auto interface = it->second;
  
  cout << "[cleanup_1d] Cleaning up problem_id " << (int)problem_id << endl;
  
  delete interface;
  OneDSolverInterface::interface_list_.erase(it);
  
  cout << "[cleanup_1d] Cleanup completed" << endl;
}