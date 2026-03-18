/**
 * @file interface.h
 * @brief svOneDSolver callable interface for 3D-1D coupling.
 */

#include <map>
#include <string>
#include <vector>

/**
 * @brief Interface class for calling svOneD from external programs
 */
class OneDSolverInterface {
 public:
  /**
   * @brief Construct a new interface object
   * @param input_file_name The 1D JSON file which specifies the model
   */
  OneDSolverInterface(const std::string& input_file_name);

  /**
   * @brief Destroy the interface object
   */
  ~OneDSolverInterface();

  /**
   * @brief Counter for the number of interfaces
   */
  static int problem_id_count_;

  /**
   * @brief List of interfaces
   */
  static std::map<int, OneDSolverInterface*> interface_list_;

  /**
   * @brief ID of current interface
   */
  int problem_id_ = 0;

  /**
   * @brief 1D input (JSON) file
   */
  std::string input_file_name_;

  /**
   * @brief Time step size of the external program (3D solver)
   *
   * This is required for coupling with a 3D solver
   */
  double external_step_size_ = 0.001;

  /**
   * @brief Number of 1D segments in the network
   */
  int num_segments_ = 0;

  /**
   * @brief Number of nodes in the network
   */
  int num_nodes_ = 0;

  /**
   * @brief Number of surfaces coupled with 3D
   */
  int num_coupled_surfaces_ = 0;

  /**
   * @brief Names of coupled segments
   */
  std::vector<std::string> coupled_segment_names_;

  /**
   * @brief Type of coupling for each surface ("DIR" or "NEU")
   */
  std::vector<std::string> coupling_types_;

  /**
   * @brief In/out signs for coupled surfaces
   */
  std::vector<double> in_out_signs_;

  /**
   * @brief IDs of coupled segments
   */
  std::vector<int> coupled_segment_ids_;

  /**
   * @brief Current flows at coupled surfaces (mL/s)
   */
  std::vector<double> current_flows_;

  /**
   * @brief Previous flows at coupled surfaces (mL/s)
   */
  std::vector<double> previous_flows_;

  /**
   * @brief Current pressures at coupled surfaces (mmHg)
   */
  std::vector<double> current_pressures_;

  /**
   * @brief Previous pressures at coupled surfaces (mmHg)
   */
  std::vector<double> previous_pressures_;

  /**
   * @brief Resistance matrix (sensitivity) dP/dQ for coupling
   */
  std::vector<std::vector<double>> resistance_matrix_;

  /**
   * @brief Current time step
   */
  int time_step_ = 0;

  /**
   * @brief Current simulation time (seconds)
   */
  double current_time_ = 0.0;

  std::string model_name_;
  double time_step_size_ = 0.0;
  int max_step_ = 0;

  // Coupling options
  std::string coupling_status_ = "OFF";     // "ON" or "OFF"
  std::string coupling_type_;       // "DIR" or "NEU"
  int coupling_substeps_ = 10;       // number of substeps
};