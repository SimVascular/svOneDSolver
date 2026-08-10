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
  * @brief Return the interface associated with a problem ID.
  * @throws cvException if the problem ID does not exist.
  */
  static OneDSolverInterface* get_interface(int problem_id);

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