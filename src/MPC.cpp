#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 10;
double dt = 0.1;

// reference velocity
double ref_v = 75;

// offsets for DD vector
size_t x_s = 0;
size_t y_s = x_s + N;
size_t psi_s = y_s + N;
size_t v_s = psi_s + N;
size_t cte_s = v_s + N;
size_t epsi_s = cte_s + N;
size_t delta_s = epsi_s + N;
size_t a_s = delta_s + N - 1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
    std::vector<double> factors;
public:
    // Fitted polynomial coefficients
    Eigen::VectorXd coeffs;

    FG_eval(Eigen::VectorXd coeffs, std::vector<double> factors) {
        this->coeffs = coeffs;
        this->factors = factors;
    }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector &fg, const ADvector &vars) {
        // TODO: implement MPC
        // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
        // NOTE: You'll probably go back and forth between this function and
        // the Solver function below.

        double cteFactor = factors[0];
        double epsiFactor = factors[1];
        double steerFactor = factors[2];
        double accFactor = factors[3];
        double sequentialSteerFactor = factors[4];
        double sequentialAccFactor = factors[5];

        std::cout << "steerFactor: " << steerFactor << std::endl;
        // initialize cost
        fg[0] = 0;

        // error costs (cte, angle, speed)
        for (size_t t = 0; t < N; t++) {
            fg[0] += cteFactor * CppAD::pow(vars[cte_s + t], 2);
            fg[0] += epsiFactor * CppAD::pow(vars[epsi_s + t], 2);
            fg[0] += CppAD::pow(vars[v_s + t] - ref_v, 2); // speed deviation
        }

        // Minimize the use of actuators.
        for (size_t t = 0; t < N - 1; t++) {
            fg[0] += steerFactor * CppAD::pow(vars[delta_s + t], 2);
            fg[0] += accFactor * CppAD::pow(vars[a_s + t], 2);
        }

        // Minimize the value gap between sequential actuations.
        for (size_t t = 0; t < N - 2; t++) {
            fg[0] += sequentialSteerFactor * CppAD::pow(vars[delta_s + t + 1] - vars[delta_s + t], 2);
            fg[0] += sequentialAccFactor * CppAD::pow(vars[a_s + t + 1] - vars[a_s + t], 2);
        }

        // setup constraints (all with offset 1 for the cost)
        fg[1 + x_s] = vars[x_s];
        fg[1 + y_s] = vars[y_s];
        fg[1 + psi_s] = vars[psi_s];
        fg[1 + v_s] = vars[v_s];
        fg[1 + cte_s] = vars[cte_s];
        fg[1 + epsi_s] = vars[epsi_s];

        for (size_t t = 1; t < N; ++t) {
            AD<double> x1 = vars[x_s + t];
            AD<double> x0 = vars[x_s + t - 1];
            AD<double> y1 = vars[y_s + t];
            AD<double> y0 = vars[y_s + t - 1];
            AD<double> psi1 = vars[psi_s + t];
            AD<double> psi0 = vars[psi_s + t - 1];
            AD<double> v1 = vars[v_s + t];
            AD<double> v0 = vars[v_s + t - 1];
            AD<double> cte1 = vars[cte_s + t];
            AD<double> cte0 = vars[cte_s + t - 1];
            AD<double> epsi1 = vars[epsi_s + t];
            AD<double> epsi0 = vars[epsi_s + t - 1];

            AD<double> delta0 = vars[delta_s + t - 1];
            AD<double> a0 = vars[a_s + t - 1];

            AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
            AD<double> f0prime = coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2);
            AD<double> psides0 = CppAD::atan(f0prime);

            fg[1 + x_s + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
            fg[1 + y_s + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
            fg[1 + psi_s + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
            fg[1 + v_s + t] = v1 - (v0 + a0 * dt);
            fg[1 + cte_s + t] = cte1 - ((f0 - y0) + v0 * CppAD::sin(epsi0) * dt);
            fg[1 + epsi_s + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
        }

    }
};

//
// MPC class definition implementation.
//
MPC::MPC(std::vector<double> factors) { this->factors = factors; }

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
    bool ok = true;
    typedef CPPAD_TESTVECTOR(double) Dvector;

    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3]; // in miles/hour  * 0.44704 to m/s
    double cte = state[4];
    double epsi = state[5];

    double steeringBounds_rad = 0.436332;
    // Set the number of model variables (includes both states and inputs).
    // For example: If the state is a 4 element vector, the actuators is a 2
    // element vector and there are 10 timesteps. The number of variables is:
    // 4 * 10 + 2 * 9
    size_t n_vars = N * 6 + (N - 1) * 2;
    // Set the number of constraints
    size_t n_constraints = N * 6;

    // Initial value of the independent variables.
    // SHOULD BE 0 besides initial state.
    Dvector vars(n_vars);
    for (size_t i = 0; i < n_vars; i++) {
        vars[i] = 0;
    }
    vars[x_s] = x;
    vars[y_s] = y;
    vars[psi_s] = psi;
    vars[v_s] = v;
    vars[cte_s] = cte;
    vars[epsi_s] = epsi;

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    // Set lower and upper limits for variables.
    for (size_t i = 0; i < delta_s; ++i) {
        vars_lowerbound[i] = -1e19;
        vars_upperbound[i] = 1e19;
    }
    for (size_t i = delta_s; i < a_s; ++i) {
        vars_lowerbound[i] = -steeringBounds_rad;
        vars_upperbound[i] = steeringBounds_rad;
    }
    for (size_t i = a_s; i < n_vars; ++i) {
        vars_lowerbound[i] = -1;
        vars_upperbound[i] = 1;
    }


    // Lower and upper limits for the constraints
    // Should be 0 besides initial state.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    constraints_lowerbound[x_s] = x;
    constraints_upperbound[x_s] = x;
    constraints_lowerbound[y_s] = y;
    constraints_upperbound[y_s] = y;
    constraints_lowerbound[psi_s] = psi;
    constraints_upperbound[psi_s] = psi;
    constraints_lowerbound[v_s] = v;
    constraints_upperbound[v_s] = v;
    constraints_lowerbound[cte_s] = cte;
    constraints_upperbound[cte_s] = cte;
    constraints_lowerbound[epsi_s] = epsi;
    constraints_upperbound[epsi_s] = epsi;

    // object that computes objective and constraints
    FG_eval fg_eval(coeffs, factors);

    //
    // NOTE: You don't have to worry about these options
    //
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    // of sparse routines, this makes the computation MUCH FASTER. If you
    // can uncomment 1 of these and see if it makes a difference or not but
    // if you uncomment both the computation time should go up in orders of
    // magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time          0.5\n";

    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;

    // solve the problem
    CppAD::ipopt::solve<Dvector, FG_eval>(
            options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
            constraints_upperbound, fg_eval, solution);

    // Check some of the solution values
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    if (!ok) {
        std::cout << "WARNING: could not optimize values!" << std::endl;
    }

    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    // Return the first actuator values. The variables can be accessed with
    // `solution.x[i]`.
    //
    // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
    // creates a 2 element double vector.

    std::vector<double> result =
            {
                    solution.x[delta_s],
                    solution.x[a_s],
                    cost
            };

    // predicted values
    for (size_t i = 0; i < N - 1; ++i) {
        result.push_back(solution.x[1 + x_s + i]);
        result.push_back(solution.x[1 + y_s + i]);

    }

    return result;
}
