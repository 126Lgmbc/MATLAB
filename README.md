## Introduction
This project contains the MATLAB implementation of the Adaptive Constraint-balancing Multi-objective Differential Evolution (ACB-MODE) Algorithm for solving Vehicle Routing Problems Considering Crowdsourced Delivery (VRP-C). The algorithm specifically addresses the complex bi-objective optimization problem of minimizing total travel distance while minimizing total driver payment, considering constraints such as vehicle capacity, time windows, and driver availability.

The ACB-MODE algorithm incorporates several innovative features including a Three-stage Progressive Initialization Strategy, Adaptive Constraint Balancing Mechanism (ACH) and Refined Local Search. The algorithm is designed to efficiently explore the Pareto front while maintaining a balance between feasible and infeasible solutions throughout the optimization process.


## File Descriptions

### 1. **`ACB_MODE_main.m`**
- **Purpose**:  Implements the main logic of the ACB-MODE algorithm.
- **Key Features**:
  - Reads Solomon-format VRP instances and initializes problem parameters.
  - Implements a three-stage initialization strategy (grouping-clustering-route building).
  - Executes the main adaptive differential evolution loop with dynamic parameter adjustment.
  - Integrates Adaptive Constraint Handling (ACH) for epsilon-based constraint management.
  - Performs hybrid DE operations with continuous mapping and decoding mechanisms.
  - Triggers local search (LNS) with adaptive probability.
  - Manages the archive of solutions and performs periodic evaluation.
  - Visualizes final results including Pareto front and optimal routes.

### 2. **`SolutionOperators.m`**
- **Purpose**: Implements solution generation, manipulation, and optimization operators.
- **Key Features**:
  - Three-stage route construction (`grouping_stage`, `cluster_stage`, `threeStage_build_routes`).
  - Hybrid differential evolution operators for solution generation (`hybridDE_and_decode`).
  - Continuous-discrete encoding mapping functions (`map_encoding_to_cont`, `decode_cont_to_encoding`).
  - Large Neighborhood Search (LNS) implementation (`local_search_triggered`).
  - Constraint repair functions (`ensureDriverConstraints`).
  - Route perturbation operators (`perturb_encoding`).

### 3. **`ACH_Selection.m`**
- **Purpose**: Implements adaptive constraint handling and selection mechanisms.
- **Key Features**:
  - ACH parameter update based on feasible solution ratio (`updateACH`).
  - ε-constraint handling technique for balancing feasible and infeasible solutions.
  - Multi-criteria selection combining feasibility, dominance, and diversity (`selection_ACB_MODE`).
  - Neighborhood selection using R(p|q) concept (`neighborhood_by_R`).
  - Parent selection mechanisms for differential evolution (`select_parents_from_neighborhood`).

### 4. **`EvaluationFunctions.m`**
- **Purpose**: Evaluates solution quality and checks constraint violations.
- **Key Features**:
  - Single solution evaluation (`calObj`) computing both objectives and constraint violations.
  - Batch evaluation of population (`calObjs`).
  - Route constraint checking functions (`checkConstraints_route`, `checkRouteTimeFeasible`).
  - Improved Clarke-Wright algorithm for initial route construction (`improvedCW_cluster`).

### 5. **`Utils.m`**
- **Purpose**: Provides utility functions for multi-objective optimization and visualization.
- **Key Features**:
  - Non-dominated sorting and crowding distance calculation (`nondominated_fronts`, `crowding_distance`).
  - Performance indicators computation (Hypervolume `calculateHV`, Distribution `calculateE`).
  - Visualization functions for routes and Pareto front (`drawRoutes`).
  - General utility functions for the optimization process.


## Usage Instructions
- Ensure MATLAB is installed and properly configured on your system.
- Input data must follow Solomon format with columns: CUST NO., XCOORD., YCOORD., DEMAND, READY TIME, DUE DATE, SERVICE TIME
- The algorithm handles both full-time (closed routes) and part-time (open routes) drivers with different constraints
- Place all program files (`.m` files) in the same working directory along with the Solomon-format VRP instance file (e.g., `rc101.txt`, `rc102.txt`).
- Run the algorithm by executing `ACB_MODE_main.m` in MATLAB. Adjust algorithm parameters in the "PARAMETERS" section of the main script as needed.
