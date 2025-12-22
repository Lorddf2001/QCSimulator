# QCSimulator-Octave ⚛️

A **didactic quantum circuit simulator** written entirely in **GNU Octave**.  

This project is designed for **clarity and learning**, not for speed or handling large numbers of qubits. It helps users understand how state vectors and unitary matrices interact to implement quantum algorithms.

## 🚀 Key Features

* **Object-Oriented Core:** Quantum register managed via the main `QCSimulator` class.  
* **Universal Gate Support:** Single-qubit gates (H, X, Y, Z, T, S) and multi-qubit gates (CNOT, Toffoli, CCZ, etc.).  
* **State Inspection:** Methods to inspect the state vector and measurement probabilities.  
* **Ready-to-Use Examples:** Scripts demonstrating foundational quantum algorithms.

## 📂 Project Structure

* `QCSimulator.m` – Main class handling simulator logic (initialization, gate application, Kronecker products).  
* `DeutschJozsa.m` – Implementation of the Deutsch-Jozsa algorithm.  
* `Grover.m` – Implementation of Grover’s search algorithm.  
* `SAT.m` – Solving Boolean satisfiability problems using quantum oracles.  

## 🛠️ Getting Started

1. Make sure [GNU Octave](https://www.gnu.org/software/octave/) is installed.  
2. Clone the repository:  
    ```bash
    git clone https://github.com/your-username/QCSimulator-Octave.git
    ```  
3. Open Octave in the project folder and run one of the example scripts.

## 💻 Usage Example

```octave
% ========================================
% Example: Quantum Circuit Simulation
% 3 Variable Qubits + 1 Ancilla
% ========================================

%% 1. Initialize Simulator
n = 4;                       % Total qubits (3 variables + 1 ancilla)
sim = QCSimulator(n);         % Create simulator instance
psi_initial = zeros(2^n, 1);  % Start in |0000>
psi_initial(1) = 1;

%% 2. Define Oracle Function
% Marks the state |101> (decimal 5) as the solution
f = @(x) (x == 5);

%% 3. Construct the Circuit
instructions = {
    % --- Single Qubit Gates ---
    {'H', [1]},   % Put Q1 in superposition
    {'H', [2]},   % Put Q2 in superposition
    {'H', [3]},   % Put Q3 in superposition
    {'T', [1]},   % Apply a 45-degree phase shift on Q1

    % --- Two-Qubit Gate ---
    {'CNOT', [1, 2]},   % Entangle Q1 and Q2

    % --- Three-Qubit Gate (Toffoli) ---
    {'CCNOT', [1, 2, 3]},   % Flip Q3 if Q1 and Q2 are |1>

    % --- Functional Oracle ---
    % Marks the solution onto the Ancilla (Q4)
    {'ORACLE', f, [1, 2, 3], 4},

    % --- Triple Qubit Phase Gate ---
    {'CCZ', [1, 2, 3]}    % Flip phase if Q1, Q2, Q3 are |1>
};

%% 4. Run the Simulation
psi_final = sim.Simulate(psi_initial, instructions);

%% 5. Analyze Results
disp('Final Quantum State:');
sim.displayState(psi_final);

%% 6. (Optional) Measure Qubits
[bits, ~] = sim.Measure(psi_final);
found_dec = sim.bits2dec(bits(1:3));
disp(['Measured value (decimal, variable qubits): ', num2str(found_dec)]);