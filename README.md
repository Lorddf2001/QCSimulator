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
% Initialize a 3-qubit simulator
numQubit = 3;
sim = QCSimulator(numQubit);

% Define a simple circuit
circuit = {
    {'X', 2},            % Flip the 2nd qubit
    {'CCZ', [1, 2, 3]},  % Flip phase of |111>
    {'X', 2}             % Flip it back
};

% Initialize state |000>
psi = zeros(2^numQubit, 1);
psi(1) = 1;

% Execute the circuit
final_psi = sim.Simulate(psi, circuit);
sim.displayState(final_psi);

% Measure
[bits, ~] = sim.Measure(final_psi);
found_dec = sim.bits2dec(bits(1:numQubit));
disp(['Measured value (decimal): ', num2str(found_dec)]);

