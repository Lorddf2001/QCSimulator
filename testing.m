%Script for testing new additions
n = 3;
sim = QCSimulator(n);

psi_0 = zeros(2^(n), 1);
psi_0(1) = 1;

circuit = {
    {'H', 1},  % Put Qubit 1 into (|0> + |1>)/sqrt(2)
    {'H', 2},  % Put Qubit 2 into (|0> + |1>)/sqrt(2)
    {'S', 1},  % Add a 90-degree phase to the |1> part of Qubit 1
    {'T', 2}   % Add a 45-degree phase to the |1> part of Qubit 2
};

psi_final = sim.Simulate(psi_0, circuit);
sim.displayState(psi_final);

% Test: Two T-gates should equal one S-gate
circuit_T_squared = {
    {'H', [1]}
    {'T', [1]},
    {'T', [1]}
};

circuit_S = {
    {'H', [1]}
    {'S', [1]}
};

psi_t2 = sim.Simulate(psi_0, circuit_T_squared);
sim.displayState(psi_t2);

psi_s = sim.Simulate(psi_0, circuit_S);
sim.displayState(psi_s);
