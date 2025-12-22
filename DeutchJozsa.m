% DEutch-Jozsa
n = 2;
sim = QCSimulator(n + 1); % 2 inputs + 1 ancilla

% Initial state: |001> (Inputs at 0, Ancilla at 1)
psi_0 = zeros(2^(n+1), 1);
psi_0(2) = 1;

% --- TEST 1: CONSTANT FUNCTION ---
f_const = @(x) 1; % Always returns 1

circuit_const = {
    {'H', 1},
    {'H', 2},
    {'H', 3},  % Step 2: Superposition
    {'ORACLE', f_const, [1, 2], 3}, % Step 3: The Oracle
    {'H', 1},
    {'H', 2}             % Step 4: Interference
};

psi_final_c = sim.Simulate(psi_0, circuit_const);
sim.displayState(psi_final_c);
[res_c, ~] = sim.Measure(psi_final_c);
fprintf('Constant Test Result (Inputs Q1, Q2): %d%d\n', res_c(1), res_c(2));

% --- TEST 2: BALANCED FUNCTION ---
f_bal = @(x) mod(x, 2); % 0 for even x, 1 for odd x

circuit_bal = {
    {'H', 1},
    {'H', 2},
    {'H', 3},
    {'ORACLE', f_bal, [1, 2], 3},
    {'H', 1},
    {'H', 2}
};

psi_final_b = sim.Simulate(psi_0, circuit_bal);
sim.displayState(psi_final_b);
[res_b, ~] = sim.Measure(psi_final_b);
fprintf('Balanced Test Result (Inputs Q1, Q2): %d%d\n', res_b(1), res_b(2));

