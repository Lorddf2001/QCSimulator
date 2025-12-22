% Grover's Search
n_inputs = 3;
n_total = n_inputs + 1;
sim = QCSimulator(n_total);

% 1. Setup: Inputs at |000>, Ancilla at |1>
psi = zeros(2^n_total, 1);
psi(2) = 1; % Index 2 is |0001>

% 2. Target: Find 5 (binary 101)
target_value = 5;
f_find = @(x) (x == target_value);

% 3. The Grover Circuit
% We define one "Step" to repeat
grover_step = {
    {'ORACLE', f_find, [1, 2, 3], 4}, % Mark target, Grover oracle
    {'H', 1},
    {'H', 2},
    {'H', 3},
    {'ORACLE', @(x) (x==0), [1, 2, 3], 4}, % Equivalent to 2|0><0| - I
    {'H', 1},
    {'H', 2},
    {'H', 3}
};

% Complete Circuit
% theta = arcsin(sqrt(1/8))) = 0.361 (20.7°) ->
% Optimal iterations = nearest_int[ pi/4*theta - 1/2 ] = 2 ->
% theta_final = 103.5°, P = sin^2(theta_final) = 0.95
full_circuit = [
    {{'H', 1},
    {'H', 2},
    {'H', 3},
    {'H', 4}}, % Init Superposition
    grover_step, % Iteration 1
    grover_step, % Iteration 2
    {{'H', 4}}   % Turn Ancilla back to |0> or |1>
];

% Execute
final_psi = sim.Simulate(psi, full_circuit);

fprintf('Simple Grover\n');
% 4. Look at the result
sim.displayState(final_psi);

[bits, ~] = sim.Measure(final_psi);
found_dec = sim.bits2dec(bits(1:n_inputs));
if f_find(found_dec) == 1
  fprintf('Search Result: %d (Binary: %s)\n', found_dec, num2str(bits(1:n_inputs)));
else
  fprintf('Grover failed\n');
end


% HARD CODED GROVER
% Gate-based Oracle for finding |101>
% We don't explicitly build the oracle, we don't use a for-loop or an if-statement here!

n_inputs = 3;
n_total = 3; % No ancilla needed if we use a CCZ gate!
sim = QCSimulator(n_total);

% --- 1. THE ORACLE (Marks |101>) ---
% We use X gates to make the target state look like |111> to the CCZ
oracle_gates = {
    {'X', 2},           % Flip the 0 to a 1
    {'CCZ', [1, 2, 3]},   % Flip phase of |111>
    {'X', 2}            % Flip it back
};

% --- 2. THE DIFFUSION (Reflection around Mean) ---
diffusion_gates = {
    {'H', 1},
    {'H', 2},
    {'H', 3},
    % 2|0><0| - I
    % ------------------
    {'X', 1},
    {'X', 2},
    {'X', 3},
    {'CCZ', [1, 2, 3]},   % Mark |000> (which is now |111> thanks to X), equivalent to 2
    {'X', 1},
    {'X', 2},
    {'X', 3},
    % -------------------
    {'H', 1},
    {'H', 2},
    {'H', 3}
};

% --- 3. FULL CIRCUIT ---
% Start in superposition
psi = zeros(2^n_total, 1);
psi(1) = 1; % Start at |000>

full_circuit = [
    {{'H', 1},
    {'H', 2},
    {'H', 3}}, % Initial Superposition
    oracle_gates,
    diffusion_gates,  % Iteration 1
    oracle_gates,
    diffusion_gates   % Iteration 2
];

% Execute
final_psi = sim.Simulate(psi, full_circuit);

fprintf('Hard coded Grover\n');
sim.displayState(final_psi);

[bits, ~] = sim.Measure(final_psi);
found_dec = sim.bits2dec(bits(1:n_inputs));
if f_find(found_dec) == 1
  fprintf('Search Result: %d (Binary: %s)\n', found_dec, num2str(bits(1:n_inputs)));
else
  fprintf('Grover failed\n');
end












