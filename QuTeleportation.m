% Quantum teleportation
n = 3;
sim = QCSimulator(n);

% Step 0: Create a unique state to teleport on Qubit 1
% Let's give it a specific phase using S and T
psi_0 = zeros(2^n, 1);
psi_0(1) = 1; % Start at |000>
setup_message = {
    {'H', [1]},
    {'S', [1]}  % Qubit 1 is now in a specific superposition
};
psi_start = sim.Simulate(psi_0, setup_message);
disp('Initial state');
sim.displayState(psi_start);

disp('Alice message:');
for i = [1 5]
    fprintf('  |%d>: %.3f + %.3fi\n', i-1, real(psi_start(i)), imag(psi_start(i)));
end

% Step 1-4: Teleportation Protocol
teleport_circuit = {
    % --- Step 1: Create Bell Pair between Q2 and Q3 ---
    {'H', [2]},
    {'CNOT', [2, 3]},

    % --- Step 2: Alice performs Bell measurement logic on Q1 and Q2 ---
    {'CNOT', [1, 2]},
    {'H', [1]},

    % --- Step 3: Bob applies corrections based on Alice's qubits ---
    % In a simulator, we can use controlled gates to represent this
    {'CNOT', [2, 3]}, % If Q2 is 1, Bob applies X
    {'CZ', [1, 3]}    % If Q1 is 1, Bob applies Z
};

psi_final = sim.Simulate(psi_start, teleport_circuit);
fprintf('\n')
disp('Final state (no collapsing):');
sim.displayState(psi_final);

% Manually extracting the "Alice saw 00" case and re-normalizing
% Since Alice has 4 possible measurement outcomes (00, 01, 10, 11) with equal
% probability, the amplitude of the recovered state is scaled by 1/sqrt(4).
% Multiplying by 2 restores the original normalization for comparison.
disp('Bob''s recovered state (measure collapsed the state):');
psi_bob = psi_final([1 2]) * 2; %if Alice measured 00
for i = [1 2]
    fprintf('  |%d>: %.3f + %.3fi\n', i-1, real(psi_bob(i)), imag(psi_bob(i)));
end

psi_message = psi_start([1 5]);   % Original amplitudes of Q1
if norm(psi_message - psi_bob) < 1e-10
    fprintf('\nTeleportation succeeded! \n');
else
    fprintf('\nTeleportation failed \n');
end

% With collapse measure
% Alice measures her two qubits (1 and 2)
[psi_bob_collapsed, alice_results] = sim.CollapseMeasure(psi_final, [1, 2]);
fprintf('\nWith random measure\n')
fprintf('Alice measured: %d%d\n', alice_results(1), alice_results(2));

% Because of Collapse_Measure, psi_final is now a "pure" state
% where Bob has the message, and Alice's qubits are fixed.
disp('Bob''s recovered state:');
sim.displayState(psi_bob_collapsed);

if norm(psi_message - psi_bob_collapsed([1 2])) < 1e-10
    fprintf('\nTeleportation succeeded! \n');
else
    fprintf('\nTeleportation failed \n');
end


% ------------------------------------------------------------
% COMMENT:
%
% This example shows quantum teleportation in two complementary ways.
%
% 1) "No-collapse" (state-vector) analysis:
%    The simulator keeps the full quantum state without performing
%    measurements. All four Bell-measurement outcomes (00, 01, 10, 11)
%    coexist as branches of the final state vector.
%
%    By manually selecting the amplitudes corresponding to Alice
%    measuring |00> and re-normalizing, we can verify that Bob's qubit
%    contains the original message state (up to a global phase).
%
% 2) "With collapse" measurement:
%    Using CollapseMeasure simulates an actual projective measurement
%    of Alice's qubits. One outcome is randomly selected according to
%    quantum probabilities, and the state collapses accordingly.
%
%    After the collapse, Bob's qubit is left in the teleported state,
%    while Alice's qubits are fixed to the measured classical bits.
%
% Together, these two approaches illustrate both the mathematical
% structure of teleportation and its operational, measurement-based
% interpretation in a quantum circuit simulator.
% ------------------------------------------------------------


