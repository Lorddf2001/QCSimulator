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
disp('Bob''s recovered state (if Alice measured 00):');
psi_bob = psi_final(1:2) * 2;
for i = [1 2]
    fprintf('  |%d>: %.3f + %.3fi\n', i-1, real(psi_bob(i)), imag(psi_bob(i)));
end

psi_message = psi_start([1 5]);   % Original amplitudes of Q1
if norm(psi_message - psi_bob) < 1e-10
    fprintf('\nTeleportation succeeded! \n');
else
    fprintf('\nTeleportation failed \n');
end




