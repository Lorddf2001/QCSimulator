function testing()

n = 3;
sim = QCSimulator(n);
initial_state_idx = 3; % Let's start with |2> (binary 010)

% Create initial state
psi_start = zeros(2^n, 1);
psi_start(initial_state_idx - 1) = 1;

% 1. Forward QFT
disp('Step 1: Applying QFT...');
psi_freq = sim.Simulate(psi_start, generate_qft_instructions(n));

% 2. Inverse QFT
disp('Step 2: Applying IQFT...');
psi_reconstructed = sim.Simulate(psi_freq, generate_invqft_instructions(n));

% 3. Check Results
disp('Reconstructed State:');
sim.displayState(psi_reconstructed);

diff = norm(psi_start - psi_reconstructed);
if diff < 1e-10
    fprintf('\nSUCCESS: IQFT(QFT(|%d>)) = |%d>\n', initial_state_idx, initial_state_idx);
else
    fprintf('\nFAILURE: State was not reconstructed. Diff: %e\n', diff);
end

psi_noise = sim.ApplyNoise(psi_reconstructed, 0.1);
sim.displayState(psi_noise);


endfunction

function insts = generate_qft_instructions(n)
    insts = {};
    for i = 1:n
        % Apply Hadamard to the current qubit
        insts{end+1} = {'H', [i]};

        % Apply controlled rotations for all qubits after it
        k = 2;
        for j = i+1:n
            insts{end+1} = {'CPHASE', [j, i], k};
            k = k + 1;
        end
    end

    % Add SWAPs to reverse the order at the end
    for i = 1:floor(n/2)
        insts{end+1} = {'SWAP', [i, n-i+1]};
    end
endfunction

function insts = generate_invqft_instructions(n)
    insts = {};

    % 1. Start with SWAPs (Inverse of the end of QFT)
    for i = floor(n/2):-1:1
        insts{end+1} = {'SWAP', [i, n-i+1]};
    end

    % 2. Reverse the QFT logic
    for i = n:-1:1
        % Apply controlled rotations in reverse order
        for j = n:-1:i+1
            insts{end+1} = {'INVCPHASE', [j, i], (j-i+1)};
        end
        % Apply Hadamard
        insts{end+1} = {'H', [i]};
    end
endfunction


