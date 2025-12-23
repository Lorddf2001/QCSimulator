%Quantum Fourier Transform
function QFT()
% 3 qubit example
n = 3;
sim = QCSimulator(n);

% Start with state |2> (binary 010)
psi_in = zeros(2^n, 1);
psi_in(3) = 1;
%sim.displayState(psi_in);

qft_instructions = {
    % Qubit 1: H and rotations
    {'H', [1]},
    {'CPHASE', [2, 1], 2}, % Rotates by 2pi / 2^2 (90°)
    {'CPHASE', [3, 1], 3}, % Rotates by 2pi / 2^3 (45°)

    % Qubit 2: H and rotations
    {'H', [2]},
    {'CPHASE', [3, 2], 2},

    % Qubit 3: H
    {'H', [3]},

    % Reversal: QFT outputs bits in reverse order
    {'SWAP', [1, 3]}
};

psi_out = sim.Simulate(psi_in, qft_instructions);
disp('QFT Result (Frequency Domain representation of |010>):');
sim.displayState(psi_out);

% 3. Verification QFT
success = verifyQFT(psi_out, 2);
if success
    fprintf('Success! All phases match the theoretical QFT formula.\n\n');
else
    fprintf('QFT failed.\n\n');
end

%Let's try bigger (and automated)
n = 5;
sim = QCSimulator(n);
x_input = 2; % Input state |2>

psi_in = zeros(2^n, 1);
psi_in(x_input + 1) = 1;

% Use your automated generator
qft_instructions = generate_qft_instructions(n);

psi_out = sim.Simulate(psi_in, qft_instructions);

disp('QFT Result:');
sim.displayState(psi_out);

% Call the verification function
success = verifyQFT(psi_out, x_input);

if success
    fprintf('Success! All phases match the theoretical QFT formula.\n');
else
    fprintf('QFT failed.\n');
end

endfunction

% --- Helper Functions Section ---

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

function success = verifyQFT(psi, x)
  N = length(psi);
  fprintf('\n--- QFT Phase Verification (Input |%d>) ---\n', x);
  success = true;

  for j = 0:N-1
      % Theoretical amplitude: (1/sqrt(N)) * exp(2*pi*i * x * j / N)
      expected_phase = (1/sqrt(N)) * exp(2j * pi * x * j / N);
      actual_phase = psi(j+1);

      % Check difference
      if abs(expected_phase - actual_phase) > 1e-10
          fprintf('Mismatch at state |%d>!\n', j);
          success = false;
      end
  end
endfunction















