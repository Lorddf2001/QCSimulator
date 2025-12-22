% NOTE:
% - This simulator works in the full Hilbert space (statevector simulation)
% - All gates are constructed as dense matrices (O(2^(2N)) memory)
% - Intended for didactic use and small numbers of qubits
classdef QCSimulator
    properties (Constant)
        H  = 1/sqrt(2) * [1, 1; 1, -1];
        X  = [0, 1; 1, 0];
        Y = [0, -1i; 1i, 0];
        Z = [1, 0; 0, -1];
        S  = [1, 0; 0, 1i];
        SDG = [1, 0; 0, -1i]; % S-Dagger (Inverse of S)
        T  = [1, 0; 0, exp(1i * pi / 4)];
        TDG = [1, 0; 0, exp(-1i * pi / 4)]; % T-Dagger (Inverse of T)
        P0 = [1, 0; 0, 0]; % Projector |0><0|
        P1 = [0, 0; 0, 1]; % Projector |1><1|
        Id = eye(2);
    endproperties

    properties
        numQubits; % Number of qubits
    endproperties

    methods
        function obj = QCSimulator(numQubits)
            obj.numQubits = numQubits;
        endfunction

        % Main function that constructs the unitary operator associated with the circuit
        % and applies it to the initial state psi_initial.
        % Instructions is a cell array describing a sequence of gates.
        % Each instruction is of the form:
        %   {'GATE', target_indices, ...}
        function psi = Simulate(obj, psi_initial, instructions)
            psi = psi_initial;

            for k = 1:length(instructions)
                inst = instructions{k};
                gateName = upper(inst{1});
                targets  = inst{2};

                if strcmp(gateName, 'ORACLE')
                  % ORACLE instruction: {'ORACLE', f_handle, input_indices, target_index}
                  f_handle = inst{2};
                  inputs = inst{3};
                  target = inst{4};
                  U = obj.build_oracle(f_handle, inputs, target);
                elseif length(targets) == 1
                    % Single Qubit Gate Logic
                    U = obj.build_single_gate(gateName, targets);
                elseif length(targets) == 2
                    % Two Qubit Gate Logic (e.g., CNOT)
                    U = obj.build_2qubit_gate(gateName, targets(1), targets(2));
                elseif length(targets) == 3
                    % CCNOT, CCZ
                    U = obj.build_3qubit_gate(gateName, targets(1), targets(2), targets(3));
                endif

                psi = U * psi;
            endfor
        endfunction

        % Helper to build U for single qubit gates
        % Build the N-qubit operator corresponding to a single-qubit gate
        % U = I ⊗ ... ⊗ G ⊗ ... ⊗ I
        function U = build_single_gate(obj, name, idx)
            % Check if the gate exists in our Constant properties
            if ~isprop(obj, name)
                error('Simulator:UnknownGate', 'Gate "%s" is not defined in QCSimulator constants.', name);
            endif

            G = obj.(name); % Dynamic access to Constant properties
            U = 1;
            for i = 1:obj.numQubits
                if i == idx
                    U = kron(U, G);
                else
                    U = kron(U, obj.Id);
                endif
            endfor
        endfunction

        % Helper to build CNOT for ANY two qubits
        function U = build_2qubit_gate(obj, name, q1, q2)
            U = zeros(2^obj.numQubits);

            % Controlled Gates: Target logic depends on name
            % General controlled gate construction
            % U = |0><0|_c ⊗ I + |1><1|_c ⊗ G_t
            if strcmp(name, 'CNOT') || strcmp(name, 'CZ')
              if strcmp(name, 'CNOT')
                targetGate = obj.X;
              elseif strcmp(name, 'CZ')
                targetGate = obj.Z;
              endif;

               % Term 1: Control is |0>, Target is Identity
               U0 = 1;
               for i = 1:obj.numQubits
                   if i == q1; U0 = kron(U0, obj.P0);
                      else U0 = kron(U0, obj.Id); endif
                  endfor

                  % Term 2: Control is |1>, Target is targetGate
                  U1 = 1;
                  for i = 1:obj.numQubits
                      if i == q1; U1 = kron(U1, obj.P1);
                      elseif i == q2; U1 = kron(U1, targetGate);
                      else U1 = kron(U1, obj.Id); endif
                  endfor
                  U = U0 + U1;

            elseif strcmp(name, 'SWAP')
                % SWAP is built by mapping 01->10 and 10->01
                % Constructed as a sum of outer products:
                % SWAP = sum_{a,b ∈ {0,1}} |b,a><a,b|
                % We iterate through all 4 basis combinations for qubits q1, q2
                % This is mathematically the cleanest way for arbitrary indices

                for b1 = 0:1
                    for b2 = 0:1
                        % In SWAP, the output bits are swapped (out1=b2, out2=b1)
                        U = U + obj.projector_2q(q1, q2, b1, b2, b2, b1);
                    endfor
                endfor
            endif
        endfunction

        % Helper to create a mapping: |b1,b2><a1,a2| at specific qubits
        function P = projector_2q(obj, q1, q2, a1, a2, b1, b2)
            P = 1;
            for i = 1:obj.numQubits
                if i == q1
                    term = zeros(2);
                    term(b1+1, a1+1) = 1; % Outer product |b1><a1|
                    P = kron(P, term);
                elseif i == q2
                    term = zeros(2);
                    term(b2+1, a2+1) = 1; % Outer product |b2><a2|
                    P = kron(P, term);
                else
                    P = kron(P, obj.Id);
                endif
            endfor
        endfunction

        function U = build_3qubit_gate(obj, gateName, c1, c2, t)
            dim = 2^obj.numQubits;
            U = eye(dim);

            for i = 0:(dim-1)
                bits = obj.dec2bits(i);

                switch gateName
                    case 'CCNOT'
                        % Flip target bit if both controls are 1
                        %  |c1,c2,t> → |c1,c2, t ⊕ (c1 ∧ c2)>
                        if bits(c1) == 1 && bits(c2) == 1
                            new_bits = bits;
                            new_bits(t) = 1 - bits(t);
                            j = obj.bits2dec(new_bits);
                            % Replace |i> → |j> in the unitary:
                            % remove identity mapping and insert permutation
                            U(i+1, i+1) = 0;
                            U(j+1, i+1) = 1;
                        endif

                    case 'CCZ'
                        % Flip phase if all three qubits are 1
                        % Note: in CCZ, t is just the 3rd control
                        % |c1,c2,t> → (-1)^{c1 c2 t} |c1,c2,t>
                        if bits(c1) == 1 && bits(c2) == 1 && bits(t) == 1
                            U(i+1, i+1) = -1;
                        endif
                endswitch
            endfor
        endfunction

        % Oracle implements:
        % |x>|y> → |x>|y ⊕ f(x)>
        % Matrix is a permutation matrix (unitary, classical reversible)
        function U = build_oracle(obj, f_handle, input_qubits, target_qubit)
          % f_handle: a function like @(x) mod(x, 2)
          % input_qubits: vector of indices, e.g., [1, 2]
          % target_qubit: scalar index, e.g., 3

          dim = 2^obj.numQubits;
          U = zeros(dim);

          for i = 0:(dim-1)
              % 1. Convert decimal index i to bits
              bits = obj.dec2bits(i);

              % 2. Extract bits for the input 'x'
              % Note: qubits are 1-indexed in our system
              x_bits = bits(input_qubits);
              x_val = obj.bits2dec(x_bits);

              % 3. Calculate f(x)
              fx = f_handle(x_val);

              % 4. Compute y XOR f(x)
              y_bit = bits(target_qubit);
              new_y_bit = mod(y_bit + fx, 2);

              % 5. Create the new bit string and convert back to decimal
              new_bits = bits;
              new_bits(target_qubit) = new_y_bit;
              j = obj.bits2dec(new_bits);

              % 6. Set the permutation in the matrix: |j><i|
              U(j+1, i+1) = 1;
          endfor
        endfunction

        % Helper methods for oracle
        % Converts decimal to a vector of bits [q1, q2, ..., qN]
        function bits = dec2bits(obj, decimal)
            str = dec2bin(decimal, obj.numQubits); %Converts a decimal integer into a string of '0's and '1's.
            bits = str - '0'; % Converts '010' string to [0, 1, 0] numeric vector
            % implicit type casting: string -> array.
            %Since '0' has an ASCII value of 48 and '1' is 49,
            %subtracting '0' from the string converts the character array
            %into a numeric array of 0s and 1s.
        endfunction

        % Converts a vector of bits back to decimal
        function decimal = bits2dec(obj, bits)
            str = char(bits + '0');
            decimal = bin2dec(str); %Converts a string of '0's and '1's back into a decimal integer.
        endfunction

        function [outcome_bits, outcome_decimal] = Measure(obj, psi)
          % 1. Calculate probabilities: P_i = |amplitude_i|^2
          probabilities = abs(psi).^2;

          % 2. Handle potential precision errors (ensure sum is exactly 1)
          %probabilities = probabilities / sum(probabilities);

          % 3. Generate a random number between 0 and 1
          r = rand();

          % 4. Determine which state was chosen (Cumulative Distribution)
          cumulative_p = 0;
          outcome_decimal = 0;
          for i = 1:length(probabilities)
              cumulative_p = cumulative_p + probabilities(i);
              if r <= cumulative_p
                  outcome_decimal = i - 1; % Convert to 0-based decimal
                  break;
              endif
          endfor

          % 5. Convert the decimal result to a bit array for the user
          outcome_bits = obj.dec2bits(outcome_decimal);
        endfunction

        function [collapsed_psi, outcome_bits] = CollapseMeasure(obj, psi, target_indices)
          % 1. Calculate probabilities for every state in the vector
          probs = abs(psi).^2;

          % 2. Pick a random state based on the probability distribution
          r = rand();
          cumulative_prob = 0;
          measured_idx = 1;
          for i = 1:length(probs)
              cumulative_prob = cumulative_prob + probs(i);
              if r <= cumulative_prob
                  measured_idx = i - 1; % Convert to 0-based decimal
                  break;
              end
          end

          % 3. Determine the bit values of the measured state for the target qubits
          full_bits = obj.dec2bits(measured_idx);
          outcome_bits = full_bits(target_indices);

          % 4. Collapse: Zero out all states that don't match the measured outcome
          collapsed_psi = zeros(size(psi));
          for i = 0:(length(psi)-1)
              current_bits = obj.dec2bits(i);
              if all(current_bits(target_indices) == outcome_bits)
                  collapsed_psi(i+1) = psi(i+1);
              end
          end

          % 5. Re-normalize the vector
          collapsed_psi = collapsed_psi / norm(collapsed_psi);
        endfunction

        function displayState(obj, psi)
          dim = length(psi);
          for i = 1:dim
              amplitude = psi(i);
              if abs(amplitude) > 1e-6 % Only show states with non-zero probability
                  bits = obj.dec2bits(i-1); % Converts Octave's 1-based index back to a 0-based decimal value
                  bit_str = num2str(bits);  % Convert the numeric array [0 1 0] to a string "0 1 0"
                  bit_str = bit_str(bit_str ~= ' '); % Remove spaces for cleaner |001>

                  % Handle complex numbers for the display
                  re = real(amplitude);
                  im = imag(amplitude);
                  if abs(im) < 1e-6
                      fprintf(' %.3f |%s> ', re, bit_str);
                  else
                      fprintf(' (%.3f + %.3fi) |%s> ', re, im, bit_str);
                  end

                  if i < dim && any(abs(psi(i+1:end)) > 1e-6)
                    fprintf('+');
                  end
              end
          end
          fprintf('\n');
        endfunction

    endmethods
endclassdef





