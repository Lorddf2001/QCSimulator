% Quantum Period Finder
function QuPeriodFinder()

% A simple function with period r_expected: f(x) = x mod r_expected
r_expected = 6
period_func = @(x) mod(x, r_expected);

m = 5; % Input register size (can find periods up to 2^m)
% Target register size: needs enough bits to hold the result of mod(x, r_expected).
n = m + ceil(log2(r_expected));
sim = QCSimulator(n);
psi_0 = zeros(2^n, 1); psi_0(1) = 1; % Start at |000>


% --- STEP 1: CREATE SUPERPOSITION ---
% Apply Hadamard gates to all input qubits (1 to m).
% This creates a state where the register holds ALL numbers from 0 to 2^m - 1  simultaneously.
setup = {};
for i = 1:m
    setup{end+1} = {'H', [i]};
end
psi = sim.Simulate(psi_0, setup);

% --- STEP 2: APPLY PERIODIC ORACLE ---
% This computes f(x) for all x in parallel.
% It transforms |x>|0> into |x>|f(x)>.
psi = sim.Simulate(psi, {{'ORACLE', period_func, [1:m], (m+1):n}});

% --- STEP 3: MEASURE TARGET REGISTER ---
% This "collapses" the input register so it only contains 'x' values that
% produced the same f(x) result. These values are exactly one period 'r' apart.
[psi, ~] = sim.CollapseMeasure(psi, (m+1):n);
fprintf('Collasped state\n')
sim.displayState(psi);

% --- STEP 4: INVERSE QFT ---
% The Inverse Quantum Fourier Transform converts the spatial period (distance r)
% into a frequency peak in the Fourier basis.
iqft_inst = generate_invqft_instructions(m);
psi_final = sim.Simulate(psi, iqft_inst);

fprintf('\nFinal state (after inv QFT)\n')
sim.displayState(psi_final);

% 3. MARGINALIZATION
    % We ignore the target qubits and sum the probabilities for the input register only.
    probs = abs(psi_final).^2;
    N_input = 2^m;
    input_probs = zeros(N_input, 1);
    for k = 0:(length(probs)-1)
        bits = sim.dec2bits(k);
        input_val = sim.bits2dec(bits(1:m));
        input_probs(input_val + 1) = input_probs(input_val + 1) + probs(k+1);
    end

    % --- PLOTTING ---
    % Visualizes the "frequency" peaks.
    figure;
    bar(0:N_input-1, input_probs);
    title(['Input Register Distribution (r=', num2str(r_expected), ')']);
    xlabel('Measured Value (k)');
    ylabel('Probability');
    grid on;

    % 4. PERIOD EXTRACTION
    % We find peaks where probability > 2%.
    found_peaks = find(input_probs > 0.02) - 1;
    fprintf('\nFound %d significant peaks.\n', length(found_peaks));
    fprintf('\nSignificant peaks at: %s\n', mat2str(found_peaks));

    for i = 1:length(found_peaks)
        k_val = found_peaks(i);
        if k_val == 0, continue; end % k=0 doesn't help find the period.

        % find_period_CF uses the Continued Fractions algorithm to guess
        % the simplest fraction p/r that matches the measured k/N.
        r_candidate = find_period_CF(k_val, N_input, period_func);

        if ~isempty(r_candidate)
            fprintf('Peak at %d gives ratio %d/%d -> Period r = %d\n', k_val, k_val, N_input, r_candidate);
            break;
        end
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

function r_found = find_period_CF(k, N, f_handle)
    % k: measured peak value
    % N: total states in input register (2^m)
    % f_handle: the periodic function to verify the candidate

    if k == 0
        r_found = []; return;
    end

    % 1. Get Continued Fraction Coefficients
    cf = [];
    val = k / N;
    for i = 1:10 % Limit iterations
        a = floor(val);
        cf(end+1) = a;
        rem = val - a;
        if rem < 1e-6, break; end
        val = 1 / rem;
    end

    % 2. Calculate Convergents (p/q)
    r_found = [];
    for i = 1:length(cf)
        % Evaluate the continued fraction up to index i
        sub_cf = cf(1:i);
        [p, q] = evaluate_cf(sub_cf);

        % Candidate period is the denominator q
        % q < N: A period cannot be larger than the number of states we measured
        if q > 0 && q < N
            % 3. VERIFICATION: Does f(0) == f(q)?
            if f_handle(0) == f_handle(q)
                r_found = q;
                return; % Found the period!
            end
        end
    end
end

function [p, q] = evaluate_cf(cf)
    % Helper to turn coefficients into a fraction p/q
    p_prev = 0; p = 1;
    q_prev = 1; q = 0;

    for i = 1:length(cf)
        a = cf(i);
        % The recurrence formula: Current = (Quotient * Previous) + One_Before_Previous
        p_next = a * p + p_prev;  % Calculate the new numerator using the previous two
        q_next = a * q + q_prev;  % Calculate the new denominator using the previous two

        p_prev = p; p = p_next;
        q_prev = q; q = q_next;
    end
    % Result is p/q (note: p and q indices shift in this standard iterative algorithm)
    % Correcting for the specific iteration logic:
    p = p; q = q;
end

