function ShorsFactorizer(N)
    N = 15;
    fprintf('Attempting to factor N = %d\n', N);

    % 1. Pick a random 'a' that is coprime to N
    a = 2; % For N=15, a=2 is a classic choice

    % 2. Quantum Order Finding
    % We need to find r where a^r mod N = 1
    fprintf('Running Quantum Order Finding for %d^r mod %d...\n', a, N);
    r = QuOrderFinder(a, N);

    if isempty(r)
        fprintf('Failed to find a valid period. Try again.\n');
        return;
    end

    fprintf('Quantum measurement successful. Found order r = %d\n', r);

    % 3. Classical Factor Extraction
    if mod(r, 2) ~= 0
        fprintf('Period r=%d is odd. Cannot proceed. Restart with different a.\n', r);
        return;
    end

    % Calculate x = a^(r/2) mod N
    x = mod(a^(r/2), N);
    if mod(x + 1, N) == 0
        fprintf('Trivial result (x = -1 mod N). Restart with different a.\n');
        return;
    end

    factor1 = gcd(x - 1, N);
    factor2 = gcd(x + 1, N);

    fprintf('\nSUCCESS! Factors found: %d and %d\n', factor1, factor2);
    fprintf('%d * %d = %d\n', factor1, factor2, factor1 * factor2);
endfunction

function r_candidate = QuOrderFinder(a, N)
    % Parameters for Order Finding: a^r mod N = 1
    %a = a; % The base
    N_mod = N; % The modulus

    % --- CLASSICAL GUARD ---
    common_factor = gcd(a, N_mod);
    if common_factor > 1
        fprintf('Error: a=%d and N=%d are not coprime (GCD = %d).\n', a, N_mod, common_factor);
        fprintf('In Shor''s algorithm, if GCD > 1, you already found a factor (%d)!\n', common_factor);
        fprintf('The quantum order finder requires a^r mod N = 1, which will never happen here.\n');
        r_candidate = common_factor;
        return; % Exit early
    endif
    % -----------------------

    % Modular exponentiation
    order_func = @(x) mod(a^x, N_mod);

    % Input register m: Should be large enough so 2^m > N_mod^2 for full accuracy,
    % but for simulation, 2^m > N_mod is usually enough to see peaks.
    m = 6;
    target_bits = ceil(log2(N_mod)); % Bits needed to hold values up to N
    total_qubits = m + target_bits;

    sim = QCSimulator(total_qubits);
    psi = zeros(2^total_qubits, 1); psi(1) = 1;

    % 1. Superposition of all possible exponents x
    setup = {};
    for i = 1:m
        setup{end+1} = {'H', [i]};
    endfor
    psi = sim.Simulate(psi, setup);

    % 2. Modular Exponentiation Oracle: |x>|1> -> |x>|a^x mod N>
    % Note: We start the target register at |1> because a^0 = 1
    target_indices = (m+1):total_qubits;
    psi = sim.Simulate(psi, {{'ORACLE', order_func, 1:m, target_indices}});

    % 3. Measure Target Register (Collapses input to periodic spikes)
    [psi, outcome] = sim.CollapseMeasure(psi, target_indices);
    %fprintf('Measured target register: a^x mod N = %d\n', sim.bits2dec(outcome));
    %fprintf('Collasped state\n')
    %sim.displayState(psi);

    % 4. IQFT to find the frequency of those spikes
    iqft_inst = generate_invqft_instructions(m);
    psi_final = sim.Simulate(psi, iqft_inst);

    %fprintf('\nFinal state (after inv QFT)\n')
    %sim.displayState(psi_final);

    % MARGINALIZATION
    % We ignore the target qubits and sum the probabilities for the input register only.
    probs = abs(psi_final).^2;
    N_input = 2^m;
    input_probs = zeros(N_input, 1);
    for k = 0:(length(probs)-1)
        bits = sim.dec2bits(k);
        input_val = sim.bits2dec(bits(1:m));
        input_probs(input_val + 1) = input_probs(input_val + 1) + probs(k+1);
    endfor

    % ORDER EXTRACTION
    found_peaks = find(input_probs > 0.02) - 1;
    for i = 1:length(found_peaks)
        k = found_peaks(i);
        if k == 0, continue; endif
    endfor

        % Use the CF algorithm to find the denominator r
        r_candidate = find_order_CF(k, 2^m, order_func);

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

function r_found = find_order_CF(k, N, f_handle)
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
          % 3. MULTIPLE-VERIFICATION: Test q and its small multiples
          % Sometimes CF finds a factor of the period (e.g., 3 instead of 6)
          for multiplier = 1:4
              test_r = q * multiplier;

              % Don't test beyond the register's physical capacity
              if test_r >= N, break; endif

              % Verification: Does f(0) == f(test_r)?
              if f_handle(0) == f_handle(test_r)
                  r_found = test_r;
                  return; % Success!
              endif
          endfor
        endif
      endfor
endfunction

function [p, q] = evaluate_cf(cf)
    % Initialization of the Fundamental Recurrence Relations
    % Seeds: p(-1)=0, p(0)=1; q(-1)=1, q(0)=0
    p_prev = 0; p = 1;
    q_prev = 1; q = 0;

    for i = 1:length(cf)
        a_i = cf(i);

        % Recurrence: p_n = a_n * p_{n-1} + p_{n-2}
        p_next = a_i * p + p_prev;
        % Recurrence: q_n = a_n * q_{n-1} + q_{n-2}
        q_next = a_i * q + q_prev;

        % Shift state for the next partial quotient
        p_prev = p; p = p_next;
        q_prev = q; q = q_next;
    end
endfunction
