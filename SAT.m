% SAT Solver (NP-complete problem) using Grover
n_inputs = 3;
n_total = 4;
sim = QCSimulator(n_total);

% Define the 3-SAT logic function
% x is the decimal representation of [v1, v2, v3]
% bitget(x,1) -> v3 (LSB)
% bitget(x,2) -> v2
% bitget(x,3) -> v1 (MSB)
%
% Clauses:
% Clause 1: ( v1 OR  v2 OR  v3)
% Clause 2: ( v1 OR  v2 OR !v3)
% Clause 3: ( v1 OR !v2 OR  v3)
% Clause 4: (!v1 OR  v2 OR  v3)
% Clause 5: (!v1 OR !v2 OR  v3)
% Clause 6: (!v1 OR  v2 OR !v3)
% Clause 7: ( v1 OR !v2 OR !v3)
% Note: Because anonymous functions in Octave are limited to one line,
% we write the conjunction compactly:
f_sat = @(x) ( ...
    ( bitget(x,3) ||  bitget(x,2) ||  bitget(x,1)) && ... % Clause 1
    ( bitget(x,3) ||  bitget(x,2) || !bitget(x,1)) && ... % Clause 2
    ( bitget(x,3) || !bitget(x,2) ||  bitget(x,1)) && ... % Clause 3
    (!bitget(x,3) ||  bitget(x,2) ||  bitget(x,1)) && ... % Clause 4
    (!bitget(x,3) || !bitget(x,2) ||  bitget(x,1)) && ... % Clause 5
    (!bitget(x,3) ||  bitget(x,2) || !bitget(x,1)) && ... % Clause 6
    ( bitget(x,3) || !bitget(x,2) || !bitget(x,1)) ...    % Clause 7
);

% SAT Truth Table
fprintf('--- SAT Truth Table ---\n');
fprintf('Dec | v1 v2 v3 | Satisfied?\n');
fprintf('---------------------------\n');
for x = 0:7
    % Check the logic for this decimal input
    satisfied = f_sat(x);

    % Convert to binary string for display
    b_str = dec2bin(x, 3);

    if satisfied
        fprintf('%d   | %s    | YES (Target)\n', x, b_str);
    else
        fprintf('%d   | %s    | No\n', x, b_str);
    end
end
fprintf('---------------------------\n\n');

% --- Grover Setup Circuit ---
fprintf('Grover\n');
psi = zeros(2^n_total, 1);
psi(2) = 1; % |0001>

grover_step = {
    {'ORACLE', f_sat, [1, 2, 3], 4},
    {'H', 1},
    {'H', 2},
    {'H', 3},
    {'ORACLE', @(y) (y==0), [1, 2, 3], 4}, %2|0><0| - I
    {'H', 1},
    {'H', 2},
    {'H', 3}
};

circuit = [
    {{'H', 1},
    {'H', 2},
    {'H', 3},
    {'H', 4}},
    grover_step,
    grover_step,
    {{'H', 4}}
];

% Execute
final_psi = sim.Simulate(psi, circuit);
sim.displayState(final_psi);

[bits, ~] = sim.Measure(final_psi);
solution = bits(1:3);
if f_sat(sim.bits2dec(solution)) == 1
  fprintf('Found SAT Assignment (v1, v2, v3): %s\n', num2str(solution));
else
  fprintf('Grover failed\n');
end
