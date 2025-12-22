n = 3;
sim = QCSimulator(n);

psi_0 = zeros(2^(n), 1);
psi_0(1) = 1;

circuit = {

    {'CCZ', [1, 2, 3]},   % Mark |000> (which is now |111> thanks to X), equivalent to 2

};

psi_final = sim.Simulate(psi_0, circuit);
sim.displayState(psi_final);
