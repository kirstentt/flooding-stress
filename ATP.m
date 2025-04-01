function rootATP = ATP(Par,AO2,CH2O,O2shoot,O2root)

% The rate of ATP production through per mole aerobic metabolism
a = 36;

% The rate of ATP production through per mole anaerobic metabolism
b = 2;

% The rate of ATP decay
d = 1;

% The rates of aerobic and anaerobic metabolisms
[dCH2O,CH2O_metabolism,CH2O_aer_root,CH2O_anr_root] = reserves(Par,AO2,CH2O,O2shoot,O2root);

% Steady-state root ATP
% It is assumed that ATP level reaches steady state instataneously
rootATP = (a * CH2O_aer_root + b * CH2O_anr_root) / d / Par.Vr;