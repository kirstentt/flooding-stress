function dROLB = ROL_barrier(Par,O2soil)

% Superimposed ROL barrier induction rate
% With anoxia the ROL barrier would be completely induced around 24 hours.
Alpha_ROLB = 1e-5; 

% Rhizosphere oxygen deficit compared to the thresold value
% At which ROL barrier induction rate is half of the maximum
K_ROLB = 0.001;

% ROL barrier induction rate
% When rhizosphere oxygen level is beyond the thresold, no ROL barrier induction
dROLB = Alpha_ROLB * max(Par.ROLB_induction - O2soil,0)^2 / (max(Par.ROLB_induction - O2soil,0)^2 + K_ROLB^2);