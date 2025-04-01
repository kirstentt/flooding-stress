function [Mshoot,Mroot] = OxygenMetabolism(Par,CH2O,O2shoot,O2root)

% Baseline respiratory rate that does not consume carbohydrates
Alpha = 0.2; 

% The effect of oxygen and carbohydrate levels on aerobic respiratory rate
% Default
Beta_O2shoot = (O2shoot^8 / (Par.hO2shoot^8 + O2shoot^8)) * (Alpha + (CH2O*Par.M) * (1 - Alpha));
Beta_O2root  = (O2root^8 / (Par.hO2root^8 + O2root^8)) * (Alpha + (CH2O*Par.M) * (1 - Alpha));
% To check relevance of strenght of non-linearity
%Beta_O2shoot = (O2shoot^2 / (Par.hO2shoot^2 + O2shoot^2)) * (Alpha + (CH2O*Par.M) * (1 - Alpha));
%Beta_O2root  = (O2root^2 / (Par.hO2root^2 + O2root^2)) * (Alpha + (CH2O*Par.M) * (1 - Alpha));

% The calculation of effective aerobic respiratory rate in shoot and root
Mshoot    = Par.mO2shoot * Beta_O2shoot * Par.shootWeight;
Mroot     = Par.mO2root * Beta_O2root * Par.rootWeight;

end