function [PhoStomata,AO2] = StomataPhotosynthesis(Par,AO2,CH2O,O2shoot,O2root)

% Root ATP level
rootATP = ATP(Par,AO2,CH2O,O2shoot,O2root);

% Relative stomatal aperture
%Default
PhoStomata = Par.BetaStomata + (1 - Par.BetaStomata) * (rootATP^8 / (rootATP^8 + Par.Kp^8));
%To check relevance of strenght of non-linearity
%PhoStomata = Par.BetaStomata + (1 - Par.BetaStomata) * (rootATP^2 / (rootATP^2 + Par.Kp^2));
%To check relevance of non-linearity altogether
%PhoStomata = Par.BetaStomata + (1 - Par.BetaStomata) *0.35*rootATP/4.5e-3;

% Photosynthetic oxygen production rate with maximum stomatal conductance
% The inhibition from carbohydrate accumulation is introduced
max_AO2 = Par.AO2 * (Par.KCH2O / (Par.KCH2O + CH2O)); 
% Check effect of removing this feedback inhibition
%max_AO2 = Par.AO2 * 0.2;%0.1 too small; 

% Effective photosynthetic oxygen production rate
AO2 = PhoStomata * max_AO2;