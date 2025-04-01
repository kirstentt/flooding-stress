function [dCH2O,CH2O_metabolism,CH2O_aer_root,CH2O_anr_root] = reserves(Par,AO2,CH2O,O2shoot,O2root)

% Aerobic respiratory rate
[Mshoot,Mroot] = OxygenMetabolism(Par,CH2O,O2shoot,O2root);

% Convert rate of photosynthetic O2 production to glucose production
CH2O_photosynthesis = AO2 / 6;

% Convert rate of respiratory O2 consumption to glucose consumption
CH2O_aer_shoot = Mshoot / 6;
CH2O_aer_root = Mroot / 6;

% The proportion of non-NSC related anaerobic metabolic rate
Alpha = 0.2; 

% Anaerobic metabolic rate in shoot
% To generate the same amount of ATP
% Carbohydrate consumption through anaerobic metabolism is 18 folds as that through aerobic metabolism
% Default
CH2O_anr_shoot = ((Par.mO2shoot / 6 * 18) * ((Par.hO2shoot^8) / (Par.hO2shoot^8 + O2shoot^8)) ...
   * Alpha + (Par.mO2shoot / 6 * 18) * ((Par.hO2shoot^8) / (Par.hO2shoot^8 + O2shoot^8)) ...
   * (1-Alpha) * (CH2O*Par.M)) * Par.shootWeight;
% To check relevance of strenght of non-linearity
%CH2O_anr_shoot = ((Par.mO2shoot / 6 * 18) * ((Par.hO2shoot^2) / (Par.hO2shoot^2 + O2shoot^2)) ...
 %   * Alpha + (Par.mO2shoot / 6 * 18) * ((Par.hO2shoot^2) / (Par.hO2shoot^8 + O2shoot^2)) ...
 %   * (1-Alpha) * (CH2O*Par.M)) * Par.shootWeight;


% Anaerobic metabolic rate in root
% Default
CH2O_anr_root = ((Par.mO2root / 6 * 18) * ((Par.hO2root^8) / (Par.hO2root^8 + O2root^8)) ...
   * Alpha + (Par.mO2root / 6 * 18) * ((Par.hO2root^8) / (Par.hO2root^8 + O2root^8)) ...
   * (1-Alpha) * (CH2O*Par.M)) * Par.rootWeight;
% To check relevance of strenght of non-linearity
%CH2O_anr_root = ((Par.mO2root / 6 * 18) * ((Par.hO2root^2) / (Par.hO2root^2 + O2root^2)) ...
 %   * Alpha + (Par.mO2root / 6 * 18) * ((Par.hO2root^2) / (Par.hO2root^2 + O2root^2)) ...
 %  * (1-Alpha) * (CH2O*Par.M)) * Par.rootWeight;

% Total metabolic rate in shoot and root
Mshoot_CH2O = CH2O_aer_shoot + CH2O_anr_shoot;
Mroot_CH2O = CH2O_aer_root + CH2O_anr_root;

% Plant-level total aerobic and anaerobic rate
CH2O_aer = CH2O_aer_shoot + CH2O_aer_root;
CH2O_anr = CH2O_anr_shoot + CH2O_anr_root;

% Plant-level total metabolic rate
CH2O_metabolism = Mshoot_CH2O + Mroot_CH2O;

% Plant-level dynamics of carbohydrate reserves per gram dry weight
dCH2O = (CH2O_photosynthesis - CH2O_metabolism) / (Par.shootWeight + Par.rootWeight);

end