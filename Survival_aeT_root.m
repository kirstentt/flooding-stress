%% 02-01-2024 Siluo Chen (Jerry)

close all
clear
clc

%% Simulation scenarios

% Aerenchyma content levels
Aerenchyma = 0:0.1:0.7;

% Rooting depths [m]
RootDepth  = 0.2:0.1:0.9;

% ROL barrier absence/presence
ROLB_on_off = [0,1];

%% Plant architectural parameters

% Shoot height [m]
Par.Zp = 0.5;

% Root radius [m]
Par.Rr = 0.01;

% shoot radius [m] 
Par.Rp = 0.01;

% Canopy radius [m]
Par.Rc = 0.2;

% Canopy thickness [m]
Par.Zc = 0.01;

% Canopy area [m2]
Par.Sc = pi * Par.Rc^2;

% Shoot volume [m3]
Par.Vp = pi * Par.Rp^2 * Par.Zp + pi * Par.Rc^2 * Par.Zc; 

% Shoot surface area
Par.Sp = (Par.Zp - Par.Zc) * (2 * pi * Par.Rp) ... %surface stem minus part covered by canopy
    + 2*Par.Sc - pi * Par.Rp^2 ...  %surface canopy + once minus cross area stem
    + 2 * pi * Par.Rc * Par.Zc; %outer edge canopy

	
% Shoot-root cross-sectional area
Par.S = pi * Par.Rr^2;

% Shoot dry weight [g]
Par.shootWeight = 10;


%% Physical parameters

% Ideal gas constant [Pa m3 mol-1 K-1]
Par.R = 8.314;

% Soil porosity [m3 m-3] 
Par.f = 0.5;

% Gas-filled soil porosity [m3 m-3]
Par.Theta = 0.3;    

% Temperature [degree C]
Par.T = 20;

% Molar volume at STP [m3 mol-1]    
Par.Vm = 22.4 / 1000 / 273.15 * (273.15 + Par.T); 

% Molar mass of glucose [g mol-1]
Par.M = 180; 

% Conversion constant from partial pressure to concentration based on Ideal gas law
Par.mol_P = Par.R * (Par.T + 273);

% Partial pressure of oxygen in air [atm]
Par.O2air   = 0.2;

% O2 diffusion coefficient in air [m2 s-1]
Par.DO2a = 2.02e-5; 

% O2 diffusion coefficient in water [m2 s-1]                  
Par.DO2l = 2.1e-9;

% O2 diffusion coefficient in soil air [m2 s-1]
% Note how the 2.5*3 effectively rescales theta from 0.3 to 0.65 and f from
% 0.5 to 0.7, thereby simulating more sandy soil, by removing the2.5*3 we
% move to loamy soil
Par.Da = 2.5*3*Par.DO2a * Par.Theta^(10/3) / Par.f^2; 

% O2 diffusion coefficient in soil water [m2 s-1]
Par.Dl = Par.DO2l;


%% Plant physiological parameters

% Default carbohydrate reserve level in plants [mol g-1 DW]
Par.CH2O        = 0.001;

% Carbohydrate reserve level at which photosynthesis is half of max
Par.KCH2O       = Par.CH2O / 5;

% Carbohydrate reserve level at which respiration is half of max
Par.hCH2O       = Par.CH2O / 2; 

% Root ATP level at which stomatal conductance is half of max
Par.Kp          = 3.2e-3;

% Minimum proportion of stomatal aperture under waterlogging
Par.BetaStomata = 0.1;

% Air-shoot oxygen diffusivity [mol m-2 s-1]
Par.DairShoot   = Par.DO2a; %1e-3;

% Max oxygen consumption rate through respiration [muMol min-1 g-1] 
mO2_mM     = 3;

% Max shoot-level oxygen consumption rate through respiration [mol s-1]
Par.mO2shoot    = mO2_mM * 1e-6 / 60;

% O2 concentration when metabolism rate is half of max [atm]
%Default
Par.hO2shoot    = 5e-2;
Par.hO2root     = 5e-2;
%to compensate in case of power 2 for O2 dependence of metabolism
%Par.hO2shoot    = 0.5*5e-2;
%Par.hO2root     = 0.5*5e-2;

% Shoot-root oxygen diffusivity [s-1]
Par.DO2plant    = 2.15*1e-8;

% Complete induction of aerenchyma increases shoot-root oxygen diffusion up to 2000 folds.
Par.Dae = 2000;

% Rhizosphere oxygen threshold level for ROL barrier induction [atm]
Par.ROLB_induction = 0.1;


%% Simulation settings
MatLayer = 0;

% The amount of water stress events
len      = 2;

dt = 5; %1;
RunTime  = 60*60*24*20 / dt;

% Matrices for data storage
CH2O_l 				= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));


%% Scenarios of ROL barrier absence/presence
for k = 1:length(ROLB_on_off)

    %% Scenarios of aerenchyma content levels 
	for m = 1:length(Aerenchyma)
		
		% Aerenchyma content level (-)
		AeT = Aerenchyma(m);
		
		%% Scenarios of Rooting depths
		for n = 1:length(RootDepth)
		
			% Each scenario takes up one vector in the matrix
			MatLayer = MatLayer + 1;

			% Rooting depth
			Par.Zr    = RootDepth(n);
			
			% Root dry weight is proportional to the rooting depth
			Par.rootWeight  = 10 * Par.Zr;
			
			% Max photosynthetic rate is rooting depth dependent to keep steady-state carbohydrate reserve constant
			% [mol m-2 s-1]
			Par.AO2         = 30e-6 * Par.Sc * (Par.shootWeight + Par.rootWeight) ...
				/ (Par.shootWeight + 4);
				
			% Root volumn [m3]
			Par.Vr          = Par.Zr * Par.S;
			
			% Root surface area [m2]
			Par.Sr          = Par.Zr * (2 * pi * Par.Rr)+ pi*Par.Rr^2;
			
			% Michaelis-Menten parameter for soil-root diffusivity
			Par.Kroot       = 0.5 * Par.Zr;
			
			% Max root-level O2 consumption rate [mol s-1];
			Par.mO2root     = mO2_mM * 1e-6 / 60;
			
			%% Water stress events
				
			% Soil water level at non-stressed condition
			Normal = 0.5 * Par.Zr;
				
			% Soil water level under waterlogging
			Waterlogging = Par.Zr;
			Hwater       = [Normal,Waterlogging];
            
            % Initial values renewed for each scenario
            AO2_ini = Par.AO2;
            O2shoot_ini = 0.2;
            O2root_ini = 0.18;
            O2soil_ini = 0.18;
            O2bulk_ini = 0.18;
            CH2O_ini = Par.CH2O;
            ROLB_ini = 0;
            
			% Vectors for data storage for each scenario
			% Renewed for each scenario
            O2shoot_t			= NaN(RunTime,1);
            O2root_t 			= NaN(RunTime,1);
            O2soil_t 			= NaN(RunTime,1);
            O2bulk_t 			= NaN(RunTime,1);
            Photo_t 			= NaN(RunTime,1);
            CH2O_t 				= NaN(RunTime,1);
            ROLB_t 				= NaN(RunTime,1);
            
			%% Simulation through different non-stress/waterlogging
            for i = 1:len
                
				% Soil water level
                Par.Hwater      = Hwater(i);
                
				% When soil water level switch from non-stress to waterlogging
				% Initial conditions becomes the steady state at non-stress condition
                if Par.Hwater == Waterlogging
                    O2soil_ini = O2soil_t(end);
                    O2bulk_ini = O2bulk_t(end);
                    O2shoot_ini = O2shoot_t(end);
                    O2root_ini = O2root_t(end);
                    AO2_ini  = Photo_t(end);
                    CH2O_ini = CH2O_t(end);
                    ROLB_ini = ROLB_t(end);
                end
                
				% After each simulation the variables start again from the initial condition
                O2soil = O2soil_ini;
                O2bulk = O2bulk_ini;
                O2shoot = O2shoot_ini;
                O2root = O2root_ini;
                AO2    = AO2_ini;
                CH2O   = CH2O_ini;
                ROLB = ROLB_ini;
				
				%% Simulation
                for t = 1:RunTime
				
					% ROL barrier induction
                    if ROLB_on_off(k) == 1
                        dROLB = ROL_barrier(Par,O2soil);
                    else
                        dROLB = 0;
                    end
					
					% Maximum ROL barrier content is set as 0.9
					ROLB_max = 0.9;
                    ROLB = min(ROLB + dROLB*dt,ROLB_max);
                    
					% Root ATP
                    rootATP = ATP(Par,AO2,CH2O,O2shoot,O2root);
					
					% Relative stomatal aperture and photosynthetic rate
                    [PhoStomata,AO2] = StomataPhotosynthesis(Par,AO2,CH2O,O2shoot,O2root);
			
					% Effective diffusion rate between the ambient environment and the shoot and root
                    [Dshoot,Droot] = OxygenDiffusion(Par,PhoStomata,O2shoot,O2root,O2soil,ROLB);
                    
					% Oxygen diffusion rate between shoot and root
                    Dplant = Par.DO2plant * (O2shoot - O2root) * (1 + Par.Dae * AeT) ...
                        * Par.S / ((Par.Zp+Par.Zr)/2);
                    
					% Oxygen consumption rate in shoot and root
                    [Mshoot,Mroot] = OxygenMetabolism(Par,CH2O,O2shoot,O2root);
					
					% Dynamics in carbohydrate reserves
                    [dCH2O,CH2O_metabolism,CH2O_aer_root,CH2O_anr_root] = reserves(Par,AO2,CH2O,O2shoot,O2root);
                    
					% Dynamics in rhizosphere and bulk soil oxygen levels
                    [dO2soil,dO2bulk] = SoilOxygen(Par,Droot,O2soil,O2bulk,PhoStomata);
					
					% Dynamics in shoot and root oxygen levels
					% The photosynthetic and respiratory rates are in the unit of [mol s-1]
					% They are thus converted to partial pressure [atm]
                    dO2shoot = ((AO2 - Mshoot) * Par.mol_P * 1e-5 + Dshoot - Dplant) / Par.Vp ;
                    dO2root  = (Droot+ Dplant - Mroot * Par.mol_P * 1e-5) / Par.Vr ;
                    
                    % Integration of carbohydrate reserves
                    CH2O    = CH2O + dCH2O * dt;
					
					% Plant survival when carbohydrate reserves stay beyond 10% of the default value
                    if  CH2O > (Par.CH2O / 10)
					
						% Integration of rhizosphere oxygen level
						O2soil=O2soil+dO2soil*dt;

						% Integration of bulk soil oxygen level
						O2bulk=O2bulk+dO2bulk*dt;

						% Integration of shoot and root oxygen levels
                        O2shoot = O2shoot + dO2shoot*dt;
                        O2root  = O2root  + dO2root*dt;
                    
					% Plant death
					else
                        CH2O = NaN;
                        O2shoot = NaN;
                        O2root = NaN;
                        PhoStomata = NaN;
                        AO2 = NaN;
                        ROLB = NaN;
                        rootATP = NaN;
                        O2soil = NaN;
                    end
                    
					% Output storage
                    CH2O_t(t)    = CH2O;
					O2shoot_t(t) = O2shoot;
					O2root_t(t)  = O2root;
					O2soil_t(t)  = O2soil;
					O2bulk_t(t)  = O2bulk;
					Photo_t(t) 	 = AO2;
					ROLB_t(t)    = ROLB;
                end
                
				% Output storage
                CH2O_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = CH2O_t;
            end
            survivals = (find(isnan(CH2O_l(:,MatLayer))) - RunTime) / (3600/dt);
            if isempty(survivals)
                survive(MatLayer) = RunTime / (3600/dt);
            else
                survive(MatLayer) = survivals(1);
            end
        end
    end
    
	
	%% Plotting
    SurvTime = reshape(survive(MatLayer-MatLayer/k+1:MatLayer),[length(RootDepth),length(Aerenchyma)]);
    figure(k)
    imagesc(RootDepth,Aerenchyma,transpose(SurvTime))
    xlabel('Rooting depth (m)');
    ylabel('Aerenchyma content');
    h = colorbar;
    title(h,'Survival time (h)')
    % axis off
    set(gca,'YDir','normal')
    set(gca,'FontSize',12)
end


%%Plotting
figure(3)
plot(RootDepth,survive(41:48),'x-b','MarkerSize',12)
hold on
plot(RootDepth,survive(105:112),'o-b')
xlabel('Rooting depth (m)');
ylabel('Survival time (h)');
set(gca,'YAxisLocation','origin');
set(gca,'FontSize',12);

figure(4)
plot(Aerenchyma,[survive(5),survive(13),survive(21),survive(29),survive(37), ...
    survive(45),survive(53),survive(61)],'x-b','MarkerSize',12)
hold on
plot(Aerenchyma,[survive(69),survive(77),survive(85),survive(93),survive(101), ...
    survive(109),survive(117),survive(125)],'o-b')
xlabel('Aerenchyma content');
ylabel('Survival time (h)');
set(gca,'YAxisLocation','origin');
set(gca,'FontSize',12);

figure(5)
plot(RootDepth,survive(1:8),'x-b','MarkerSize',10)
xlabel('Rooting depth (m)');
ylabel('Survival time (h)');
set(gca,'YAxisLocation','origin');
set(gca,'FontSize',12);
set(gca,'YAxisLocation','origin');