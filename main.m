%% 02-01-2024 Siluo Chen

close all
clear
clc


%% Simulation scenarios

% Aerenchyma content levels
Aerenchyma = [0,0.66,0.5,0.2];

% Rooting depths [m]
RootDepth  = [0.8,0.6,0.3];

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
Par.Dl = Par.DO2l; %1.29e-10;

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
Par.DairShoot   = Par.DO2a;  %diffusion of oxygen in air

% Max oxygen consumption rate through respiration [muMol min-1 g-1]
mO2_mM     = 3;

% Max shoot-level oxygen consumption rate through respiration [mol s-1]
Par.mO2shoot    = mO2_mM * 1e-6 / 60;

% O2 concentration when metabolism rate is half of max [atm]
Par.hO2shoot    = 5e-2;
Par.hO2root     = 5e-2;

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

dt = 5;
RunTime  = 60*60*24*20 / dt;

% Matrices for data storage
O2shoot_l 			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
O2root_l 			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
O2soil_l 			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
O2bulk_l 			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
Droot_l  			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
Mshoot_l  			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
Mroot_l  			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
Dplant_l  			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
Photo_l  			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
PhoStomata_l 		= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
rootATP_l 				= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
CH2O_l 				= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
CH2O_metabolism_l 	= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
CH2O_aer_root_l 			= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
CH2O_anr_root_l 	= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));
ROLB_l 				= NaN(len*RunTime,length(Aerenchyma)*length(RootDepth)*length(ROLB_on_off));

%% Scenarios of aerenchyma content levels
for m = 1:length(Aerenchyma)

    % Aerenchyma content level (-)
    AeT = Aerenchyma(m);

    %% Scenarios of Rooting depths
    for n = 1:length(RootDepth)

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


        %% Scenarios of ROL barrier absence/presence
        for k = 1:length(ROLB_on_off)

            % Each scenario takes up one vector in the matrix
            MatLayer = MatLayer + 1;

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
            Droot_t 			= NaN(RunTime,1);
            Mshoot_t 			= NaN(RunTime,1);
            Mroot_t 			= NaN(RunTime,1);
            Dplant_t 			= NaN(RunTime,1);
            Photo_t 			= NaN(RunTime,1);
            PhoStomata_t 		= NaN(RunTime,1);
            rootATP_t 			= NaN(RunTime,1);
            CH2O_t 				= NaN(RunTime,1);
            Mroot_CH2O_t 		= NaN(RunTime,1);
            CH2O_metabolism_t 	= NaN(RunTime,1);
            CH2O_aer_root_t 	= NaN(RunTime,1);
            CH2O_anr_root_t 	= NaN(RunTime,1);
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
                    O2shoot_t(t) = O2shoot;
                    O2root_t(t) = O2root;
                    O2soil_t(t) = O2soil;
                    O2bulk_t(t) = O2bulk;
                    Droot_t(t) = Droot;
                    Mshoot_t(t) = Mshoot;
                    Mroot_t(t) = Mroot;
                    Dplant_t(t) = Dplant;
                    Photo_t(t) = AO2;
                    PhoStomata_t(t) = PhoStomata;
                    rootATP_t(t) = rootATP;
                    CH2O_t(t) = CH2O;
                    CH2O_metabolism_t(t) = CH2O_metabolism;
                    CH2O_aer_root_t(t) = CH2O_aer_root;
                    CH2O_anr_root_t(t) = CH2O_anr_root;
                    ROLB_t(t) = ROLB;
                end

                % Output storage
                O2shoot_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = O2shoot_t;
                O2root_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = O2root_t;
                O2soil_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = O2soil_t;
                O2bulk_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = O2bulk_t;
                Droot_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = Droot_t;
                Mshoot_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = Mshoot_t;
                Mroot_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = Mroot_t;
                Dplant_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = Dplant_t;
                Photo_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = Photo_t;
                PhoStomata_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = PhoStomata_t;
                rootATP_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = rootATP_t;
                CH2O_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = CH2O_t;
                CH2O_metabolism_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = CH2O_metabolism_t;
                CH2O_aer_root_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = CH2O_aer_root_t;
                CH2O_anr_root_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = CH2O_anr_root_t;
                ROLB_l((i-1)*(RunTime)+1:i*(RunTime),MatLayer) = ROLB_t;
            end
        end
    end
end

% The ratio between root aerobic to root total metabolic rate
% As an indicator for plotting
CH2O_Mroot_l = CH2O_aer_root_l + CH2O_anr_root_l;
CH2O_ratio_l = CH2O_aer_root_l ./ CH2O_Mroot_l;

% Convert the units of time points from second to hour for plotting
t_s = 1:len*(RunTime);
t_h = t_s / (3600/dt);
t_h_plot = t_h - RunTime / (3600/dt);


figure(1)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,3),'x-k', 'MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,4),'o-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,4)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,9),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,9)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,10),'o-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,10)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,21),'x-g','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,21)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,22),'o-g','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,22)))
xlim([-40 RunTime / (3600/dt)])
ylim([0 1.2])
set(gca,'YAxisLocation','origin')
title('Ratio of root aerobic/total in metabolic rate')
ylabel('Ratio')
xlabel('Time after waterlogging (h)')
set(gca,'FontSize',12)

figure(3)
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,3),'x-k','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,4),'o-k','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,4)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,9),'x-r','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,9)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,10),'o-r','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,10)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,21),'x-g','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,21)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,22),'o-g','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,22)))
hold on
ylabel('Partial pressure in air saturation')
xlabel('Time after waterlogging (h)')
title('Root oxygen')
xlim([-40 RunTime / (3600/dt)])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)
handaxes2 = axes('position',[0.52,0.52,0.36,0.36]);
plot(t_h_plot(RunTime:end),ROLB_l(RunTime:end,4),'o-k','MarkerIndices',1:(3600/dt)*24:length(ROLB_l(RunTime:end,4)))
hold on
plot(t_h_plot(RunTime:end),ROLB_l(RunTime:end,10),'o-r','MarkerIndices',1:(3600/dt)*24:length(ROLB_l(RunTime:end,10)))
hold on
plot(t_h_plot(RunTime:end),ROLB_l(RunTime:end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(ROLB_l(RunTime:end,16)))
hold on
plot(t_h_plot(RunTime:end),ROLB_l(RunTime:end,22),'o-g','MarkerIndices',1:(3600/dt)*24:length(ROLB_l(RunTime:end,22)))
ylim([0 1])
xlim([0 120])
ylabel('Content level')
xlabel('Time after waterlogging (h)')
title('ROL barrier content')
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)


figure(4)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,3),'x-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,4),'o-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,4)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,9),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,9)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,10),'o-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,10)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,21),'x-g','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,21)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,22),'o-g','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,22)))
hold on
line([-40 480],[0.0001 0.0001],'LineStyle','--')
xlim([-40 RunTime / (3600/dt)])
ylim([0 2.5e-3])
ylabel('(mol/g DW)')
xlabel('Time after waterlogging (h)')
title('Carbohydrate reserves')
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(5)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,5)))
hold on
xlim([-40 250])
xlabel('Time after waterlogging (h)')
ylabel('Rate (mol glucose s-1)')
title('Root aerobic metabolism');
set(gca,'FontSize',12)
set(gca,'YAxisLocation','origin')


figure(6)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,13),'x-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,13)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,14)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,17),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,17)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,18)))
hold on
xlim([-40 RunTime / (3600/dt)])
ylim([0 1.2])
xlabel('Time after waterlogging (h)');
title('Ratio of root aerobic/total in metabolic rate')
ylabel('Ratio')
set(gca,'FontSize',12)
set(gca,'YAxisLocation','origin')


figure(7)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_ratio_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_ratio_l(RunTime-24*(3600/dt):end,5)))
hold on
xlim([-40 250])
xlabel('Time after waterlogging (h)')
title('Ratio of root aerobic/total in metabolic rate')
ylabel('Ratio')
set(gca,'FontSize',12)
set(gca,'YAxisLocation','origin')

figure(8)
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,13),'x-k','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,13)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,14)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,17),'x-r','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,17)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,18)))
title('Root oxygen')
xlabel('Time after waterlogging (h)')
ylabel('Partial pressure in air saturation')
xlim([-40 RunTime / (3600/dt)])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)
handaxes2 = axes('position',[0.5,0.5,0.38,0.38]);
plot(t_h_plot(RunTime:end),ROLB_l(RunTime:end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(ROLB_l(RunTime:end,14)))
hold on
plot(t_h_plot(RunTime:end),ROLB_l(RunTime:end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(ROLB_l(RunTime:end,16)))
hold on
plot(t_h_plot(RunTime:end),ROLB_l(RunTime:end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(ROLB_l(RunTime:end,18)))
ylim([0 1])
xlim([0 120])
title('ROL barrier')
ylabel('Content level')
xlabel('Time after waterlogging (h)')
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(9)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,13),'x-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,13)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,14)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,17),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,17)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,18)))
ylabel('(mol/g DW)')
xlabel('Time after waterlogging (h)')
title('Carbohydrate reserves')
line([-40 480],[0.0001 0.0001],'LineStyle','--')
xlim([-40 RunTime / (3600/dt)])
ylim([0 2.5e-3])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(10)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,13),'x-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,13)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,14)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,17),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,17)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_aer_root_l(RunTime-24*(3600/dt):end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_aer_root_l(RunTime-24*(3600/dt):end,18)))
hold on
xlim([-40 RunTime / (3600/dt)])
xlabel('Time after waterlogging (h)');
ylabel('Rate (mol glucose s-1)');
title('Root aerobic metabolism');
set(gca,'FontSize',12)
set(gca,'YAxisLocation','origin')

figure(11)
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,3),'^-k','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,4),'o-k','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,4)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,9),'^-r','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,9)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,10),'o-r','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,10)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,15),'^-b','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,21),'^-g','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,21)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,22),'o-g','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,22)))
title('Stomatal dynamics')
ylabel('Relative stomatal aperture')
xlabel('Time after waterlogging (h)')
xlim([-40 RunTime / (3600/dt)])
ylim([0 1.2])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(12)
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,13),'^-k','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,13)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,14)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,15),'^-b','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,17),'^-r','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,17)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,18)))
title('Stomatal dynamics')
ylabel('Relative stomatal aperture')
xlabel('Time after waterlogging (h)')
xlim([-40 RunTime / (3600/dt)])
ylim([0 1.2])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(13)
plot(t_h_plot(RunTime-24*(3600/dt):end),O2soil_l(RunTime-24*(3600/dt):end,23),'x-r','MarkerIndices',1:(3600/dt)*24:length(O2soil_l(RunTime-24*(3600/dt):end,23)))
hold on
% Plot referencing experimental data
PlotX = [0,12,24,43,53.3,61.7,72.4,77.8,81.8,93.9,104.6,106,120];
PlotY = [18.8,11.7,9.55,6.35,4.62,3.16,2.40,1.63,0.94,0.59,0.035,0,0];
plot(PlotX,PlotY/100,'s-r')
title({'Soil oxygen';'Experimental data and model output'})
ylabel('Partial pressure in air saturation')
xlabel('Time after waterlogging (h)')
xlim([-40 RunTime / (3600/dt)])
legend('A_{0.2}R-, Root 0.3','Experimental data')
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(14)
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,23),'x-r','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,23)))
hold on
PlotX = [0,2,4,6,9,11,13,14,16] * 24;
PlotY = [1,0.74,0.16,0.25,0.088,0.17,0.15,0.088,0.088];
plot(PlotX,PlotY,'s-r')
title({'Stomatal dynamics';'Experimental data and model output'})
ylabel('Relative stomatal aperture')
xlabel('Time after waterlogging (h)')
legend('A_{0.2}R-, Root 0.3','Experimental data')
xlim([-40 RunTime / (3600/dt)])
ylim([0 1.2])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

%% Only shoot/root ratio
figure(15)
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,5)))
title('Root oxygen')
xlabel('Time after waterlogging (h)')
ylabel('Partial pressure in air saturation')
xlim([-40 250])
ylim([0 0.25])
set(gca,'FontSize',12)
set(gca,'YAxisLocation','origin')


%%FigS2d
figure(2)
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2root_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(O2root_l(RunTime-24*(3600/dt):end,5)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2soil_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(O2soil_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2soil_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(O2soil_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),O2soil_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(O2soil_l(RunTime-24*(3600/dt):end,5)))
title('Root & rhizosphere oxygen')
xlabel('Time after waterlogging (h)')
ylabel('Partial pressure in air saturation')
xlim([-40 250])
ylim([0 0.25])
set(gca,'FontSize',12)
set(gca,'YAxisLocation','origin')



figure(16)
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),CH2O_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(CH2O_l(RunTime-24*(3600/dt):end,5)))
title('Carbohydrate reserves')
xlabel('Time after waterlogging (h)')
ylabel('(mol/g DW)')
xlim([-40 250])
line([-40 250],[0.0001 0.0001],'LineStyle','--')
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(17)
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),PhoStomata_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(PhoStomata_l(RunTime-24*(3600/dt):end,5)))
title('Stomatal dynamics')
ylabel('Relative stomatal aperture')
xlabel('Time after waterlogging (h)')
xlim([-40 RunTime / (3600/dt)])
ylim([0 1.2])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(20)
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,1),'x-k','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,3),'x-b','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,5),'x-r','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,5)))
title('Root ATP')
ylabel('Concentration (mol L^{-1})')
xlabel('Time after waterlogging (h)')
xlim([-40 RunTime / (3600/dt)])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(18)
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,3),'^-k','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,4),'o-k','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,4)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,9),'^-r','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,9)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,10),'o-r','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,10)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,15),'^-b','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,21),'^-g','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,21)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,22),'o-g','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,22)))
title('Root ATP')
ylabel('Concentration (mol L^{-1})')
xlabel('Time after waterlogging (h)')
xlim([-40 RunTime / (3600/dt)])
% ylim([0 5e-3])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(19)
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,13),'^-k','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,13)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,14)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,15),'^-b','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,17),'^-r','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,17)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),rootATP_l(RunTime-24*(3600/dt):end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(rootATP_l(RunTime-24*(3600/dt):end,18)))
xlabel('Time after waterlogging (h)')
title('Root ATP')
ylabel('Concentration (mol L^{-1})')
xlim([-40 RunTime / (3600/dt)])
% ylim([0 5e-3])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(21)
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,1),'o-k','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,1)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,3),'o-b','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,3)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,5),'o-r','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,5)))
title('Shoot-root oxygen diffusion')
ylabel('Rate (mol s-1)')
xlabel('Time after waterlogging (h)')
xlim([-40 RunTime / (3600/dt)])
set(gca,'YAxisLocation','origin')
set(gca,'FontSize',12)

figure(22)
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,13),'x-k','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,13)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,14),'o-k','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,14)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,15),'x-b','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,15)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,16),'o-b','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,16)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,17),'x-r','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,17)))
hold on
plot(t_h_plot(RunTime-24*(3600/dt):end),Dplant_l(RunTime-24*(3600/dt):end,18),'o-r','MarkerIndices',1:(3600/dt)*24:length(Dplant_l(RunTime-24*(3600/dt):end,18)))
xlim([-40 RunTime / (3600/dt)])
title('Shoot-root oxygen diffusion')
ylabel('Rate (mol s-1)')
xlabel('Time after waterlogging (h)')
set(gca,'FontSize',12)
set(gca,'YAxisLocation','origin')