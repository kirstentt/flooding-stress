# flooding-stress
A mechanistic model for simulating plant responses to flooding stress. The model describes oxygen and carbon dynamics of plants subjected to flooding stress, and their dependence on plant rooting depth, aerenchyma levels and ROL barrier induction.

Code accompanies the article:
A minimal mechanistic model for plant stress and acclimation responses under waterlogging

Overview
This project was designed to mechanistically model plant responses to waterlogging, with a particular focus on plant fitness under different acclimation strategies, such as aerenchyma formation and radial oxygen loss (ROL) barrier induction, in plants with varying rooting depths. The code was written in MATLAB, containing 18 files, including seven function files, five main MATLAB files corresponding to different settings of the model that can be run to generate the figures of the main text, and six MATLAB files that can be directly run for sensitivity analysis.

Key features:
The file main.m simulates the temporal dynamics of plant response variables, including root oxygen levels, carbohydrate reserves, stomatal aperture, ROL barrier induction, and others, under specific scenarios of rooting depths, aerenchyma content levels, and the presence or absence of ROL barriers. The dynamics are simulated for a single specific parameter setting, which can be specified by the user.
The four main MATLAB files described next are used to simulate plant dynamics for a two dimensional matrix of parameter conditions, varying two model parameters over a range of values while keeping other model parameters constant. Survival_aeR_root.m simulates plant survival across a range of rooting depths (0.2 m to 0.9 m) and aerenchyma content levels (0 to 0.7), and was run under both the presence and absence of a ROL barrier. When present, the maximum ROL barrier content was set at 0.9. Survival_induction_root.m simulates plant survival across the previously mentioned rooting depth range (0.2 m to 0.9 m) and a range of maximum ROL barrier induction levels (0 to 0.9). Similarly, Survival_induction_aerenchyma.m simulates plant survival across the specified range of aerenchyma content levels (0 to 0.7) and maximum ROL barrier induction levels (0 to 0.9). Finally Survival_InductionTiming.m simulates plant survival across the specified rooting depth (0.2 m to 0.9 m) and aerenchyma content level (0 to 0.7) ranges, and was applied under 3 different settings to compare the effects of early, intermediate, and late timings of ROL barrier induction.
The code also contains 7 function files. ATP.m simulates the dynamics of root ATP concentration. OxygenDiffusion.m simulates the oxygen diffusion rates between air and shoot, and between soil and root. OxygenMetabolism.m simulates the oxygen consumption rate in shoot and root through aerobic respiration. reserves.m simulates the dynamics of plant-level carbohydrate reserves, which are generated through photosynthesis and consumed through aerobic and anaerobic metabolisms. ROL_barrier.m simulates the dynamic induction of ROL barriers. StomataPhotosynthesis.m simulates the stomatal aperture and photosynthetic rate. SoilOxygen.m simulates the dynamics of soil oxygen, which were divided into rhizosphere oxygen and bulk soil oxygen. As soil oxygen dynamics is not the primary focus of this study, it was treated as a boundary condition for plant behavior. 
The 6 sensitivity analysis files examine how sensitive the model is to some of the key parameters.  Sen_canopy_root.m and Sen_canopy_aeT.m perform sensitivity analysis on canopy area, where five canopy area values centered around the reference value used in the main text were selected. Sen_canopy_root.m examines the plant survival time across the specified range of rooting depths, and Sen_canopy_aeT.m examines the plant survival time across the specified range of aerenchyma content levels. Sen_Kp_root.m and Sen_Kp_aeT.m perform the sensitivity analysis on the value of K_p, which refers to the sensitivity of stomatal activities to the change of root ATP level. Again five K_p values centered around the reference value in the main text were selected, and the two files examine the plant survival time across the specified range of rooting depths and aerenchyma content levels respectively. Sen_Dplant_root.m and Sen_Dplant_aeT.m perform the sensitivity analysis on the baseline shoot-root oxygen diffusivity in absence of aerenchyma. Still, five diffusivity values centered around the reference value in the main text were selected. Sen_Dplant_root.m examines the plant survival time across the specified range of rooting depths, and Sen_Dplant_aeT.m examines the plant survival time across the specified range of aerenchyma content levels.


Prerequisites
Before running this code, please make sure you have MATLAB with version R2018a or higher.

Usage
All figures in the manuscript can be generated through the five main files and the six sensitivity analysis files, which can be run directly on MATLAB. Below a detailed mapping from the figure number in the manuscript to the MATLAB file and the figure number that was used to generate it are provided.

Figure mapping 
Article figure number => script name & figure number
Fig 2a => Survival_aet_root Fig 1
Fig 2b => Survival_aet_root Fig 2
Fig 2c=> Survival_induction_root Fig1
Fig 2d=> Survival_induction_aerenchyma Fig1
Fig3a. => main.m Fig15
Fig3b => main.m Fig7
Fig3c => main.m Fig16 
Fig3d=> Survival_aet_root Fig 5
Fig 4a=> main.m Fig8
Fig 4b=> main.m Fig6
Fig 4c=> main.m Fig9
Fig 4d=> Survival_aet_root Fig 3
Fig 5a => main.m Fig3 
Fig 5b=> main.m Fig1 
Fig 5c=> main.m Fig4
Fig 5d=> Survival_aet_root.m Fig 4  

Fig S1a=> main.m Fig13
Fig S1b=> main.m Fig14 
Fig S2 => Survival_aet_root.m Fig1 & Fig2; requiring adjustments in reserves.m; Oxygenmetabolism.m and StomataPhotosynthesis.m (see comments in code)
Fig S3=> Survival_aet_root.m Fig1 & Fig2; requiring adjustments in StomataPhotosynthesis.m (see comments in code)
Fig S4a => main.m Fig5
Fig S4b=> main.m Fig20 
Fig S4c=> main.m Fig17
Fig S4d=>main.m Fig2
Fig S5a=>main.m Fig22
Fig S5b=> main.m Fig10
Fig S5c=>main.m Fig19
Fig S5d=> main.m Fig12
Fig S6a=> Sen_canopy_root Fig1 
Fig S6b => Sen_canopy_aeT Fig1 
Fig S7a => Sen_Kp_root Fig1 
Fig S7b=> Sen_Kp_aeT Fig1 
Fig S8a => Sen_Dplant_root Fig1 
Fig S8b => Sen_Dplant_aeT Fig1 
Fig S9a=> Survival_InductionTiming Fig1
Fig S9b=> Survival_InductionTiming Fig2
Fig S9c=> Survival_InductionTiming Fig3
Fig S10a=> Survival_InductionTiming Fig4
Fig S10b=> Survival_InductionTiming Fig4




