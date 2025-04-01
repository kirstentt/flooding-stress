function [dO2soil,dO2bulk] = SoilOxygen(Par,Droot,O2soil,O2bulk,PhoStomata)

% Rhizosphere radius [m]
Par.R_rhi   = 2 * Par.Rr;

% Bulk soil radius [m]
Par.R_bulk  = 4 * Par.R_rhi; 

% Rhizosphere depth [m]
Par.Z_rhi   = Par.Zr + Par.Rr;

% Bulk soil depth [m]
Par.Z_bulk  = Par.Z_rhi + Par.R_rhi; 

% Rhizosphere volume [m3]
Par.V_rhi   = pi * Par.R_rhi^2 * Par.Z_rhi - Par.Vr;

% Bulk soil volume [m3]
Par.V_bulk  = pi * Par.R_bulk^2 * Par.Z_bulk - Par.Vr - Par.V_rhi;

% The interface area between rhizosphere and bulk soil
Par.S_rhi   = 2 * pi * Par.R_rhi * Par.Z_rhi + pi * Par.R_rhi^2;

% Soil water level at which oxygen diffusivity at the air-soil interface is half of maximum
Par.Ksoil = Par.Kroot;

Dsoil= Par.Da * (Par.Ksoil^10/(Par.Hwater^10 + Par.Ksoil^10)) ...
        + Par.Dl * (Par.Hwater^10/(Par.Hwater^10 + Par.Ksoil^10));

% Oxygen levels at which rhizosphere and bulk soil respiratory rate is half of max [atm]
hO2soil = 0.01;
hO2bulk = 0.01;

% Superimposing the maximum oxygen consumption rate in rhizosphere and bulk soil
% With only oxygen consumption process the rhizosphere oxygen goes to anoxia within 2 days.
Par.m_rhi  = 0.0032; 
Par.m_bulk = 7.6e-5; 

% Oxygen consumption rates through rhizosphere and bulk soil respiration
rhi_csm = Par.m_rhi * O2soil^2 / (O2soil^2 + hO2soil^2);
bulk_csm = Par.m_bulk * O2bulk^2 / (O2bulk^2 + hO2bulk^2);

% Dynamics of rhizosphere oxygen and bulk soil oxygen
dO2soil = Dsoil* ((Par.O2air - O2soil)/(0.5*Par.Z_rhi) )* (pi * (Par.R_rhi^2 - Par.Rr^2))/ Par.V_rhi ...
    - Droot/Par.V_rhi ...
    - Dsoil * ((O2soil - O2bulk)/(4*Par.Rr)) * Par.S_rhi/Par.V_rhi ...
    - rhi_csm*Par.V_rhi;

dO2bulk = Dsoil * ((Par.O2air - O2bulk)/(0.5*Par.Z_bulk) )* (pi * (Par.R_bulk^2 - Par.R_rhi^2))/ Par.V_bulk ...
    + Dsoil * ((O2soil - O2bulk)/(4*Par.Rr)) * Par.S_rhi / Par.V_bulk ...
    - bulk_csm *Par.V_bulk;


