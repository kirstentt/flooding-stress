function [Dshoot,Droot] = OxygenDiffusion(Par,PhoStomata,O2shoot,O2root,O2soil,ROLB)

% Effective plant-environment diffusivity to oxygen in shoot and root
Par.Ksoil = Par.Kroot;
Dsoil= Par.Da * (Par.Ksoil^10/(Par.Hwater^10 + Par.Ksoil^10)) ...
        + Par.Dl * (Par.Hwater^10/(Par.Hwater^10 + Par.Ksoil^10));

% The effective air-shoot and soil-root diffusion rates
%This one is incorrect in terms of dimensions but could be seen (given that
%Zc=0.01) as a 200 times smaller effective diffusion
%We use it to speed up simulations, as the normal effective diffusion
%requires a 5 times smaller timestep in order to avoid numerical instabilities
%Results were checked between the original and sped up diffusion and do not differ
%This can be understood from the fact that shoot oxygen hardly changes and 
%soil and intraplant diffusion are limiting anyways,irrespective of this
%200 fold slowdown or not.
Dshoot    = Par.DairShoot * PhoStomata * (Par.O2air - O2shoot) * Par.Sp;
%The actual dimensionally correct equation (requires setting dt=1 instead of 5):
%Dshoot    = (Par.DairShoot/(Par.Zc/2)) * PhoStomata * (Par.O2air - O2shoot) * Par.Sp;
Droot     = Dsoil * (1 - ROLB) * ((O2soil - O2root)/(1.5*Par.Rr)) * Par.Sr;

end