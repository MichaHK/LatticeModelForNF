function PressureLattice(protein)
NaCl_mM=70;
MES_mM=70;
pH=7;
IncludePhos=1;
CylRadius=5;    % in nm!!
D=20;          % in nm!!
% NaCl_mM=150;
% MES_mM=100;
aa=0.6;             % amino acid size. 
kuhn=5;             % kuhn length in aa's. 
TwoCoronaDistance=45; % in nm


[Charge, PartialCharge, ~, ~, ~, ~, ~]=MCharge(protein,pH,IncludePhos,1,'Tail');
% alpha_B=-PartialCharge;
alpha_B=-Charge/length(protein.Tail); %>0
C_H=mM_to_nm((10^-pH)*10^3);
C_s=mM_to_nm(NaCl_mM+MES_mM);
alpha=zeros(length(D),1);
% alpha2=zeros(length(D),1);
f_ion=alpha; c_p=zeros(length(D),1); S=zeros(length(D),1); deriv_S=zeros(length(D),1); f_ion=zeros(length(D),1); f_conc=zeros(length(D),1); deriv_c_p=zeros(length(D),1); f_conf=zeros(length(D),1);


% V_aa=4*pi/3*(aa^3);                   % amino acid volume in nm^3
V_aa=3*(aa^3)/4;
N_L=length(protein.Tail);
VolTail=N_L*V_aa;
Lambda_Debye=0.304/sqrt((NaCl_mM+MES_mM)*10^-3);  % consider changing to c_s (instead of C_s)
% Nu=V_aa+alpha_B*(Lambda_Debye)^3;
Nu=V_aa;

  



P_conc=(4096*N_L*Nu)./((sqrt(3)*D^2 - 2*pi*CylRadius^2)^2*TwoCoronaDistance^2);
P_ion=-2*(C_H + C_s) +  2*sqrt((C_H + C_s)^2 + (1024*Charge^2)./((sqrt(3)*D.^2 - 2*pi*CylRadius^2).^2*TwoCoronaDistance^2));
P_conf=-((16*sqrt(3)*(D - 2*CylRadius))./(aa^2*D*N_L*TwoCoronaDistance));
P=SetUnits(P_conc+P_ion+P_conf)
5;

    function out=SetUnits(Pressure) %from nm^-3 to Pa
        Pressure=Pressure*10^27;
        Kb= 1.38065*10^-23; % Boltzmann constant's: J/K  (joules per kelvin) 
        T= 298.15; % Temperature: Kelvins. 
        out=Pressure*Kb*T; 
    end


  function out=mM_to_nm(input)
        N_a=6.0221413*10^23; % avogadro's number
        out=input*N_a*10^-27;
    end

    function [S,deriv_S]=Surface(D,CylRadius)
        S=(sqrt(3)*D^2)/4-(pi*CylRadius^2)/2;
        deriv_S=(sqrt(3)*D)/2;
    end

end