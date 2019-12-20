function []=main_lattice(NaCl_mM,MES_mM,D)
% main_lattice(protein,70,80)
% Interesting: high pressure behaviour not much affected by salt
% concentrations in NFL. 
% Maybe alpha allows a little more M and H, which makes the compression
% harder???
clf; 

ProteinList={'INA','NEFL','NEFM','NEFH','VIM','DES','PRPH','GFAP'};
isPhos=[0 1 1 0 0 0 0 0];
ProteinStoichiometry=[0 3 1 0 0 0 0 0];
ProteinRatios=ProteinStoichiometry/sum(ProteinStoichiometry);
Fraction_of_tails_in_contact=0.2;   % should be higher for mixes!!

pH=7;
aa=0.6; 
CylRadius=5;

V_aa=4/3*aa^3;
Nu=V_aa;
lambda_Debye=0.304/sqrt((NaCl_mM+MES_mM)*10^-3);
Num_of_tails_in_Vol=32/4;

C_H=mM_to_nm((10^-pH)*10^4);
C_S=mM_to_nm(NaCl_mM+MES_mM);

clear protein;
protein=struct;
load('ProteinData.mat');
LegendNames=[];

S=@(D)((sqrt(3)*D.^2)/4-(pi*CylRadius^2)/2);
H=@(D)((D-2*CylRadius)/2);
V=@(D)(S(D)*45);



for i=1:length(ProteinList)
    if ProteinRatios(i)>0
        whichIs=find(cellfun(@(c)~isempty(strfind(c, ProteinList{i})), {protein.Gene}));
%         if i==3
%             1;
%         end
        if isPhos(i)==1
%             PhosState='Phos';
            PhosIndices=protein(whichIs).tail.Sites__Phosphoserine;
%             PhosIndices=PhosIndices(PhosIndices>0);
        else
%             PhosState='NoPhos';
            PhosIndices=[];
        end
        
        
        sequence=protein(whichIs).tail.Sequence;
        Length(i)=length(sequence);
        
        [~,~,PartialCharge,~, ~, ~, ~]=...
            Charge(sequence,PhosIndices,pH);

        alpha_B(i)=abs(PartialCharge);
        c_p{i}=@(D)(ProteinRatios(i)*Num_of_tails_in_Vol*Length(i)./V(D));

        ChargeDensity{i}=alpha_B(i)*c_p{i}(D);
        
        P_conf{i}=-Fraction_of_tails_in_contact*SetUnits(ProteinRatios(i)*Num_of_tails_in_Vol*3.*H(D)/Length(i)/aa^2.*(1/45/sqrt(3)./D));
        
%         P_total(i)=P_conf(i)+P_ion(i)+P_mono;
        
%         scatter(D,P_total)
%         LegendNames=[LegendNames ProteinList{i}];
%         hold on;
    else
        c_p{i}=zeros(1,length(D));
        alpha_B(i)=0;
        ChargeDensity{i}=zeros(1,length(D));
        P_conf{i}=zeros(1,length(D));
        Length(i)=0;
    end
end


% Values_ChargeDensity=cellfun(@(x)x(D),ChargeDensity,'UniformOutput',false);    % cells with values, not functions as in ChargeDensity
% Sum_ChargeDensity=sum(cellfun(@(x)x,Values_ChargeDensity,'UniformOutput',false))
c_p_tot=@(D)(Num_of_tails_in_Vol*dot(Length,ProteinRatios)/sum(ProteinRatios)/45./S(D));
sum_ChargeDensity=sum(vertcat(ChargeDensity{:})); %1 value for each D value.

sum_P_conf=sum(vertcat(P_conf{:})); %Units were already set!
P_ion=SetUnits(sqrt((sum_ChargeDensity).^2+4*(C_S+C_H)^2)-2*(C_S+C_H));

P_mono=SetUnits(Nu*c_p_tot(D).^2);

P_total=sum_P_conf+P_ion+P_mono;

% c_p_tot_compare=sum(cellfun(@(y)y(D),c_p));

% min(sum_ChargeDensity.^-(1/3))          % minimun (in D) average distance between polymer charges
% lambda_Debye
[Ratio_dist_To_lambda,I]=min(1./(sum_ChargeDensity.^(1/3))/lambda_Debye); % charged monomers do not ES feel each other very much, they are much more distanced the l_d.

disp(['The minimal ratio between "charged monomer distance" to lambda_{D} is ' num2str(Ratio_dist_To_lambda) ' at ' num2str(D(I))])

NFL_Pressure=importdata('NFL.csv');
NFL_Pressure(:,1)=NFL_Pressure(:,1)/10; 
scatter(NFL_Pressure(:,1),NFL_Pressure(:,2));
set(gca,'yscale','log');
hold on;
scatter(D,P_total);



% legend([LegendNames 'Experiment'])



    function out=mM_to_nm(input)
        N_a=6.0221413*10^23; % avogadro's number
        out=input*N_a*10^-27;
    end

    function out=SetUnits(Pressure) %from nm^-3 to Pa
        Pressure=Pressure*10^27;
        Kb= 1.38065*10^-23; % Boltzmann constant's: J/K  (joules per kelvin) 
        T= 298.15; % Temperature: Kelvins. 
        out=Pressure*Kb*T; 
    end


end