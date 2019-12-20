function [alpha,D]=IonizationMatrix(protein,NaCl_mM,MES_mM)
pH=7;
IncludePhos=1;
CylRadius=5;    % in mn!!
D=30:2:150;          % in mn!!
% NaCl_mM=150;
% MES_mM=100;
aa=0.3;             % amino acid size



[~, PartialCharge, ~, ~, ~, ~, ~]=MCharge(protein,pH,IncludePhos,1,'Tail');
alpha_B=-PartialCharge;
C_H=mM_to_nm((10^-pH)*10^4);
C_s=mM_to_nm(NaCl_mM+MES_mM);
alpha=zeros(length(D),1);
alpha2=zeros(length(D),1);
f_ion=alpha; c_p=zeros(length(D),1); S=zeros(length(D),1); deriv_S=zeros(length(D),1); f_ion=zeros(length(D),1); f_conc=zeros(length(D),1); deriv_c_p=zeros(length(D),1); f_conf=zeros(length(D),1);
x0 = 0.5;

% V_aa=4*pi/3*(aa^3);                   % amino acid volume in nm^3
V_aa=(aa^3);
N_L=length(protein.Tail);
VolTail=N_L*V_aa;
Lambda_Debye=0.304/sqrt((NaCl_mM+MES_mM)*10^-3);
Nu=V_aa+alpha_B*(Lambda_Debye)^3;
% Nu=V_aa;


for i=1:length(D)
    %%
%     
    [S(i),deriv_S(i)]=Surface(D(i),CylRadius);
    c_p(i)=16/45/S(i)*N_L; % density of tail monomers in S
    deriv_c_p(i)=-(16/45)*deriv_S(i)*N_L/(S(i)^2);
    %     fun=@(x)(alpha_B*(1-x)/(x*(1-alpha_B))-x*(1-alpha_B)/(alpha_B*(1-x))-x*c_p(i)/(C_H+C_s));
    %     Zhulina model from 1995. obtains alpha>alpha_B which makes no sense
    lb=0;
    ub=1;
    rng default
    % fbnd=@(x)(1-alpha_B)/(alpha_B*(1-x))-sqrt((x*c_p(i)/(C_H+C_s)^2)+1)+x*c_p(i)/(C_H+C_s);
    fbnd=@(x)(alpha_B*(1-x)/(x*(1-alpha_B))-x*(1-alpha_B)/(alpha_B*(1-x))-x*c_p(i)/(C_H+C_s));
    [alpha(i),~] = lsqnonlin(fbnd,x0,lb,ub);
    

    
end

%
plot(D,alpha)
ylabel('alpha- fraction of aa charged')
xlabel('D [nm] interfilament distance')
% legend(['Charged fraction:',10,protein.Name])



    function out=mM_to_nm(input)
        N_a=6.0221413*10^23; % avogadro's number
        out=input*N_a*10^-27;
    end

    function [S,deriv_S]=Surface(D,CylRadius)
        S=(sqrt(3)*D^2)/4-(pi*CylRadius^2)/2;
        deriv_S=(sqrt(3)*D)/2;
    end


end