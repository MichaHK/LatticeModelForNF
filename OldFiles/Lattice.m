function [Dmin,f_conc_min,f_conf_min,f_i_min,ChargeMonomers_mM]=Lattice(protein,NaCl_mM,MES_mM)
% Plot also c, cp to see if this is the salt regime (where the protein is
% pH insensitive). Use the term for L=Na: -N*(grafting density)*ln(1-L/Na)
% with a the kuhn length.
% Actually: "Added salt does not affect the brush unless the ionic strength of
% the solution approaches the level of the ionic strength inside the brush (Fig. 6d). In that
% case the prefactor is an inverse cubic root function of the external salt concentration and
% the grafting density (so called salted brush regime)." from "Responsive
% Polymer Brushes" by SERGIY MINKO. It will fix the low salt case (less
% intersting)

% !!The picture I made on the phase transition. The explination is OK with
% the fact it is only observed to be the same for high salts. The
% transition is explained by the contacts. Ask Roy for reversibilty. 
% The Collapse: the collapse is independent of salt and polymer type (=charge) because it is controlled by F_conc and not F_ion 
% !! Fix debye length to workby c_s and not C_s
% !! Have the energy plotted around the minD. Could be handy for
% understaing stability in the later L-H case (where H is outside). 
% model (f_ion, f_conc -> De-3; f_conf-> D

% Switch to Kd instead of alpha_B which was a mistake. Continue with:
% http://www.tau.ac.il/~andelman/reprints/092_jpc_2000_104_11027.pdf and
% perhaps with http://pubs.acs.org/doi/pdf/10.1021/jp990095v
% Consider plotting for the f_total(D) curve at different salts, around their minimum. Maybe the
% minimum is not very deep at some salts- this will 
N_a=6.0221413*10^23; % avogadro's number
pH=6.8;
IncludePhos=0;
CylRadius=5;    % in nm!!
D=15:0.2:100;          % in nm!!
% NaCl_mM=150;
% MES_mM=100;
aa=0.3;             % amino acid size. !!!!!!!!!!!!!!!
kuhn=2;             % kuhn length in aa's. !!!!!!!!!!!!!!!!



[Charge, PartialCharge, ~, ~, ~, ~, ~]=MCharge(protein,pH,IncludePhos,1,'Tail');
% alpha_B=-PartialCharge;
alpha_B=-Charge/length(protein.Tail) %>0
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

for i=1:length(D)
    %%
%     
    [S(i),deriv_S(i)]=Surface(D(i),CylRadius);
    c_p(i)=16/45/S(i)*N_L; % density of tail monomers in S
    deriv_c_p(i)=-(16/45)*deriv_S(i)*N_L/(S(i)^2);
%     %     fun=@(x)(alpha_B*(1-x)/(x*(1-alpha_B))-x*(1-alpha_B)/(alpha_B*(1-x))-x*c_p(i)/(C_H+C_s));
%     %     Zhulina model from 1995. obtains alpha>alpha_B which makes no sense
%     lb=0;
%     ub=1;
%     rng default
%     % fbnd=@(x)(1-alpha_B)/(alpha_B*(1-x))-sqrt((x*c_p(i)/(C_H+C_s)^2)+1)+x*c_p(i)/(C_H+C_s);
%     fbnd=@(x)(alpha_B*(1-x)/(x*(1-alpha_B))-x*(1-alpha_B)/(alpha_B*(1-x))-x*c_p(i)/(C_H+C_s));
%     [alpha(i),~] = lsqnonlin(fbnd,x0,lb,ub);
%     
    %%
    f_ion(i)=(deriv_S(i)*(45/16)*(sqrt((alpha_B*c_p(i))^2+4*(C_s+C_H)^2)-2*(C_s+C_H))+...
        S(i)*(45/16)*alpha_B*c_p(i)*deriv_c_p(i)/sqrt((alpha_B*c_p(i))^2+4*(C_s+C_H)^2));   
    %%
%     Nu=alpha_B^2/c_p(i);
    f_conc(i)=Nu*N_L*deriv_c_p(i);  %<0
    %%
    f_conf(i)=3/2*(D(i)-2*CylRadius)/N_L/aa^2/kuhn; %>0   % 4.25! Should have been 3. 
    
    
end
[~,I]=min(abs(f_conc+f_ion+f_conf));
Dmin=D(I);
f_conc_min=f_conc(I);
f_conf_min=f_conf(I);
f_i_min=f_ion(I);
f_total=f_conc+f_conf+f_i_min;
c_p_min=c_p(I);
ChargeMonomers_mM=alpha_B*c_p_min/N_a/(10^-27);


       

% plot(D,f_total)
% ylim([-0.1 0.1])
%%
% figure
% subplot(1,3,1)
% plot(D,alpha)
% hold on
% ylabel('alpha- fraction of aa charged')
% xlabel('D [nm] interfilament distance')
% legend(['Charged fraction:',10,protein.Name])
% subplot(1,3,2)
% plot(D,S)
% legend('Surface')
% subplot(1,3,3)
% plot(D,c_p)
% legend('c_p')
% %%
% figure
% subplot(1,3,1)
% plot(D,f_ion)
% hold on
% plot(D,f_conc,'r')
% plot(D,f_conf,'c')
% ylabel('f_ion')
% xlabel('D [nm] interfilament distance')
% legend('f_{ion}','f_{conc}','f_{conf}')
% % ylim([-0.1 0.1])
% % legend(['f_{ion}:',10,protein.Name(1:15)])
% subplot(1,3,2)
% plot(D,f_conc+f_ion+f_conf)
% legend(['Total energy. Min at:',num2str(D(I)),'nm'])
% ylim([-0.5 0.5])
% subplot(1,3,3)
% plot(D,c_p)
% legend('c_p')


    function out=mM_to_nm(input)
        N_a=6.0221413*10^23; % avogadro's number
        out=input*N_a*10^-27;
    end

    function [S,deriv_S]=Surface(D,CylRadius)
        S=(sqrt(3)*D^2)/4-(pi*CylRadius^2)/2;
        deriv_S=(sqrt(3)*D)/2;
    end


end