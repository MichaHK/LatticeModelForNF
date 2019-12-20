function [Dmin,f_conc_min,f_conf_min,f_i_min]=Lattice(protein,NaCl_mM,MES_mM)
% PROBLEM: f_conc is comparable to f_ion, unlike what Zhulina stated in her
% paper. 
% Why does Zhulina take kuhn lengtn to be 1 aa??
% Switch to Kd instead og alpha_B which was a mistake. Continue with:
% http://www.tau.ac.il/~andelman/reprints/092_jpc_2000_104_11027.pdf and
% perhaps with http://pubs.acs.org/doi/pdf/10.1021/jp990095v

pH=7;
IncludePhos=1;
CylRadius=5;    % in mn!!
D=30:0.2:100;          % in mn!!
% NaCl_mM=150;
% MES_mM=100;
aa=0.3;             % amino acid size
kuhn=1;             % kuhn length in aa's. 



[~, PartialCharge, ~, ~, ~, ~, ~]=MCharge(protein,pH,IncludePhos,1,'Tail');
alpha_B=-PartialCharge;
C_H=mM_to_nm((10^-pH)*10^3);
C_s=mM_to_nm(NaCl_mM+MES_mM);
alpha=zeros(length(D),1);
alpha2=zeros(length(D),1);
f_ion=alpha; c_p=zeros(length(D),1); S=zeros(length(D),1); deriv_S=zeros(length(D),1); f_ion=zeros(length(D),1); f_conc=zeros(length(D),1); deriv_c_p=zeros(length(D),1); f_conf=zeros(length(D),1);
x0 = 0.5;

% V_aa=4*pi/3*(aa^3);                   % amino acid volume in nm^3
V_aa=(aa^3);
N_L=length(protein.Tail);
VolTail=N_L*V_aa;
Lambda_Debye=0.304/sqrt((NaCl_mM+MES_mM)*10^-3);  % consider changing to c_s (instead of C_s)
Nu=V_aa+alpha_B*(Lambda_Debye)^3;
% Nu=V_aa;


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
    f_ion(i)=deriv_S(i)*(45/16)*(sqrt((alpha_B*c_p(i))^2+4*C_s^2)-2*C_s); 
    %%
%     Nu=alpha_B^2/c_p(i);
    f_conc(i)=-Nu*N_L*deriv_c_p(i);
    %%
    f_conf(i)=-(D(i)-2*CylRadius)/N_L*aa^2/kuhn;
    
    
end
[~,I]=min(abs(f_conc+f_ion+f_conf));
Dmin=D(I);
f_conc_min=f_conc(I);
f_conf_min=f_conf(I);
f_i_min=f_ion(I);
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