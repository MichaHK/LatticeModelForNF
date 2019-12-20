function ProteinFractionCharge(nfl)
%must plot x to undersatnd changes in pH

CylRadius=5; % in nm!!
protein=nfl;
Seq=protein.Tail;
pH=6.8;
D=12; 
N_L=length(Seq);

pKset = struct('N_term',8.6,'K',10.8,'R',12.5,'H',6.5,'D',3.9,'E',4.1,...
    'C',8.5,'Y',10.1,'C_term',3.6); %maybe place Y=0; K=1; R=1. 
acids='DEY'; AcidFractions=zeros(1,length(acids));
bases='KRH';BasicFractions=zeros(1,length(bases));

NaCl_mM=70;
MES_mM=80;
C_H_M=10^-pH;

C_H=mM_to_nm((10^-pH)*10^3);
C_s=mM_to_nm(NaCl_mM+MES_mM);




S=zeros(length(D),1); c_p=zeros(length(D),1); S=zeros(length(D),1); deriv_S=zeros(length(D),1);

    for i=1:length(acids)
        AcidFractions(i)=AminoAcidFractions(Seq,acids(i));
    end
    for j=1:length(bases)
        BasicFractions(j)=AminoAcidFractions(Seq,bases(j));
    end


i=1;

[S(i),deriv_S(i)]=Surface(D(i),CylRadius);
c_p(i)=16/45/S(i)*N_L; % density of tail monomers in S
global alpha_E alpha_D alpha_Y alpha_R alpha_H alpha_K x tau
syms alpha_E alpha_D alpha_Y alpha_R alpha_H alpha_K x tau
acid_alphas=[alpha_D,alpha_E,alpha_Y];  %cysteine missing here!
basic_alphas=[alpha_K,alpha_R,alpha_H]; 

% tau=c_p(i)*(acid_alphas*AcidFractions'-basic_alphas*BasicFractions');
% tau=sym(c_p(i)*(acid_alphas*AcidFractions'-basic_alphas*BasicFractions'));
% tau=sym(c_p(i)*(acid_alphas*AcidFractions'));
    tau=sym(c_p(i)*(alpha_D*AcidFractions(1)+alpha_E*AcidFractions(2)+alpha_Y*AcidFractions(3)-alpha_H*BasicFractions(3)-alpha_R*BasicFractions(2)-alpha_K*BasicFractions(1)));
%             tau=sym(c_p(i)*(-alpha_K*BasicFractions(1)));
% x=((tau/C_s)+sqrt((tau/C_s)^2+4))/2;
x=sym(((tau/C_s)+sqrt((tau/C_s)^2+4))/2);
% AcidX=@(z)(1-z)/z/C_H_M;
% BasicX=@(z)z/(1-z)/C_H_M;

AcidX=@(z)(1-z)/z;
BasicX=@(z)z/(1-z);

assume(0 <= alpha_D <= 1)
assume(0 <= alpha_E <= 1)
assume(0 <= alpha_Y <= 1)
assume(0 <= alpha_K <= 1)
assume(0 <= alpha_R <= 1)
assume(0 <= alpha_H <= 1)


%         init_guess=[0 1];
%             [ValAlpha_K]=...
%                 vpasolve(x==10^-10.8*BasicX(alpha_K),...
%                 alpha_K,init_guess)

%               init_guess=[0 1];
%             [ValAlpha_K]=...
%                 solve(x==10^-10.8*BasicX(alpha_K),...
%                 alpha_K)
            
% [ValAlpha_D,ValAlpha_E,ValAlpha_Y,ValAlpha_K,ValAlpha_R,ValAlpha_H]=...
%     solve(x==10^-3.9*AcidX(alpha_D),x==10^-4.1*AcidX(alpha_E),x==10^-10.1*AcidX(alpha_Y),...
%     x==10^-10.8*BasicX(alpha_K),x==10^-12.5*BasicX(alpha_R),x==10^-6.5*BasicX(alpha_H),...
%     alpha_D,alpha_E,alpha_Y,alpha_K,alpha_R,alpha_H)

% init_guess=[0 1; 0 1;0 1; 0 1;0 1; 0 1];

% [ValAlpha_D,ValAlpha_E,ValAlpha_Y,ValAlpha_K,ValAlpha_R,ValAlpha_H]=...
%     vpasolve([x==10^-3.9*AcidX(alpha_D),x==10^-4.1*AcidX(alpha_E),x==10^-10.1*AcidX(alpha_Y),...
%     x==10^-10.8*BasicX(alpha_K),x==10^-12.5*BasicX(alpha_R),x==10^-6.5*BasicX(alpha_H)],...
%     [alpha_D,alpha_E,alpha_Y,alpha_K,alpha_R,alpha_H],init_guess)
old = digits;
digits(32)
s=0.000001;
l=0.99999;
        init_guess=[s l;s l;s l;s l; s l; s l];
    [ValAlpha_D,ValAlpha_E,ValAlpha_H,ValAlpha_Y,ValAlpha_R,ValAlpha_K]=...
        vpasolve([x==10^(pH-3.9)*AcidX(alpha_D),x==10^(pH-4.1)*AcidX(alpha_E),x==10^(pH-6.5)*BasicX(alpha_H),alpha_Y==0.000001,alpha_R==0.99999,alpha_K==0.99999],...
        [alpha_D,alpha_E,alpha_H,alpha_Y,alpha_R,alpha_K],init_guess)

    TotalCharge=-[ValAlpha_D,ValAlpha_E,ValAlpha_Y]*(AcidFractions*N_L)'+[ValAlpha_K,ValAlpha_R,ValAlpha_H]*(BasicFractions*N_L)'
    
%       init_guess=[0 1; 0 1;0 1];
%     [ValAlpha_D,ValAlpha_E,ValAlpha_H,ValAlpha_Y]=...
%         solve(x==10^(pH-3.9)*AcidX(alpha_D),x==10^(pH-4.1)*AcidX(alpha_E),x==10^(pH-6.5)*BasicX(alpha_H),x==10^(pH-10.1)*AcidX(alpha_Y),...
%         alpha_D,alpha_E,alpha_H,alpha_Y)


    function out=mM_to_nm(input)
        N_a=6.0221413*10^23; % avogadro's number
        out=input*N_a*10^-27;
    end

    function [S,deriv_S]=Surface(D,CylRadius)
        S=(sqrt(3)*D^2)/4-(pi*CylRadius^2)/2;
        deriv_S=(sqrt(3)*D)/2;
    end

    function AAfracion=AminoAcidFractions(Seq,letter) %validated against uniprot for nfl tail
        if ~ischar(letter)
            disp('error: AAfraction function got a letter which is not an AA')       
        else
            Seq=upper(Seq);
            count=length(strfind(Seq,letter));
            AAfracion=count/length(Seq);
        end
    end
end