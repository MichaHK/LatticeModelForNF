function handles=main_lattice(handles)
% 1) http://pubs.acs.org/doi/pdf/10.1021/ma1024413 -- The langevin function
% and its Cohen approximation!!!! However, I think the authors lied with
% the A factor they put in... nontheless A~1.6 so no so bad... I think there
% model is very similar to Zhulina: I do not understand why figure 4
% comparing the two is SO different (OK, its zhulina's SCF). 

%% init
% bsxfun(@times...)

D=[19:1:30 31:0.2:90]; 
handles.D=D;
axes(handles.axes1)
cla

pH=6.8;
V_H20=0.03;
N_a=6.0221413*10^23;

persistent protein;
persistent FusedRod;
    
if (isempty(protein))
    protein=struct;
    FusedRod=[];
    load('CustomDataForLatticeModelBovine.mat');
    load('FusedRod.mat');
end

ProteinList=        {'INA','NEFL','NEFM','Nefh','VIM','DES','PRPH','GFAP'}; handles.ProteinList=ProteinList;
SurfaceChargeList=  [0 ,-12.6,-13.9,-6.7,0,0,0,0];

phi_star=0;
handles.aa=str2double(get(handles.aa_edit,'String')); 
% d_a=handles.aa;


handles.CylRadius=str2double(get(handles.cyl_radius_edit,'String'));
chi=str2double(get(handles.chi_flory_edit,'String'));
kuhn=str2double(get(handles.kuhn_length_edit,'String'));
H=@(D)((D-2*handles.CylRadius)/2);

% Hexagonal lattic
S=@(D)((sqrt(3)*D.^2)/4-(pi*handles.CylRadius^2)/2);
S_diff=@(D)((sqrt(3)*2*D)/4);
Num_of_tails_in_Vol=32/2;

% Cylinder lattice
% S=@(x)(pi*(H(x)+handles.CylRadius).^2-(handles.CylRadius).^2);
% S_diff=@(x)(2*pi*H(x));
% Num_of_tails_in_Vol=32; % CHECK!!

% cubic lattice
% S=@(D)(D.^2-handles.CylRadius^2);
% S_diff=@(D)(2*D);
% Num_of_tails_in_Vol=32; % CHECK!!

V=@(D)(S(D)*45);
V_diff=@(x)(S_diff(x)*45);
handles.V_diff=V_diff;
%% Load Zhulina and Experimental files 
s=50; 
FileList=get(handles.listbox1,'String');
Zhu_FileList=get(handles.listbox2,'String');
filename=['ExperimentalData\' FileList{get(handles.listbox1,'Value')}];
Zhu_FileList=get(handles.listbox2,'String');
Zhu_filename=['\ZhuData\' Zhu_FileList{get(handles.listbox2,'Value')}];

Orig_Exp_Pressure=importdata(filename);
Exp_Pressure=Orig_Exp_Pressure;
Exp_Pressure(:,1)=Exp_Pressure(:,1)/10; %%% switching from experimental A to model nm
Zhu_P=importdata(Zhu_filename);
Zhu_P(:,1)=Zhu_P(:,1)/10; 
%% Load Composition and salt concentrations from experimental data 
CompositionFromFile=[];
if get(handles.Auto_checkbox,'Value')
    CompositionFromFile=SetRatiosFromFile(filename,handles);
end
%% Read GUI data for calculations
    isPhos=[get(handles.A_phos,'value') get(handles.L_phos,'value') get(handles.M_phos,'value') get(handles.H_phos,'value') 0 0 0 0];
    ProteinStoichiometry=getProteinStoichiometries(handles);
    ProteinRatios=ProteinStoichiometry/sum(ProteinStoichiometry);
    NaCl_mM=0;
    MES_mM=str2double(get(handles.Salt_edit,'String'));
    lambda_Debye=0.304/sqrt((NaCl_mM+MES_mM)*10^-3);
%     V_aa=4*pi/3*r_a^3; % 0.134 is the value according to partial volume
    V_aa=0.134;
    % V_aa=aa^3;
    % V_elec=2*pi*0.7*lambda_Debye^2;
    % V_aa=V_aa+0.28*V_elec;
    d_a=sqrt(4*V_aa/handles.aa/pi); % reciprocal to the peptide chain, assuming these are cylinders. 
    a=handles.aa;
%     r_a=0.5*handles.aa;
    r_a=0.5*d_a;                    % reciprocal to the peptide chain, assuming these are cylinders.
    
    set(handles.aa_vol_text,'String',['V_aa=' num2str(V_aa,5) ' nm^3']);
    
    hydration=0; % in nm. http://www.ncbi.nlm.nih.gov/pmc/articles/PMC19315/pdf/pq002267.pdf. Probably best to do two layers: inner as hard sphere and outer soft.
    Nu=4*V_aa;
    f_NN=Nu;
    l_B=0.7; % Bjerrum length;
%   f_el=2/3*pi*(aa)^3+2*pi*0.7*lambda_Debye^2./(1+aa/lambda_Debye).*exp(aa/lambda_Debye);        %     model='EM_ex'
    f_el=Nu+2*pi*0.7*lambda_Debye^2;        %     model='EM_ex'
    f_rama=Nu+2*pi*0.7*lambda_Debye^2*4*pi*r_a*log(4*lambda_Debye/r_a*log(lambda_Debye/r_a));        %     model='Ramanathan'
%   f_tom==2/3*pi*(aa)^3+2*pi*0.7*lambda_Debye^2;
    
    f_ee=Nu+2*pi*0.7*lambda_Debye^2;
    f_pe=Nu-2*pi*0.7*lambda_Debye^2;
    f_phph=Nu+2*4*pi*0.7*lambda_Debye^2;
    f_pph=Nu-2*2*pi*0.7*lambda_Debye^2;
    f_eph=Nu+2*2*pi*0.7*lambda_Debye^2;

    
    
    set(handles.debye_text,'String',['Debye=' num2str(lambda_Debye) ' nm']);
    
    C_H=mM_to_nm((10^-pH)*10^4);
    C_S=mM_to_nm(NaCl_mM+MES_mM);
%% Calculate model densities
for i=1:length(ProteinList)
    if ProteinRatios(i)>0
        whichIs=find(cellfun(@(c)~isempty(strfind(c, ProteinList{i})), {protein.Gene}));

        if isPhos(i)==1
            handles.PhosState='Phos';
            PhosIndices=protein(whichIs).tail.Sites__Phosphoserine;
        else
            handles.PhosState='NoPhos';
            PhosIndices=[];
        end
        if strcmp(ProteinList{i},'NEFL')
            NumTruncatedAAs=str2double(get(handles.L_truncate,'string'));
        else
            NumTruncatedAAs=0;
        end
        
        1;
        sequence=protein(whichIs).tail.Sequence;
        sequence=sequence(1:end-NumTruncatedAAs);
        Length(i)=length(sequence);
        
%         PartialCharge=protein(whichIs).tail.(handles.PhosState).PartialCharge;
        AAChargeList=protein(whichIs).tail.(handles.PhosState).AACharges;
        AAChargeList=AAChargeList(1:end-NumTruncatedAAs);
        1;
        
        TotalChargePerTail(i)=sum(AAChargeList);
        PartialCharge(i)=TotalChargePerTail(i)/Length(i);

        alpha_B(i)=abs(PartialCharge(i));
        c_p{i}=(ProteinRatios(i)*Num_of_tails_in_Vol*Length(i)./V(D));
        
        ChargeDensity{i}=alpha_B(i)*c_p{i};
        TailFirstOsmoticTerm{i}=c_p{i}./Length(i);
%         P_conf{i}=-SetUnits(ProteinRatios(i)*Num_of_tails_in_Vol*sqrt(3)*H(D)/Length(i)/d_a^2.*(1/45./D)); %% Fix to work by V_diff!
        
        if ~handles.New_P_Conf_check_Callback
            P_conf{i}=-SetUnits(ProteinRatios(i)*Num_of_tails_in_Vol*...
                3/2*H(D)/Length(i)/a^2./V_diff(D));
        else
            P_conf{i}=-SetUnits(ProteinRatios(i)*Num_of_tails_in_Vol^2*...
                3/2*H(D)/Length(i)/a^2./V_diff(D));
        end
        
        Num_of_Dia_neg{i}=sum(AAChargeList<-1.5);
        c_of_Dia_neg{i}=(Num_of_Dia_neg{i}*ProteinRatios(i)*Num_of_tails_in_Vol./V(D));
        
        Num_of_Mono_neg{i}=sum((AAChargeList<-0.9).*(AAChargeList>-1.1)); 
        c_of_Mono_neg{i}=(Num_of_Mono_neg{i}*ProteinRatios(i)*Num_of_tails_in_Vol./V(D));
        
        Num_of_pos{i}=sum(AAChargeList>0.7);
        c_of_pos{i}=(Num_of_pos{i}*ProteinRatios(i)*Num_of_tails_in_Vol./V(D));
        
    else
        c_p{i}=zeros(1,length(D));
        TailFirstOsmoticTerm{i}=zeros(1,length(D));
        alpha_B(i)=0;
        ChargeDensity{i}=zeros(1,length(D));
        P_conf{i}=zeros(1,length(D));
        Length(i)=0;
        
        c_of_pos{i}=zeros(1,length(D));
        c_of_Mono_neg{i}=zeros(1,length(D));
        c_of_Dia_neg{i}=zeros(1,length(D));
    end
end



% % % % % % c_s=sqrt((sum_ChargeDensity).^2+4*(C_S+C_H)^2)
% % % % % % Fixed_lambda_Debye=

c_p_tot=@(D)(Num_of_tails_in_Vol*dot(Length,ProteinRatios)/sum(ProteinRatios)./V(D)); %number density of monomers. 

c_p_tot_diff=@(D)(-Num_of_tails_in_Vol*dot(Length,ProteinRatios)/sum(ProteinRatios)/45./((S(D)).^2).*S_diff(D));

sum_ChargeDensity=sum(vertcat(ChargeDensity{:})); %1 c_el: value for each D value.
c_el=sum_ChargeDensity;

c_of_Mono_neg_tot = sum(vertcat(c_of_Mono_neg{:}));
c_of_Dia_neg_tot = sum(vertcat(c_of_Dia_neg{:}));
c_of_pos_tot = sum(vertcat(c_of_pos{:}));
c_of_Neut_tot     = c_p_tot(D)-c_of_Mono_neg_tot-c_of_Dia_neg_tot-c_of_pos_tot;
Num_of_mono_in_vol=Num_of_tails_in_Vol*dot(Length,ProteinRatios)/sum(ProteinRatios);

handles.EffAlpha=sum(alpha_B.*ProteinRatios);
% handles.EffAlpha
1;
%% Calculate model selected polymer model
P_conf_prefactor=str2double(get(handles.P_conf_edit,'String')); 
P_ion_prefactor=str2double(get(handles.P_ion_edit,'String')); 
P_mono_prefactor=str2double(get(handles.P_mono_edit,'String')); 

ModelNames=get(handles.Model_listbox,'String');
ModelNum=get(handles.Model_listbox,'Value');
ModelName=ModelNames{ModelNum};
% TailFirstOsmoticTerm_tot=sum(vertcat(TailFirstOsmoticTerm{:}));
TailFirstOsmoticTerm_tot=Num_of_tails_in_Vol./V(D);         % simpler calculation... 
switch ModelName
%     case h.Children(1)
%         model='Some_spring'
    case 'Flory'          % Flory
        model='Flory'
        G_FH_diff_V=@(x)(...
            -c_p_tot(x).^2.*V_aa*chi...
            -(1./V_aa).*log(1-c_p_tot(x).*V_aa)...
            -c_p_tot(x)...
            +Num_of_tails_in_Vol./V(D)...
            );
        P_mono=SetUnits(G_FH_diff_V(D).*P_mono_prefactor);
    case 'Mean Spherical Approximation'         % excluded
%         model='Steric_ex'
        model='Mean Spherical Approximation'
        % following DOI: 10.1007/s10953-007-9131-8. This overshhots because
        % unlike Rubinstein, the first term in the osmotic pressure
        % expansion is not divided by the number of monomers! Therefore,
        % the first term dominates and is huge! 
        eta=c_p_tot(D)*V_aa;
        DisIncludeMonomerEntropy=1;
        P_MSA=c_p_tot(D).*((1+eta+eta.^2-eta.^3)./((1-eta).^3)-DisIncludeMonomerEntropy);
        P_mono=P_mono_prefactor*SetUnits(P_MSA);
    case 'HS theta solvent'    %
        model='HS_theta'
        P_third=10*c_p_tot(D).^3.*V_aa^2;
        P_el2=0;
%         P_el2=(2*pi*0.7*lambda_Debye^2)*c_el.^2;  % already included in
%         the osmotic pressure. 
        P_mono=P_mono_prefactor*SetUnits(P_third+P_el2+TailFirstOsmoticTerm_tot/kuhn);
    case 'Excluded volume with 4*V_aa'   
        model='HS_second'
        P_2=Nu*c_p_tot(D).^2;
%         P_el2=(2*pi*0.7*lambda_Debye^2)*c_el.^2; % already included in
%         the osmotic pressure. 
        P_el2=0;
        P_mono=P_mono_prefactor*SetUnits(P_2+P_el2+TailFirstOsmoticTerm_tot);
    case 'Electrostatic Excluded Volume'             % Electrostatic excluded
        model='Electrostatic Excluded Volume'
        P_el=f_el*c_el.^2+f_NN*(c_p_tot(D)-c_el).^2+2*c_el.*(c_p_tot(D)-c_el)*f_NN;       
        P_mono=P_mono_prefactor*SetUnits(P_el+TailFirstOsmoticTerm_tot);
    case 'Multicomponent Electrostatic'              % multicomponent
        model='Multicomponent Electrostatic'
%         Factor 0.5 only to same ion species, so the factor 2 here is
%         correct. 18, 1426 (1950); doi: 10.1063/1.1747506
%          I also considered 10.1021/ma021456z and  doi: 10.1063/1.472362,
%          but it seems that after derivation of the pressure-  there is
%          very little difference (if any, less pressure). 
% Best explinations: McQuarrie book has it. see page 327
        P_el=f_ee*c_of_Mono_neg_tot.^2+...                
            f_ee*c_of_pos_tot.^2+...
            2*f_pe*c_of_Mono_neg_tot.*c_of_pos_tot+...
            f_phph*c_of_Dia_neg_tot.^2+...
            2*f_pph*c_of_Dia_neg_tot.*c_of_pos_tot+...
            2*f_eph*c_of_Dia_neg_tot.*c_of_Mono_neg_tot+...
            f_NN*c_of_Neut_tot.^2+...
            2*f_NN*c_of_Neut_tot.*(c_of_pos_tot+c_of_Mono_neg_tot+c_of_Dia_neg_tot);
        P_mono=P_mono_prefactor*SetUnits(P_el+TailFirstOsmoticTerm_tot);
    case 'Polyion (Ramanathan)'              % multicomponent
        model='Polyion (Ramanathan)'
        P_el=f_rama*c_el.^2+f_NN*(c_p_tot(D)-c_el).^2+2*c_el.*(c_p_tot(D)-c_el)*f_NN; 
        P_mono=P_mono_prefactor*SetUnits(P_el+TailFirstOsmoticTerm_tot);     
    case 'General Exponent'              
        model='ExponentOnly'
        exponent=str2double(get(handles.Ex_exp_edit,'String'));
        Prefactor_n=str2double(get(handles.A_n_edit,'String'));
        P_mono=P_mono_prefactor*SetUnits(Prefactor_n*c_p_tot(D).^exponent);
        P_conf_prefactor=0; 
        P_ion_prefactor=0; 
    case 'Odijk cylinder'
        model='OdijkCyl'
        %0.1111/j.1749-6632.1949.tb27296.x equation 33 (Onsager 1949 paper)
        l=handles.aa*kuhn;
        A_cyl=0.25*pi*d_a*(l^2+0.5*(pi+3)*l*d_a+0.25*pi^2*d_a^2); 
        P_mono=P_mono_prefactor*SetUnits((c_p_tot(D)/kuhn).^2*A_cyl+TailFirstOsmoticTerm_tot);
    case 'Charged cylinder'
        model='OdijkCylES'
        l=handles.aa*kuhn;
        Vol_Cyl=l*pi*(d_a/2)^2;
        A_cyl=0.25*pi*d_a*(l^2+0.5*(pi+3)*l*d_a+0.25*pi^2*d_a^2);
%         P_cyl2=65*(c_p_tot(D)/kuhn).^3*A_cyl; % This is how the third
%         virial coefficient should have looked like. 
%         P_cyl2=0.422*d_a^2*l^4*(c_p_tot(D)/kuhn).^3; % WRONG...
%         P_cyl2=A_cyl*0.15*d_a*l^2*(c_p_tot(D)/kuhn).^3; % Third Virial Coefficient for the Gas of Long Rods, eq.24 for A=0
        P_cyl1=(c_p_tot(D)/kuhn).^2*A_cyl;
        P_cyl2=2/3*0.64*(c_p_tot(D)/kuhn).^3*l^4*d_a^2;
%         P_cyl2=2/3*2.6841*(c_p_tot(D)/kuhn).^3*l^3*d_a^3;
%         figure
%         plot(D,P_cyl2,D,P_cyl2b)
%         legend('1st','2nd')
%         P_cyl2=0;
%         P_el=A_el*(c_el/kuhn).^2+A_cyl*((c_p_tot(D)-c_el)/kuhn).^2+2*A_cyl*c_el.*(c_p_tot(D)-c_el)*kuhn^(-2);
%         P_el2=(2*pi*0.7*lambda_Debye^2)*c_el.^2; % already included in
%         the osmotic pressure. 
%         P_mono=P_mono_prefactor*SetUnits(P_el2+P_cyl); 
                P_mono=P_mono_prefactor*SetUnits(P_cyl2+P_cyl1+TailFirstOsmoticTerm_tot); 

    case 'Hard-Fused cylinder'
        model='HardFusedCyl'
        % PHYSICAL REVIEW E 91, 042134 (2015)
%         diameter=d_a;
        l=handles.aa*kuhn;
        diameter=2*(3*V_aa/4/pi)^(1/3);
        m_array=2:10;
        L_array=[0.05 0.1 0.25 0.5 0.75 0.9 1 1.25 1.5 1.75 2.0];
        m=kuhn;
        L=1;                % the number of diameters between the centers of two adjacent sphers. L=1 they are tight, L<1 overlap, L=2 a sphere can pass in between two spheres.
        m_ind=(m_array==m);
        L_ind=(L_array==L);
        if (sum(m_ind)~=1 || sum(L_ind)~=1)
            disp('Error in Hard-fused cylinder: L or m value is not available in the simulation data');
        end
        B_2_star=FusedRod(m_ind,L_ind);
%         B_2_star=((11*m-3)+0.25*3.53390*pi*(m-1)^2)/(2*m);  % Only for L=1. D. C. Williamson and G. Jackson, J. Chem. Phys. 108, 10294 (1998).
        if L<0.25
            Cyl_vol=pi*diameter^3*(1+(m-1)*(3*L/2-0.5*L^3))/6; 
        else
            Cyl_vol=1/6*m*pi*diameter^3; 
        end
        B_2=B_2_star*Cyl_vol;
        a_fused=m*L*diameter;
        rho_star=a_fused^(-3);                  % critical number density of rods
        rho_star_mono=m*a_fused^(-3);           % critical number density of mono
        phi_star=rho_star_mono*V_aa;            % critical volume fraction
        P_cyl=B_2*(c_p_tot(D)/m).^2;
%         P_el2=(2*pi*0.7*lambda_Debye^2)*c_el.^2; % already included in
        P_cyl2=B_2*0.15*d_a*l^2*(c_p_tot(D)/kuhn).^3; % Third Virial Coefficient for the Gas of Long Rods, eq.24 for A=0
%         the osmotic pressure. 
        P_el2=0;
        P_mono=P_mono_prefactor*SetUnits(P_el2+P_cyl+P_cyl2);
%         c_p_tot(40)
%         sum_ChargeDensity((D==40))
    case 'Charged Rubinstein Polyelectrolyte'
        model='RubPolyelectro'
% Cite: "Theory of polyelectrolytes in solutions and at surfaces"
        b=a;
        u=l_B/b;                            % the ratio of the Bjerrum length lB to the bond size b
        f_star=sum(alpha_B.*ProteinRatios); % Average partial charge (e/ per monomer)
        c_s=C_S;
%         A_Rub=1/sum(alpha_B.*ProteinRatios);% Average monomer number between charges. 
%         P_el_Rub=((c_p_tot(D)).^2)./(4*A_Rub^2*C_S+A_Rub*c_p_tot(D)); % 1995 paper term for ionic osmotic pressure? If so, remove from here. 
        xi_zero=b*(u*f_star^2)^(-1/7)*(c_p_tot(D)*V_aa).^(-1/2);    % un-numbered equation right after equation 3.17. Exponent -1/6 or -1/7: see supplemental material of our paper. 
        xi=xi_zero.*(1+(2*c_s)./(c_p_tot(D)*f_star)).^(1/4);        % Eq. 3.18
%         xi_thermal=ones(length(xi),1)'*b*kuhn; % see Colby, Rubinstein book eqs. 4.76 and 4.73 for N_B one monomer only. I checked this since I noticed that the Rubinstein P_mono term decreases with incresing salts, even beyond 150mM and in contrast to our experiments. Reason: the thermal correlation length is bigger than the electrostatic. 
%         P_mono_thermal=P_mono_prefactor*SetUnits(xi_thermal.^-3);
        x_ideal=b*(u*f_star^2)^(-1/6)*(c_p_tot(D)*V_aa).^(-1/2).*(1+(2*c_s)./(c_p_tot(D)*f_star)).^(1/4);

        P_mono_real=P_mono_prefactor*SetUnits(xi.^-3);
        P_mono_ideal=P_mono_prefactor*SetUnits(x_ideal.^-3);
        P_mono=P_mono_ideal;
end
handles.ExportModelName=model;
%% Calculate model P_conf and P_ion
sum_P_conf=sum(vertcat(P_conf{:}))/kuhn; %Units were already set!
handles.sum_P_conf=sum_P_conf;
FractionOfRodsOnOutside=0.5;
% TotalSurfaceChargeInVol=0;
TotalSurfaceChargeInVol=-sum(ProteinRatios.*SurfaceChargeList)*Num_of_tails_in_Vol*FractionOfRodsOnOutside;

P_ion=P_ion_prefactor*SetUnits(sqrt((sum_ChargeDensity+TotalSurfaceChargeInVol./V(D)).^2+4*(C_S+C_H)^2)-2*(C_S+C_H));

% Num_Water_sites=@(x)(N_a*V(x)*55./10^(24)); 
% Num_Mono_sites=@(x)(V(x)/V_aa);

if (get(handles.P_minimize_checkbox,'Value')) %&& ~strcmp(ModelName,'General Exponent'))
    P_conf_prefactor=Minimize_P_conf(Exp_Pressure,sum_P_conf,P_ion,P_mono,D,handles);
end
P_conf_value=P_conf_prefactor*sum_P_conf;
P_total=P_conf_value+P_ion+P_mono;

axes(handles.axes1)
%% Plot GUI figure
% this is a total mess. I sould really rewrite it at this point...
x=D; 
y=P_total;
y_exp=Exp_Pressure(:,2);
x_exp=Exp_Pressure(:,1);
handles.IsX_axis_water_volume=(get(handles.popupmenu_X,'Value')==2 || get(handles.popupmenu_X,'Value')==3);
handles.IsX_SlopePlot=(get(handles.popupmenu_X,'Value')==4);
handles.IsLinLinPlot=(get(handles.popupmenu_X,'Value')==5);
handles.IsRubinstein=(get(handles.popupmenu_X,'Value')==6);

if handles.IsX_axis_water_volume
    %     D_to_vol_per_aa=@(y)(45*pi*(y/2-handles.CylRadius).^2)/Num_of_mono_in_vol;
    D_to_vol_per_aa=@(D)(V(D)/Num_of_mono_in_vol);
    x_exp=D_to_vol_per_aa(Exp_Pressure(:,1));
    x=D_to_vol_per_aa(D);
    PlotModelCurves(y,P_ion,P_mono,P_conf_value,handles,x,CompositionFromFile)
    xlim([D_to_vol_per_aa(15),D_to_vol_per_aa(90)])
    ylim([10,10^7])
    xlabel('Volume per monomer [nm^-3]')
    if (get(handles.popupmenu_X,'Value')==3); 
        set(gca,'XScale','log');
    else
        set(gca,'XScale','linear');
    end      
else if handles.IsX_SlopePlot
        [x,y]=TransformSlopePlot(x,y,c_p_tot,V_aa);
        [~,P_ion]=TransformSlopePlot(x,P_ion,c_p_tot,V_aa);
        [~,P_mono]=TransformSlopePlot(x,P_mono,c_p_tot,V_aa);
        [~,P_conf_value]=TransformSlopePlot(x,P_conf_value,c_p_tot,V_aa);
        [x_exp,y_exp]=TransformSlopePlot(x_exp,y_exp,c_p_tot,V_aa);
        PlotModelCurves(y,P_ion,P_mono,P_conf_value,handles,x,CompositionFromFile)
%         xlim auto
%         ylim auto
        Small_c_p_tot=0.03;
        Large_c_p_tot=0.53;
        xhigh=2*(V_aa*Small_c_p_tot)^-1;
        xlow=0.5*(V_aa*Large_c_p_tot)^-1;
        xlim([xlow,xhigh])
        ylim([3*10^-7,1])
%         set(gca,'ytick',10.^(-6:1),'yticklabel',(-6:1));
%         set(gca,'ytick',10.^(-6:1));

        set(gca,'XScale','log');
        xlabel('$\phi^{-1}$','FontSize',20,'Interpreter','Latex')
        ylabel('$\frac{\Pi\nu}{K_{B}T}$','FontName','Times New Roman','FontSize',20,'Interpreter','Latex');
        
%         set(t,'Interpreter','Latex');
 
        1;
    else
        PlotModelCurves(y,P_ion,P_mono,P_conf_value,handles,x,CompositionFromFile)
        handles.x=x; 
        handles.y=y;
        xlim([10,90])
        ylim([10,10^7])
        1;
%         set(gca,'ytick',10.^(1:7),'yticklabel',10.^(1:7));
        set(gca,'XScale','linear');
    end
end
scatter(x_exp,y_exp,'^','k','fill');
handles.x_exp=x_exp;
handles.y_exp=y_exp;

if ~handles.IsLinLinPlot
    set(gca,'yscale','log');
else
    set(gca,'yscale','linear'); % Otherwise the numbers on the axis come out funny. 
    ylim ([0 max(y_exp)*1.1])
        ylim ([0 2*10^5])
end



hold on;
if get(handles.PlotZhu_checkbox,'Value')
        scatter(Zhu_P(:,1),Zhu_P(:,2),s,'s','r');
end
1;
[~,ind_salt_equal_charge_density]=min(abs(C_S-sum_ChargeDensity));
scatter(x(ind_salt_equal_charge_density),y(ind_salt_equal_charge_density),140,'x','k','LineWidth',4)
%% Export model fit
if get(handles.ExportModel,'Value')
    %     ExportModel(y,P_ion,P_mono,sum_P_conf,handles,x,MES_mM,ProteinRatios);
    ExportModel(x,y,x_exp,y_exp,P_ion,P_mono,P_conf_value,MES_mM,ProteinRatios,handles);
end

% if get(handles.ExportATable,'Value')
%     %     ExportModel(y,P_ion,P_mono,sum_P_conf,handles,x,MES_mM,ProteinRatios);
%     ExportATable(D,V_diff,model,EffAlpha,handles);
% end
%% Plot external figure: concentrations
if (0)
    h1=figure(2);
    1;
    clf;
    plot(x,c_p_tot(D),'b');
    hold on;
    plot(x,sum_ChargeDensity,'r');
    C_plot=ones(length(x),1);
    C_plot=C_S*C_plot;
    plot(x,C_plot,'g');
    ylim([0 0.5])
    legend('c_p monomers','\alpha c_p charged monomers','Bulk salt c^+ = c^-')
end
%% Command Prompt notifications
if (0)
ExecuteLine=get(handles.Command_edit,'String');
eval(ExecuteLine);

disp(['Volume Fraction 20nm: ' num2str(V_aa.*c_p_tot(20)) ' and 40nm: ' num2str(V_aa.*c_p_tot(40))])
disp(['Distance Between Monomers 20nm: ' num2str(c_p_tot(20)^-(1/3)) ' and 40nm: ' num2str(c_p_tot(40)^-(1/3))])
disp(['Monomer density at 20nm: ' num2str(c_p_tot(20)) ' and 40nm: ' num2str(c_p_tot(40))])
disp(['Polymer charge density at 20nm: ' num2str(c_el(D==20)) ' and 40nm: ' num2str(c_el(D==40))])

disp(['Bulk ion concentration is: ' num2str(2*C_S) '. The ion density at 20nm is: ' num2str(sqrt((sum_ChargeDensity((D==20))).^2+4*(C_S+C_H)^2))...
    ' and at 40nm: ' num2str(sqrt((sum_ChargeDensity((D==40))).^2+4*(C_S+C_H)^2))]); % 2*C_S since anion and cation

temp=sum(vertcat(c_p{:}));
temp(D(D==20));

set(handles.Tail_Pull_text,'String',['Force_{tail}=?' ' N']);
end
% legend([LegendNames 'Experiment'])
1;
guidata(handles.output, handles);

end
%% Aux functions
function nu_cyl=CylVirialExcludedVolume(Length,Radius)
% L. Onsager. The effects of shape on the interaction of colloidal
% particles. Ann. NY Acad. Sci., 51, 627 (1949). We need \beta_1, which
% according to 20, 21 is the 2nd virial, with \Pi~0.5\beta_1. 
V_cyl=pi*Radius^2*Length;
f=0.5*4.856*V_cyl;
end
function [x,y]=TransformSlopePlot(x,y,c_p_tot,V_aa)
% make a dimensioless pressure and volume fraction. 
    Kb= 1.38065*10^-23;
    T= 298.15;
    x=(c_p_tot(x)*V_aa).^(-1);
    Vol_nm_to_m=@(x)(x*10^(-27));
    y=y*Vol_nm_to_m(V_aa)/Kb/T;
end
function CompositionFromFile=SetRatiosFromFile(filename,handles)
FileTokens=regexp(filename,'(?<CompositionFromFile>[\w-]*)_(?<SaltFromFile>\w*)mM(?<suffix>\w*)','names');
% TruncationToken=regexp(FileTokens.CompositionFromFile,'\w*L(?<TruncatedAA>\d+)(?<Suffix2>\w*)','names');
set(handles.L_phos,'Value',1); set(handles.M_phos,'Value',1); set(handles.H_phos,'Value',1);
a_ratio=0; l_ratio=0; m_ratio=0; h_ratio=0;
Truncated_val=0;
1;
switch FileTokens.CompositionFromFile
    case 'L'
        l_ratio=100;
    case 'L5'
        l_ratio=100;
        Truncated_val=5;
    case 'L11'
        l_ratio=100;
        Truncated_val=11;
    case 'LM'
        l_ratio=70;
        m_ratio=30;
    case 'LH'
        l_ratio=80;
        h_ratio=20;
    case 'AM'
        a_ratio=70;
        m_ratio=30;
    case 'AH'
        a_ratio=80;
        h_ratio=20;
    case 'LMH'
        l_ratio=10/15*100;
        m_ratio=3/15*100;
        h_ratio=2/15*100;
    case 'de-LMH'
        l_ratio=10/15*100;
        m_ratio=3/15*100;
        h_ratio=2/15*100;
        set(handles.L_phos,'Value',0);
        set(handles.M_phos,'Value',0);
        set(handles.H_phos,'Value',0);
    case 'L5M'
        l_ratio=70;
        m_ratio=30;
        Truncated_val=5;
    case 'L11M'
        l_ratio=70;
        m_ratio=30;
        Truncated_val=11;
    case 'L-deH'
        l_ratio=80;
        h_ratio=20;
        set(handles.L_phos,'Value',0);
        set(handles.H_phos,'Value',0);
    case 'L-deM'
        l_ratio=70;
        m_ratio=30;
        set(handles.L_phos,'Value',0);
        set(handles.M_phos,'Value',0);
    otherwise
        disp('ERROR: Unkown composition for filename')
end
set(handles.Salt_edit,'String',FileTokens.SaltFromFile);
set(handles.A_Stoich_edit,'String',num2str(a_ratio));
set(handles.L_Stoich_edit,'String',num2str(l_ratio));
set(handles.M_Stoich_edit,'String',num2str(m_ratio));
set(handles.H_Stoich_edit,'String',num2str(h_ratio));
set(handles.L_truncate,'String',num2str(Truncated_val));
CompositionFromFile=FileTokens.CompositionFromFile;

end

function ExportModel(x,y,x_exp,y_exp,P_ion,P_mono,sum_P_conf,MES_mM,ProteinRatios,handles)
% NOT ALL EXPORT MODELS ARE WORKING WELL. Reason: x is already the data for
% export, and here the function ruins it...
ExportInAngstorm=1;
ProteinComposition=[];
for j=1:numel(ProteinRatios) % make filename with protein ratios
    if ProteinRatios(j)~=0
        if (strcmp(handles.ProteinList{j},'NEFL') && str2double(get(handles.L_truncate,'string'))~=0)
            StringTrucAAs=get(handles.L_truncate,'string');
            TruncationSuffix=['t' StringTrucAAs];
        else
            TruncationSuffix=[];
        end
        ProteinComposition=[ProteinComposition num2str(ProteinRatios(j)*100)  handles.ProteinList{j} TruncationSuffix];
    end
    
end
ext='.out';
aa_string=strrep(num2str(handles.aa), '.', 'p');
Cyl_string=strrep(num2str(handles.CylRadius), '.', 'p');
Prefac_string=strrep(get(handles.P_conf_edit,'String'), '.', 'p');
exponent=strrep(get(handles.Ex_exp_edit,'String'),'.', 'p');

filename=([handles.ExportModelName '_' handles.PhosState '_' ProteinComposition '_' num2str(MES_mM) 'Mono_' 'aa' aa_string '_Cyl' Cyl_string '_Prefac'  Prefac_string]);
currentFolder = pwd;
NewFolder=[currentFolder '\' handles.ExportModelName '\'];
if ~exist(NewFolder,'file')
    mkdir(NewFolder)
end
if ExportInAngstorm
    x=10*x;                  % Model data
    x_exp=10*x_exp;          % Experimental data
end

if handles.IsX_axis_water_volume
    % This has to be fixed. The exported x, x_exp is wrong!
    filename=[NewFolder 'Water_volume_' filename ext];
    ExportMat=[x; y]';
    NaN_vec=nan(numel(x)-numel(x_exp),1);  % also include the experimental data!
    ExportMat=[ExportMat [x_exp; NaN_vec] [y_exp; NaN_vec]];
else if handles.IsX_SlopePlot
        filename=[NewFolder 'NoDimPlot_' filename '_A_n' exponent ext];
        ExportMat=[x; y]';
        NaN_vec=nan(numel(x)-numel(x_exp),1);  % also include the experimental data!
        ExportMat=[ExportMat [x_exp; NaN_vec] [y_exp; NaN_vec]];
    else
        filename=[NewFolder filename ext];
        ExportMat=[x; y; P_ion; P_mono; -sum_P_conf]';
    end
end
%         fid = fopen(filename, 'w');
%         fprintf(fid, [FitsHeader(1:end-1) '\n']); % I should set headers:
%         Model Volume, Model Pressure, "Exp" Volume, exp Pressure
%         fclose(fid)
%
fID=fopen(filename,'w');
csvwrite(filename,ExportMat)
fclose(fID);
disp(['Saved in: ' filename])
end

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

function ProteinStoichiometry=getProteinStoichiometries(handles)
    ProteinStoichiometry=[  str2double(get(handles.A_Stoich_edit,'String'))... 
                            str2double(get(handles.L_Stoich_edit,'String')) ...
                            str2double(get(handles.M_Stoich_edit,'String')) ...
                            str2double(get(handles.H_Stoich_edit,'String')) ...
                            0 0 0 0];
    
    % Hint: get(hObject,'Value') returns toggle state of P_total_checkbox
%        str2double(get(hObject,'String')) returns contents of Salt_edit as a double

end

function PlotModelCurves(P_total,P_ion,P_mono,sum_P_conf,handles,x,CompositionFromFile)
s=10;     
    if get(handles.P_total_checkbox,'Value')
        plot(gca,x,P_total,'r','LineWidth',2);
    end
    if get(handles.P_conf_checkbox,'Value')
        %         scatter(x,-sum_P_conf,'o','m');
        plot(x,-sum_P_conf,'m--','LineWidth',1.25);      
    end
    if get(handles.P_ion_checkbox,'Value')
        %         scatter(x,P_ion,'o','r');
        plot(x,P_ion,'r:','LineWidth',1.25);
    end
    if get(handles.P_mono_checkbox,'Value')
        %         scatter(x,P_mono,'o','g');
        plot(x,P_mono,'g-.','LineWidth',1.25);

    end
    xlabel('Distance [nm]','FontName','Times New Roman')
    ylabel('Pressure[Pa]','FontName','Times New Roman')
    set(gca,'FontSize',18);
    A_n=get(handles.A_n_edit,'String');
    exponent_A=get(handles.Ex_exp_edit,'String');
    A_P_conf_prefactor=get(handles.P_conf_edit,'String');
    Salt_mM=get(handles.Salt_edit,'String');
%     legend(['NF:' CompositionFromFile ' , A_{' exponent_A '} = ' A_n],'Fontsize',20)
    legend(['NF:' CompositionFromFile ' , A=' A_P_conf_prefactor ' , ' Salt_mM 'mM'])%,'FontName','Times New Roman','Fontsize',20)
end

function P_conf_prefactor=Minimize_P_conf(Exp_Pressure,sum_P_conf,P_ion,P_mono,D,handles)
    % change pressure to peg. 
    % take points below 2% (not inc.)
    % minimize. D(PEG)
    %     Exp_Pressure(:,1),Exp_Pressure(:,2)
%     PEG_Data=abs((log10((Exp_Pressure(:,1)-130).*10^-0.57)/2.75).^(1/0.21));
%     FittingInds=(PEG_Data<2);
    UpperPressureLimitForFit=5000; % [Pa]
    PressureValues_forFit=Exp_Pressure((Exp_Pressure(:,2)<UpperPressureLimitForFit),2);
    Exp_dist_forFit=Exp_Pressure((Exp_Pressure(:,2)<UpperPressureLimitForFit),1);
%     P_conf_prefactor_array=[0:0.01:0.15 0.2:0.05:1];
    P_conf_prefactor_array=[0.01:0.01:7];
    deviation=zeros(numel(P_conf_prefactor_array),1);
    for i=1:numel(P_conf_prefactor_array)
        P_total=P_conf_prefactor_array(i)*sum_P_conf+P_ion+P_mono;
        model_dist_forFit=interp1(P_total,D,PressureValues_forFit);               %vq = interp1(x,v,xq)
        deviation(i)=sum(abs(model_dist_forFit-Exp_dist_forFit));
    end
  [~,I]=min(deviation);
  P_conf_prefactor=P_conf_prefactor_array(I);
  set(handles.P_conf_edit,'String',num2str(P_conf_prefactor));  
end
%% Old drafts
% Dual: With H2O sites:
% x=D; 
% G_FH_diff_V=@(x)(...
%     -1./V_diff(x).*Num_Water_sites(x).*c_p_tot_diff(x).*V_aa...
%             .*(-1-log(1-c_p_tot(x).*V_aa)+chi.*(1-2*c_p_tot(x).*V_aa))...
%     -(N_a*55*10^(-24))...
%             .*((1-c_p_tot(x).*V_aa).*log(1-c_p_tot(x).*V_aa)+chi*c_p_tot(x).*V_aa.*(1-c_p_tot(x).*V_aa))...
%     -1./V_diff(x).*Num_Water_sites(x).*c_p_tot_diff(x)*V_H20...
%             .*(1+log(c_p_tot(x).*V_aa))...
%     -(N_a*55*10^(-24))...
%             .*V_H20.*c_p_tot(x).*log(c_p_tot(x).*V_aa)...
%     );
% 
% Dual: With monomer sites:
% x=D; 
% G_FH_diff_V=@(x)(...
%     -1./V_diff(x).*Num_Water_sites(x).*c_p_tot_diff(x).*V_aa...
%             .*(-1-log(1-c_p_tot(x).*V_aa)+chi.*(1-2*c_p_tot(x).*V_aa))...
%     -(N_a*55*10^(-24))...
%             .*((1-c_p_tot(x).*V_aa).*log(1-c_p_tot(x).*V_aa)+chi*c_p_tot(x).*V_aa.*(1-c_p_tot(x).*V_aa))...
%     -1./V_diff(x).*Num_Water_sites(x).*c_p_tot_diff(x)*V_H20...
%             .*(1+log(c_p_tot(x).*V_aa))...
%     -(N_a*55*10^(-24))...
%             .*V_H20.*c_p_tot(x).*log(c_p_tot(x).*V_aa)...
%     );

% % With monomer sites: - DO NOT USE. USE THE EASY CALCULATION
% G_FH_diff_V=@(x)(...
%     -1./V_diff(x).*Num_Mono_sites(x).*c_p_tot_diff(x).*V_aa...
%             .*(-1-log(1-c_p_tot(x).*V_aa)+chi.*(1-2*c_p_tot(x).*V_aa))...
%     -(1./V_aa)...
%             .*((1-c_p_tot(x).*V_aa).*log(1-c_p_tot(x).*V_aa)+chi*c_p_tot(x).*V_aa.*(1-c_p_tot(x).*V_aa))...
%     );

% % With H2O sites:
% G_FH_diff_V=@(x)(...
%     -1./V_diff(x).*Num_Water_sites(x).*c_p_tot_diff(x).*V_aa...
%             .*(-1-log(1-c_p_tot(x).*V_aa)+chi.*(1-2*c_p_tot(x).*V_aa))...
%     -(N_a*55*10^(-24))...
%             .*((1-c_p_tot(x).*V_aa).*log(1-c_p_tot(x).*V_aa)+chi*c_p_tot(x).*V_aa.*(1-c_p_tot(x).*V_aa))...
%     );