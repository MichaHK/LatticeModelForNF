function []=main_lattice(handles)
% Hint: get(hObject,'Value') returns toggle state of P_total_checkbox
%        str2double(get(hObject,'String')) returns contents of Salt_edit as a double

D=[19:1:30 31:0.2:90]; 
axes(handles.axes1)
cla
chi=0.5; 
N_a=6.0221413*10^23;

persistent protein;
    
if (isempty(protein))
    protein=struct;
    load('CustomDataForLatticeModelBovine.mat');
end

ProteinList={'INA','NEFL','NEFM','NEFH','VIM','DES','PRPH','GFAP'};
isPhos=[0 1 1 1 0 0 0 0];

ProteinStoichiometry=getProteinStoichiometries(handles);
ProteinRatios=ProteinStoichiometry/sum(ProteinStoichiometry);
NaCl_mM=0; 
MES_mM=str2double(get(handles.Salt_edit,'String'));

P_conf_prefactor=str2double(get(handles.P_conf_edit,'String')); % should be higher for mixes!! It is related to the fraction of tails that are in contact
P_ion_prefactor=str2double(get(handles.P_ion_edit,'String')); 
P_mono_prefactor=str2double(get(handles.P_mono_edit,'String')); 


% Fraction_of_tails_in_contact=0.2;   
pH=6.8;
V_H20=0.03;
aa=str2double(get(handles.aa_edit,'String')); 
% b=get(handles.radio_Polynom);
% % % ModelType=str2double(get(handles.NewValue,'String'));
% 1;
% handles.ButtonSelection
CylRadius=str2double(get(handles.cyl_radius_edit,'String'));

V_aa=4/3*pi*aa^3;
% V_aa=aa^3;
% V_aa=0.3; 
Nu=V_aa;
lambda_Debye=0.304/sqrt((NaCl_mM+MES_mM)*10^-3);
Num_of_tails_in_Vol=32/2;


C_H=mM_to_nm((10^-pH)*10^4);
C_S=mM_to_nm(NaCl_mM+MES_mM);




S=@(D)((sqrt(3)*D.^2)/4-(pi*CylRadius^2)/2);
S_diff=@(D)((sqrt(3)*2*D)/4);
H=@(D)((D-2*CylRadius)/2);
V=@(D)(S(D)*45);
V_diff=@(x)(S_diff(x)*45);


for i=1:length(ProteinList)
    if ProteinRatios(i)>0
        whichIs=find(cellfun(@(c)~isempty(strfind(c, ProteinList{i})), {protein.Gene}));

        if isPhos(i)==1
            PhosState='Phos';
            PhosIndices=protein(whichIs).tail.Sites__Phosphoserine;
        else
            PhosState='NoPhos';
            PhosIndices=[];
        end
        
        sequence=protein(whichIs).tail.Sequence;
        Length(i)=length(sequence);
        
        PartialCharge=protein(whichIs).tail.(PhosState).PartialCharge;
%         [~,~,PartialCharge,~, ~, ~, ~]=...
%             Charge(sequence,PhosIndices,pH);

        alpha_B(i)=abs(PartialCharge);
        c_p{i}=@(D)(ProteinRatios(i)*Num_of_tails_in_Vol*Length(i)./V(D));

        ChargeDensity{i}=alpha_B(i)*c_p{i}(D);
        
        P_conf{i}=-P_conf_prefactor*SetUnits(ProteinRatios(i)*Num_of_tails_in_Vol*sqrt(3)*H(D)/Length(i)/aa^2.*(1/45./D));
        
    else
        c_p{i}=zeros(1,length(D));
        alpha_B(i)=0;
        ChargeDensity{i}=zeros(1,length(D));
        P_conf{i}=zeros(1,length(D));
        Length(i)=0;
    end
end

1;

Num_of_each_tail=ProteinRatios*Num_of_tails_in_Vol;
Num_of_monomers_each_tail=Length.*Num_of_each_tail;

phi_each_tail=Num_of_monomers_each_tail*V_aa./V(D);
phi_each_tail_diff=-V_diff(D).*Num_of_monomers_each_tail*V_aa./(V(D)).^2;


c_p_tot=@(D)(Num_of_tails_in_Vol*dot(Length,ProteinRatios)/sum(ProteinRatios)/45./S(D)); %number density of monomers. 

c_p_tot_diff=@(D)(-Num_of_tails_in_Vol*dot(Length,ProteinRatios)/sum(ProteinRatios)/45./((S(D)).^2).*S_diff(D));

sum_ChargeDensity=sum(vertcat(ChargeDensity{:})); %1 value for each D value.

sum_P_conf=sum(vertcat(P_conf{:})); %Units were already set!
P_ion=P_ion_prefactor*SetUnits(sqrt((sum_ChargeDensity).^2+4*(C_S+C_H)^2)-2*(C_S+C_H));


Num_Water_sites=@(x)(N_a*V(x)*55./10^(24));%%%%%%%%% 
Num_Mono_sites=@(x)(V(x)/V_aa);



% % With monomer sites:
G_FH_diff_V=@(x)(...
    -1./V_diff(x).*Num_Mono_sites(x).*c_p_tot_diff(x).*V_aa...
            .*(-1-log(1-c_p_tot(x).*V_aa)+chi.*(1-2*c_p_tot(x).*V_aa))...
    -(1./V_aa)...
            .*((1-c_p_tot(x).*V_aa).*log(1-c_p_tot(x).*V_aa)+chi*c_p_tot(x).*V_aa.*(1-c_p_tot(x).*V_aa))...
    );

% % % With monomer sites++:
% G_FH_diff_V=@(x)(...
%     -1./V_diff(x).*Num_Mono_sites(x).*c_p_tot_diff(x).*V_aa...
%             .*(-1-log(1-c_p_tot(x).*V_aa)+chi.*(1-2*c_p_tot(x).*V_aa)+%%%%%%%%%%%%%%)...
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


P_mono=SetUnits(G_FH_diff_V(D).*P_mono_prefactor);

% P_mono=P_mono_prefactor*SetUnits(Nu*c_p_tot(D).^2);


% P_mono=P_mono_prefactor*SetUnits(Nu*c_p_tot(D).^2);

P_total=sum_P_conf+P_ion+P_mono;

[Ratio_dist_To_lambda,I]=min(1./(sum_ChargeDensity.^(1/3))/lambda_Debye); % charged monomers do not ES feel each other very much, they are much more distanced the l_d.

% disp(['The minimal ratio between "charged monomer distance" to lambda_{D} is ' num2str(Ratio_dist_To_lambda) ' at ' num2str(D(I))])

axes(handles.axes1)


FileList=get(handles.listbox1,'String');
filename=FileList{get(handles.listbox1,'Value')};


Exp_Pressure=importdata(filename);
Exp_Pressure(:,1)=Exp_Pressure(:,1)/10; 
scatter(Exp_Pressure(:,1),Exp_Pressure(:,2));
set(gca,'yscale','log');
hold on;
PlotModelCurves(P_total,P_ion,P_mono,sum_P_conf,handles,D)



xlim([10,90])
ylim([10,10^7])

ExecuteLine=get(handles.Command_edit,'String');
eval(ExecuteLine);
1;
disp('Volume Fraction 20 and 40')
V_aa.*c_p_tot(D(D==20))
V_aa.*c_p_tot(D(D==40))
V_aa
disp('Mono pressure')
P_mono(D(D==40))
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

function ProteinStoichiometry=getProteinStoichiometries(handles)
    ProteinStoichiometry=[  str2double(get(handles.A_Stoich_edit,'String'))... 
                            str2double(get(handles.L_Stoich_edit,'String')) ...
                            str2double(get(handles.M_Stoich_edit,'String')) ...
                            str2double(get(handles.H_Stoich_edit,'String')) ...
                            0 0 0 0];
    
    % Hint: get(hObject,'Value') returns toggle state of P_total_checkbox
%        str2double(get(hObject,'String')) returns contents of Salt_edit as a double

end

function PlotModelCurves(P_total,P_ion,P_mono,sum_P_conf,handles,D)
    if get(handles.P_total_checkbox,'Value')
        scatter(D,P_total);
    end
    if get(handles.P_conf_checkbox,'Value')
        scatter(D,-sum_P_conf);
    end
    if get(handles.P_ion_checkbox,'Value')
        scatter(D,P_ion);
    end
    if get(handles.P_mono_checkbox,'Value')
        scatter(D,P_mono);
    end
    
end



%%R Old drafts
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