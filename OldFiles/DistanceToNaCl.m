NaCl_mM=0:40:420;
MES_mM=80;
protein=alpha;
% figure
Dmin=zeros(1,length(NaCl_mM)); f_conc_min=Dmin; f_i_min=Dmin; f_conf_min=Dmin; ChargeMonomers_mM=Dmin;
for i=1:length(NaCl_mM)
    [Dmin(i),f_conc_min(i),f_conf_min(i),f_i_min(i),ChargeMonomers_mM(i)]=Lattice(protein,NaCl_mM(i),MES_mM);
%     hold on;
end


figure
subplot(3,1,1)
plot(NaCl_mM,f_conc_min./f_i_min)
ylabel('f_{conc min}/f_{ion min}')
xlabel('NaCl [mM]')

subplot(3,1,2)
plot(Dmin,ChargeMonomers_mM)
ylabel('Monomer Charge mM')
xlabel('Distance [nm]')

% plot(NaCl_mM,f_i_min)
% hold on;
% plot(NaCl_mM,f_conf_min,'r')
% ylabel('f_i_min')
% xlabel('NaCl [mM]')
% legend('f_{ion min}','f_{conf min}')

subplot(3,1,3)
plot(NaCl_mM,Dmin)
hold on; 
% plot(NaCl_mM,76*NaCl_mM.^-(1/6),'r')
ylabel('Distance [nm]')
xlabel('NaCl [mM]')