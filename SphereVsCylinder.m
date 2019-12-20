a=0.34;
Vol=4*pi*a^3/3;
Length=a*(1:10);
B_2=0.5*pi*a*(Length.^2+(pi+3)*Length*a+pi*a^2);
% A_cyl=0.25*pi*Diameter_aa*(l^2+0.25*(pi+3)*l*Diameter_aa+0.25*pi^2*Diameter_aa^2); 


B_2_sphere=4*Vol;
clf
plot(Length/a,B_2./(Length/a),'r')
hold on
plot(Length/a,B_2_sphere.*ones(numel(Length),1))
legend('B_{cyl}','B_{sphere}')
xlabel('l_k / a','Fontsize',14)
ylabel('B_2 per amino acid nm^3','Fontsize',14)
set(gca,'Fontsize',14)
