Salt_mM=90;                   % mM
a=0;                       % "colloid" radius [nm]
l_b=0.7;                       % Bejurrm [nm]
epsilon=80;                    % dielctric constant
K=8.9875517873681764*10^9;      % Coulomb constant.
KT=4.1164*10^-21;               % K_bT in Joules. 
l_d=0.301/sqrt(Salt_mM/1000);  % Debye length [nm]
e=-1.6021766*10^(-19);         % electron charge [C]  (coulombs)
% phi_point=@(r)(K*e/epsilon./(r*10^-9).*exp(-r/l_d));  

phi_point=@(r)(K*e/epsilon./(r*10^-9)./(1+a/l_d).*exp(-(r-a)/l_d));  % Potential


virial_integrad=@(r)((1-exp(-e*phi_point(r)/KT))*4*pi.*r.^2);  % Virial integrad
Volume=0.5*integral(virial_integrad,a,inf)                      % A_2 virial

2*pi*l_b*l_d^2                                                  % A_2 virial with phi<<KT approximation
