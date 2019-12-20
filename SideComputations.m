% heavily on http://pubs.acs.org/doi/pdf/10.1021/ma051353r

syms x

a=0.38;
b=a;
l_B=0.7;
alpha=0.24;
f=alpha;
salt_M=0.15;
lambda_debye=0.301/sqrt(salt_M);
Radius=0.38;
Diameter=2*Radius;
Length=1.3;
Vol=Length*pi*Radius^2;
u=l_B/a;

% D_e=vpasolve(x == a*(u*alpha^2)^(-1/3)*(log(lambda_debye/x))^(-1/3), x, [0 Inf])
g=@(x)(x/a)^2;
D_e=vpasolve(3*x/(g(x)*a^2)...
    -5/x...
    -(l_B*f^2*(g(x))^2)/(x^2)*log(lambda_debye/x)==0, x, [0 Inf]);

f^2*l_B/a/(a/lambda_debye)
% l_B*alpha^2

% a/lambda_debye

B_el_cyl=(4*pi*l_B*f^2)*lambda_debye^2 % 2nd virial eq. 28,29 from 10.1088/0953-8984/21/42/424112. Limit good for high salts. 

l_k_OSF=(lambda_debye^2*l_B*alpha^2)/(4*a^2)
l_k_DJ=0.32*l_B/a*alpha^2/lambda_debye
% l_k_DJ_WC=-lambda_debye/log(D_e/lambda_debye)
l_k_DJ_WC=double(0.32*l_B*f^2*g(D_e)/D_e*lambda_debye)
lambda_cr=a*exp(1/(u*alpha^2))  %%%%%%%%%%%%%%%%%

% Onsager excluded volume
beta=0.5*pi*Diameter*(Length^2+0.5*(pi+3)*Length*Diameter+0.25*pi*Diameter^2)*Vol;
B_0=beta/2
