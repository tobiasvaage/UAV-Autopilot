function [fa, ma] = aero(Va,AOA,beta,angrate,delta_e,delta_a,flag)
rho = 1.225;
S = 0.8;
b = 2.1;
c = 0.357;
AR = (b^2)/S;
e = 0.9;    % Oswald efficiency factor
p = angrate(1);
q = angrate(2);
r = angrate(3);

% Aerodynamic coefficients (WT data from Gryte (2018))
% C_L_AOA = (pi*AR)/(1+sqrt(1+(AR/2)^2)); % Small aircraft approx.
CL_0 = 0.0867;  CL_AOA = 4.02; CL_elev = 0.278;   CL_q = 3.87;
CD_0 = 0.0197;  CD_AOA1 = 0.0791;   CD_AOA2 = 1.06; CD_elev = 0.0633;
CD_beta2 = 0.148;   CD_beta1 = -0.00584; CD_q = 0;   Cm_0 = 0.0302;
Cm_AOA = -0.126;    Cm_elev = -0.206;   Cm_q = -1.3;    CY_0 = 0.00316;
CY_beta = -0.224;   CY_ail = 0.0433;    CY_p = -0.137;  CY_r = 0.0839;
Cl_0 = 0.00413; Cl_beta = -0.0849;  Cl_ail = 0.12;  Cl_p = -0.404;
Cl_r = 0.0555;  Cn_0 = -0.000471;   Cn_beta = 0.0283;   Cn_ail = -0.00339;
Cn_p = 0.00437; Cn_r = -0.012;

if flag == 1      % Nonlinear
    CD = CD_0 + CD_AOA1*AOA + CD_AOA2*AOA^2 + CD_elev*delta_e^2 + ...
        CD_q*(c*q)/(2*Va) + CD_beta2*beta^2;
    CY = CY_0 + CY_beta*beta + CY_p*(b*p)/(2*Va) + CY_r*(b*r)/(2*Va) ...
        + CY_ail*delta_a;
    CL = CL_0 + CL_AOA*AOA + CL_q*(c*q)/(2*Va) + CL_elev*delta_e;
    Cl = Cl_0 + Cl_beta*beta + Cl_p*(b*p)/(2*Va) + Cl_r*(b*r)/(2*Va) + ...
        Cl_ail*delta_a;
    Cm = Cm_0 + Cm_AOA*AOA + Cm_q*(b*q)/(2*Va) + Cm_elev*delta_e;
    Cn = Cn_0 + Cn_beta*beta + Cn_p*(b*p)/(2*Va) + Cn_r*(b*r)/(2*Va) + ...
        Cn_ail*delta_a;
    
    D = (1/2)*rho*Va^2*S*CD;
    Y = (1/2)*rho*Va^2*S*CY;
    L = (1/2)*rho*Va^2*S*CL;
    l = (1/2)*rho*Va^2*S*b*Cl;
    m = (1/2)*rho*Va^2*S*c*Cm;
    n = (1/2)*rho*Va^2*S*b*Cn;
    
    % wind to body rotation matrix
    R = [cos(beta)*cos(AOA), -sin(beta)*cos(AOA), -sin(AOA);
         sin(beta)         , cos(beta)          , 0;
         cos(beta)*sin(AOA), -sin(beta)*sin(AOA), cos(AOA)];
     
    fa = R*[-D, Y, -L]';
    ma = [l, m, n]';
else
    C_L = C_L0 + C_L_AOA * AOA;
    C_D = C_D0 + C_D_AOA * AOA;


end