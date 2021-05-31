function [A_pr] = A_pr(model_param, X_HAT, ts)
% RE-LINEARIZATIONG AT CURRENT TIME STEP 
% computes A-prime matrix at current time step given previous estimates

% INPUTS: 
%	model_param {struct}: A struct containing model parameters, 
%   X_HAT {vector}: A Posteriori State Estimation [K-1|K-1]
%   ts: sampling interval 

% OUTPUTS:
%   A_pr[K]: re-linearized DT matrix A given previous a posteriori
%   estimate

% ----------------------------------------------------------------------

% DECOMPOSE into components
C = model_param.C;         % Pac. Tire Hyperparam.
B = model_param.B;         % Pac. Tire Hyperparam.
E = model_param.E;         % Pac. Tire Hyperparam.
r_e = model_param.r_e;     % Effective Tire Radius [m]; 0.37338;
J = model_param.J;         % Wheel Rotational Inertia [kg-m^2]
m = model_param.m;         % Vehicle Mass [kg]
Fz = model_param.Fz;       % Tire Normal Force [N]

U_HAT = X_HAT(1);
W_HAT = X_HAT(2); 
MU_HAT = X_HAT(3);

S_HAT = ( (r_e * W_HAT) / U_HAT ) - 1;
T_HAT = B*S_HAT - B*E*atan(B*S_HAT)*S_HAT;

% Build the matrix
alp = 4 / m; 
beta = -r_e / J;
gam = C * Fz * MU_HAT;

del1 = cos(C * atan(T_HAT) );
del2 = B*r_e*( E*W_HAT*atan(B*S_HAT) - W_HAT ) + (B^2)*E*r_e*W_HAT*S_HAT / (1 + (B*S_HAT)^2);
del3 = (1 + T_HAT^2) * U_HAT^2;
del = del1 * del2 / del3;

eps1 = del1; 
eps2 = B*r_e*( E*atan(B*S_HAT) - 1 ) + (B^2)*E*r_e*S_HAT / (1 + (B*S_HAT)^2);
eps3 = del3 / U_HAT;
eps = eps1 * eps2 / eps3;

eta = Fz*sin( C*atan( B*(1 - E)*S_HAT + E*atan(B*S_HAT) ) );

% choose convenient constants
CONST1 = alp*del - beta*eps; 
CONST2 = exp(gam*ts*CONST1);

% Matrix Components
A11 = gam*( alp*del*CONST2 - beta*eps );
A12 = alp*gam*eps*( 1 - CONST2 );
A13 = -beta*gam*eta*( 1 - CONST2 );
A21 = -beta*gam*del*( 1 - CONST2 );
A22 = gam*( alp*del - CONST2*beta*eps );
A23 = -alp*gam*eta*( 1 - CONST2 );
A33 = gam*CONST1; 

A_pr = (1 / A33) * [A11 A12 A13; A21 A22 A23; 0 0 A33];

end

