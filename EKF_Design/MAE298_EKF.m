%% Non-linear System
syms U w mu % States
syms T % Inputs, which can be replaced with actual values later
syms C B E % Macro-parameters, which can be replaced with actual values later
syms m Fz J re % Tire parameters, which can be replaced with actual values later
s = (re*w/U)-1;
Fx = Fz*mu*sin(C*atan(B*s-E*(B*s*atan(B*s))));
f = [Fx/(m/4); (T-re*Fx)/J; 0];
%% Continuous-Time system
A_c = sym(zeros(3,3)); x = [U;w;mu];
for i = 1:3
    for j = 1:3
    A_c(i,j) = diff(f(i),x(j));  
    end
end; clear i j
B_c = [0 1/J 0];
C_c = [1 0 0; 0 1 0];
D_c = [0;0];
%% Discrete-Time system
%%%% After applying actual values, the discrete-time system can be found
T_s = 1;% Sampling period, which can be replaced with actual values later
%% Noise
Q_k = ;
R_k = ;
%% EKF
t = ; % Time
x_est = zeros(1,length(t));
for i = 2:length(t)
%%%% Put the evaluation process using function subs()
F = expm(A_c*dt);
G = integral(@(t) expm(A_c.*t),0,dt, 'ArrayValued', true); G = G*B_c;
end
