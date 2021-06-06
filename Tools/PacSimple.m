function [F_x,F_y,G_xalpha,G_ykappa] = PacSimple(Tire,F_z,alpha,kappa)
%PACSIMPLE Implementation of the Simplified Pacejka Tire Model
% Inputs:
%   Tire.mu_max: maximum coefficient of friction
%   Tire.B_x: Longitudinal Stiffness Factor
%   Tire.C_x: Longitudinal Shape Factor
%   Tire.E_x: Longitudinal Curvature Factor
%   Tire.B_xalpha: Combined Stiffness Factor
%   Tire.C_xalpha: Combined Shape Factor
%   Tire.E_xalpha: Combined Curvature Factor
%   Tire.B_y: Lateral Stiffness Factor
%   Tire.C_y: Lateral Shape Factor
%   Tire.E_y: Lateral Curvature Factor
%   Tire.B_ykappa: Combined Stiffness Factor
%   Tire.C_ykappa: Combined Shape Factor
%   Tire.E_ykappa: Combined Curvature Factor
%   F_z: Normal Tire Force [N]
%   alpha: Lateral Slip Angle [rad]
%   kappa: Longitudinal Slip Ratio [-1,1]
% Outputs:
%   F_x: Longitudinal Force [N]
%   F_y: Lateral Force [N]
%   G_xalpha: Longitudinal Weighting Factor (for Debugging/parameter fit)
%   G_ykappa: Lateral Weighting Factor (for Debugging/paramter fit)

% mu_max = Tire.mu_max;
D_x = Tire.D_x; D_y = Tire.D_y;
B_x = Tire.B_x;C_x = Tire.C_x;E_x = Tire.E_x;
B_xalpha = Tire.B_xalpha;C_xalpha = Tire.C_xalpha; E_xalpha = Tire.E_xalpha;
B_y = Tire.B_y;C_y = Tire.C_y;E_y = Tire.E_y;
B_ykappa = Tire.B_ykappa;C_ykappa = Tire.C_ykappa;S_Hykappa = Tire.S_Hykappa;

% Longitudinal Slip:
% kappa = -(v_wx - omega*r_c)/v_wx;
% Nominal Longitudinal Tire/Road Contact Force:
F_x0 = D_x*sin(C_x*atan(B_x*kappa - E_x*(B_x*kappa - atan(B_x*kappa))));
% Longitudinal Force Weighting Function:
G_xalpha = cos(C_xalpha*atan(B_xalpha*alpha - E_xalpha*(B_xalpha*alpha - atan(alpha))));
% Longitudinal Tire/Road Contact Force:
F_x = G_xalpha*F_x0;
% Nominal Lateral Tire/Road Contact Force:
F_y0 = D_y*sin(C_y*atan(B_y*alpha - E_y*(B_y*alpha - atan(B_y*alpha))));
% Lateral Force Weighting Function:
G_ykappa = cos(C_ykappa*atan(B_ykappa*(kappa + S_Hykappa)))/...
            cos(C_ykappa*atan(B_ykappa*S_Hykappa));
% Lateral Tire/Road Contact Force:
F_y = G_ykappa*F_y0;
end

