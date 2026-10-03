function [Kframe, Gframe, vp, vs, ro2] = newcem(k, mu, x, ro1, P_cem)
%CONTACT_CEMENT_MODEL - Effective elastic moduli using Dvorkin's Contact Cement Model
%
% [Kframe, Gframe, vp, vs, ro2] = contact_cement_model(k, mu, x, ro1, P_cem)
%   k:      Bulk moduli vector. k(1) = Grain (Host), k(2) = Cement.
%   mu:     Shear moduli vector. mu(1) = Grain (Host), mu(2) = Cement.
%   x:      Fraction of phases. x(1) = solid fraction, x(2) = initial porosity (PhiC).
%   ro1:    Initial density of the dry, uncemented rock.
%   P_cem:  Fraction of pore space filled by cement.
%
% Note: Cementation scheme is fixed to Scheme 1 (cement at grain contacts).
% Gassmann fluid substitution is excluded (remaining pores remain blank/vacuum).

%% Extract parameters from inputs
K = k(1);
G = mu(1);
nu = (3*K - 2*G) / (2*(3*K + G)); % Poisson's ratio of host grain

Kc = k(2);
Gc = mu(2);
nuC = (3*Kc - 2*Gc) / (2*(3*Kc + Gc)); % Poisson's ratio of cement

PhiC = x(2);              % Critical/Initial porosity before cementation
Phi0 = PhiC - P_cem;      % Current porosity after cement is added
C = 8.5;                  % Constant coordination number

% Ensure cement volume is positive
if P_cem <= 0
    warning('P_cem is zero or negative. Contact cement model requires cement > 0.');
    P_cem = 1e-5; 
    Phi0 = PhiC - P_cem;
end

%% Scheme 1: Cement at Contacts
a = 2 * (((PhiC - Phi0) / (3 * C * (1 - PhiC)))^0.25);

% Capital Lambdas
alam = (2/pi) * (Gc/G) * (1-nu) * (1-nuC) / (1-2*nuC);
alamtau = (1/pi) * (Gc/G);

%% Calculate Effective Bulk Modulus (Kframe)
r1_K = Kc + (4 * Gc / 3);
r2_K = C * (1 - PhiC) / 6;
r3_K = -0.024153 * (alam^(-1.3646)) * (a^2) + 0.20405 * (alam^(-0.89008)) * a + 0.00024649 * (alam^(-1.9864));
Kframe = r1_K * r2_K * r3_K;

%% Calculate Effective Shear Modulus (Gframe)
r1_G = Gc;
r2_G = 3 * C * (1 - PhiC) / 20;

% Polynomial coefficients for shear calculation
a1t = -0.01 * (2.2606*nu^2 + 2.0696*nu + 2.2952);
a2t = 0.079011*nu^2 + 0.17539*nu - 1.3418;
b1t = 0.05728*nu^2 + 0.09367*nu + 0.20162;
b2t = 0.027425*nu^2 + 0.052859*nu - 0.87653;
c1t = 0.0001 * (9.6544*nu^2 + 4.9445*nu + 3.1008);
c2t = 0.018667*nu^2 + 0.4011*nu - 1.8186;

% Evaluate shear intermediate parameter r3
r3_G = (a1t * (alamtau^a2t)) * a^2 + (b1t * (alamtau^b2t)) * a + c1t * (alamtau^c2t);
Gframe = 0.6 * Kframe + r1_G * r2_G * r3_G;

%% Density and Velocities
ro_calcite = 2710; 
ro2 = ro1 + (P_cem * ro_calcite); % New density includes the added cement

% Calculate dry velocities
vp = sqrt((Kframe + (4/3)*Gframe) / ro2);
vs = sqrt(Gframe / ro2);

end