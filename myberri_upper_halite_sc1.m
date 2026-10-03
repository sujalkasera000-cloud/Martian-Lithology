function out = myberri_upper_halite_sc1(theta, H)
% MYBERRY - Objective function combining Berryman SC and Contact Cement Models
%
% INPUTS:
%   theta(1): Aspect ratio (for Berryman)
%   theta(2): Porosity (x_phi)
%   theta(3): P_gas_raw (Raw MCMC guess for gas fraction)
%   theta(4): Host Rock Bulk Modulus (GPa)
%   theta(5): Host Rock Shear Modulus (GPa)
%   theta(6): Host Rock Mineral Density (kg/m^3)
%   theta(7): P_cem_raw (Raw MCMC guess for cement fraction)
%   H:        Selection vector for outputs [Vp, Vs, Rho]

%% 1. Unpack Inversion Parameters (theta)
asp = [1 theta(1)];
x_phi = theta(2);
rock_vol = 1 - x_phi;
x = [rock_vol, x_phi]; 

mineral_density = theta(6);
rock_density = mineral_density * rock_vol; % Dry frame density of host rock

% --- THE NORMALIZATION TRICK ---
% Extract the raw MCMC guesses
raw_P_gas = theta(3);
raw_P_cem = theta(7);

% Force them to sum exactly to 1 before running the physics
total_fill = raw_P_gas + raw_P_cem;

% Safety check to prevent divide-by-zero if both happen to hit exactly 0
if total_fill == 0
    P_gas = 0.98;
    P_cem = 0.02;
else
    P_gas = raw_P_gas / total_fill;
    P_cem = raw_P_cem / total_fill;
end
% -------------------------------

% Fluid substitution parameters for Berryman
gas_density = 0.020 * x_phi;        
rhob1 = rock_density + gas_density; % Bulk density for Berryman input


%% 2. Run Berryman SC Model
% K and Mu inputs for Berryman (Solid Host and Fluid=0)
k_berry_in  = [theta(4)*1e9, 0];
mu_berry_in = [theta(5)*1e9, 0];

[kbr_berry, mubr_berry, ~, ~, ro2_berry, k2_berry] = berryscm_upper_halite_sc1(k_berry_in, mu_berry_in, asp, x, rhob1, P_gas);


%% 3. Run Contact Cement Model (Dvorkin - Scheme 1)
% Calcite Properties
K_halite = 25.2e9; 
G_halite = 15.3e9; 

% K and Mu inputs for Cement Model (Solid Host and Solid Cement)
k_cem_in  = [theta(4)*1e9, K_halite];
mu_cem_in = [theta(5)*1e9, G_halite];

% Pass rock_density as the initial dry frame density (ro1)
[Kframe_cem, Gframe_cem, ~, ~, ro2_cem] = cem_upper_halite_sc1(k_cem_in, mu_cem_in, x, rock_density, P_cem);


%% 4. Average the Moduli and Density
kbr  = (kbr_berry + Kframe_cem) / 2;
mubr = (mubr_berry + Gframe_cem) / 2;
k2   = (k2_berry + Kframe_cem) / 2;     % Fluids do not affect dry cement frame
ro2  = (ro2_berry + ro2_cem) / 2;


%% 5. Compute Final Averaged Velocities
vp = sqrt((k2 + (4/3)*mubr) ./ ro2);
vs = sqrt(mubr ./ ro2);
rhob = ro2; % Averaged bulk density for output


%% 6. Format Output based on 'H'
if sum(H) == 3
    out = [[vp; vs]./1e3; rhob];
elseif sum(H) == 2
    if H(1)==1 && H(2) == 1     %% vp and vs
         out = [vp; vs]./1e3;
    elseif H(1)==1 && H(3) == 1 %% vp and rho
        out = [vp/1e3; rhob];
    elseif H(2)==1 && H(3) == 1 %% vs and rho   
        out = [vs/1e3; rhob];
    else
        error('Invalid H vector configuration.')
    end
elseif sum(H) == 1
    if H(1) == 1
        out = vp/1e3;
    elseif H(2) == 1
       out = vs/1e3;
    elseif H(3) == 1
        out = rhob;
    else
        error('Invalid H vector configuration.')
    end
else
    error('H vector must contain at least one ''1''.')
end

end