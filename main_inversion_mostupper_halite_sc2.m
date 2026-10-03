clearvars
close all
clc
Colors = brewermap(8,'Dark2');

%% ------------------------------------------------------------------------
%% Data
%% ------------------------------------------------------------------------
d =    [3.15; 1.71;   2300];  % data: vp, vs, rho_m             %given data
s =    [0.3;  0.05;    130];  % standard deviation of the noise 
H = [1;1;1];                 % select which data to use
nData = sum(H);
assigndatarun % Assuming this is a script you have in your path

%% ------------------------------------------------------------------------
%% Model (7 Parameters)
%% ------------------------------------------------------------------------
% Bounds: To get 0-2% cement and the rest gas, we let the sampler explore
% raw values, which will be mathematically normalized to 1 inside the models.
%     asp     phi     P_gas_raw   K      mu     rho_min  P_cem_raw
lb = [0.03    0.10    0.95        12.0   6.0    2677     0.00001   ]';  % Clay-like limits
ub = [0.99    0.50    1.00        37.0   44.0   2943     0.05  ]';  % Quartz-like limits
n = length(ub);                                               
 
%% ------------------------------------------------------------------------
%% Set up the MCMC and inversion
%% ------------------------------------------------------------------------
% Make a function for the log posterior 
% (Ensure this matches your actual function name, e.g., mylogpi_upper_calcite_sc1)
logpi = @(x)mylogpi_mostupper_halite_sc2(x,lb,ub,d,s,H);

%% Run the emcee hammer
Ne = 3*n;
Xo = zeros(n,Ne);
numModels = 0;

% Initialize with parameters that satisfy all constraints
go = 1;
counter = 0;
while go == 1
    xo = lb+(ub-lb).*rand(n,1);
    
    if isfinite(logpi(xo))
        counter = counter+1;
        Xo(:,counter) = xo;
        if counter>=Ne
            go = 0;
        end
    end
end

% ---------------- Cold start ----------------
Nsteps = 5e2;
[X, ~, LogPi,~] = myHammer(Nsteps,Xo,2.6,logpi,H);
Xrs = reshape(X,[n,size(X,2)*size(X,3)]);
LogPirs = reshape(LogPi,[1,size(LogPi,1)*size(LogPi,2)]);
RMSE = sqrt(2*LogPirs/length(d));
Xrs = Xrs(:,RMSE<3);
Xo = Xrs(:,randi(length(Xrs),Ne,1));

% ---------------- Warm start ----------------
Nsteps = 1e5;
[X, D, LogPi, AccRatio] = myHammer(Nsteps,Xo,2.6,logpi,nData);

% ---------------- MCMC aftermath ----------------
BurnIn = 1e4; 
% Safety check: If BurnIn is >= Nsteps, take the second half of the chain
if BurnIn >= Nsteps
    BurnIn = floor(Nsteps / 2); 
end

X = X(:,BurnIn:end,:);
D = D(:,BurnIn:end,:);
LogPi = LogPi(:,BurnIn:end);

% =========================================================================
% ENFORCE NORMALIZATION ON THE POSTERIOR ENSEMBLE
% This ensures your plotted results match the exact normalized physics 
% calculations performed inside myBerry.
% =========================================================================
sum_pores = X(3,:,:) + X(7,:,:);
X(3,:,:) = X(3,:,:) ./ sum_pores;  % Normalized P_gas
X(7,:,:) = X(7,:,:) ./ sum_pores;  % Normalized P_cem
% =========================================================================

Xrs = reshape(X,[n,size(X,2)*size(X,3)]);
Drs = reshape(D,[nData,size(D,2)*size(D,3)]);
LogPirs = reshape(LogPi,[1,size(LogPi,1)*size(LogPi,2)]);
RMSE = sqrt(2*LogPirs/length(d));

%% ------------------------------------------------------------------------
%% Plots (Your Custom Triangle/Histogram script)
%% ------------------------------------------------------------------------
newplotscript % Ensure your newplotscript has 7 labels for the TrianglePlot

%% ------------------------------------------------------------------------
%% Display numerical results
%% ------------------------------------------------------------------------
m = mean(Xrs,2);
sp = std(Xrs,[],2);
disp(' '), disp(' ')
fprintf('--- Inversion Results (Mean +/- Std) ---\n')
fprintf('Asp. ratio: %g +/- %g \n',m(1),sp(1))
fprintf('Porosity:   %g +/- %g \n',m(2),sp(2))
fprintf('P_gas:      %g +/- %g \n',m(3),sp(3))
fprintf('K (GPa):    %g +/- %g \n',m(4),sp(4))
fprintf('mu (GPa):   %g +/- %g \n',m(5),sp(5))
fprintf('rho_min:    %g +/- %g \n',m(6),sp(6))
fprintf('P_cem:      %g +/- %g \n',m(7),sp(7))
fprintf('----------------------------------------\n\n')

%% ------------------------------------------------------------------------
%% 1-D Ensemble Plot for Vp, Vs, and Density (Peak/Mode Plotting)
%% ------------------------------------------------------------------------
% 1. Extract the forward-modeled ensemble data
Vp_mod  = Drs(1, :);
Vs_mod  = Drs(2, :);
Rho_mod = Drs(3, :);

% 2. Identify the "Peak" Model (Mode of the PDF)
% This finds the highest peak of the probability distribution for the red line
[N_vp, edges_vp] = histcounts(Vp_mod, 100);
[~, max_idx_vp]  = max(N_vp);
peak_Vp          = (edges_vp(max_idx_vp) + edges_vp(max_idx_vp+1)) / 2;

[N_vs, edges_vs] = histcounts(Vs_mod, 100);
[~, max_idx_vs]  = max(N_vs);
peak_Vs          = (edges_vs(max_idx_vs) + edges_vs(max_idx_vs+1)) / 2;

[N_rho, edges_rho] = histcounts(Rho_mod, 100);
[~, max_idx_rho]   = max(N_rho);
peak_Rho           = (edges_rho(max_idx_rho) + edges_rho(max_idx_rho+1)) / 2;

% 3. Plotting Setup & Subsampling
% Subsample to 3000 models so the plot doesn't turn into a solid block of color
num_plot_models = min(3000, length(RMSE));
plot_indices = randperm(length(RMSE), num_plot_models);

% Sort by RMSE descending so the worst fits are in the background
plot_rmse = RMSE(plot_indices);
[~, sort_order] = sort(plot_rmse, 'descend');
sorted_indices = plot_indices(sort_order);

% Setup Colormap
cmap = colormap('parula');
min_c = min(plot_rmse);
max_c = max(plot_rmse);
getColor = @(rmse_val) cmap(max(1, round(255 * (rmse_val - min_c) / (max_c - min_c)) + 1), :);

% Depth constraints
z_top = 0;
z_bot = 2.1;
z_layer = [z_top, z_bot];
line_alpha = 0.75; 

% Initialize the Figure
figure('Name', '1-D Ensemble Plot', 'Position', [100, 100, 1300, 600], 'Color', 'w');

% =========================================================
% Subplot 1: P-Wave Velocity (Vp)
% =========================================================
subplot(1, 3, 1); hold on; grid on; box on;
for i = 1:length(sorted_indices)
    idx = sorted_indices(i);
    c = getColor(RMSE(idx));
    plot([Vp_mod(idx), Vp_mod(idx)], z_layer, 'Color', [c, line_alpha], 'LineWidth', 1, 'HandleVisibility', 'off');
end
% Plot the PEAK Model as Red
plot([peak_Vp, peak_Vp], z_layer, 'r-', 'LineWidth', 3, 'DisplayName', 'Peak Probability Model');
plot([d(1), d(1)], z_layer, 'k--', 'LineWidth', 2, 'DisplayName', 'Given Data');

yline([0.0-0.0,0.+0.0], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
yline([2.1-0.9,2.1+0.9], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
set(gca, 'YDir', 'reverse'); 
ylim([0, 11.5]);
xlabel('V_p (km/s)', 'FontWeight', 'bold');
ylabel('Depth (km)', 'FontWeight', 'bold');
title('1-D V_p Ensemble');
legend('Location', 'northeast');

% =========================================================
% Subplot 2: S-Wave Velocity (Vs)
% =========================================================
subplot(1, 3, 2); hold on; grid on; box on;
for i = 1:length(sorted_indices)
    idx = sorted_indices(i);
    c = getColor(RMSE(idx));
    plot([Vs_mod(idx), Vs_mod(idx)], z_layer, 'Color', [c, line_alpha], 'LineWidth', 1, 'HandleVisibility', 'off');
end
% Plot the PEAK Model as Red
plot([peak_Vs, peak_Vs], z_layer, 'r-', 'LineWidth', 3, 'DisplayName', 'Peak Probability Model');
plot([d(2), d(2)], z_layer, 'k--', 'LineWidth', 2, 'DisplayName', 'Given Data');

yline([0.0-0.0,0.+0.0], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
yline([2.1-0.9,2.1+0.9], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
set(gca, 'YDir', 'reverse');
ylim([0, 11.5]);
xlabel('V_s (km/s)', 'FontWeight', 'bold');
title('1-D V_s Ensemble');

% =========================================================
% Subplot 3: Density (Rho)
% =========================================================
subplot(1, 3, 3); hold on; grid on; box on;
for i = 1:length(sorted_indices)
    idx = sorted_indices(i);
    c = getColor(RMSE(idx));
    plot([Rho_mod(idx), Rho_mod(idx)], z_layer, 'Color', [c, line_alpha], 'LineWidth', 1, 'HandleVisibility', 'off');
end
% Plot the PEAK Model as Red
plot([peak_Rho, peak_Rho], z_layer, 'r-', 'LineWidth', 3, 'DisplayName', 'Peak Probability Model');
plot([d(3), d(3)], z_layer, 'k--', 'LineWidth', 2, 'DisplayName', 'Given Data');

yline([0.0-0.0,0.+0.0], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
yline([2.1-0.9,2.1+0.9], 'k:', 'LineWidth', 1, 'HandleVisibility', 'off'); 
set(gca, 'YDir', 'reverse');
ylim([0, 11.5]);
xlabel('Density, \rho (kg/m^3)', 'FontWeight', 'bold');
title('1-D Density Ensemble');