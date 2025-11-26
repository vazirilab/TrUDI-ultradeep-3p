%% Monte-Carlo simulation to calculate 3p signal at focus as a function of depth and NA 
% Ballistic photons contributes to 3p signal 
% Author: Kalyan Santra

%%
clc;
rng(0)

%% Init parameters 
n0 = 1.33;
ls = 300e-6;
lambda = 1300e-9;
f_obj = 12.5e-3;
NA = [0.2, 0.4, 0.6, 0.8, 1.0];
w_BFA_list = f_obj * NA;
w0_focus = lambda ./ NA;

% Get percentage of ballistic photons at focus via MC
z0_list = [500, 1000, 1500, 2000, 2500, 3000] * 1e-6;
for NN = 1:length(w_BFA_list)
    for KK = 1:length(z0_list)
        z0 = z0_list(KK);
        w_beam = w_BFA_list(NN);
        Focal_percent(KK, NN) = MC_transport2focus(z0, w_beam, ls);
    end
end

%%
F_3p_ideal = (((exp(-z0_list'/ls) .* exp(+z0_list'/ls))).^3) * (1 ./ w0_focus.^2); % compensated 3p signal over depth
Norm_factor = F_3p_ideal(1, end);
F_3p_ideal_norm = F_3p_ideal ./ Norm_factor; % normalized 3p signal; without marginal attenuation effect

F_3p = (((Focal_percent .* exp(+z0_list'/ls))).^3) ./ w0_focus.^2;
F_3p_norm = F_3p ./ Norm_factor;


%% Make plots 
% 3p signal vs depth: without marginal attenuation effect
figure(2)
plot(z0_list'./1e-6, F_3p_ideal_norm)
xlabel('Focal depth (\mum)')
ylabel('3P signal (a.u.)')
legendStrings = "NA = " + string(NA);
legend(legendStrings, Box = "off")
title("3p signal without marginal beam attenuation")
ylim([0, 1.1])

% 3p signal vs depth: with marginal attenuation effect
figure(3)
plot(z0_list'./1e-6, F_3p_norm)
xlabel('Focal depth (\mum)')
ylabel('3P signal (a.u.)')
legend(legendStrings, Box = "off")
title("3p signal based on Monte-Carlo simulation")

if true
    filename = ['MC_parameter_and_results', datestr(now, '_yyyymmdd_HHMMSS')]; 
    save(filename)
end 

%%
function Focal_percent = MC_transport2focus(z0, w_beam, ls)
    % Parameters
    % NA = 0.25;
    % n0 = 1.33;
    % lambda = 1300e-9;
    f_obj = 12.5e-3;
    %z0 = 1000e-6;
    %g = 0.9;
    %ls = 300e-6;
    N_photons = 1e7; % photon packet 
    %w_beam = 2.0e-3;
    % Histogram grid at focal plane
    xbins = linspace(-1e-6, 1e-6, 100);
    ybins = linspace(-1e-6, 1e-6, 100);

    photon_hits = zeros(N_photons, 2); % x, y positions
    photon_count = false(N_photons, 1); % valid photon mask

    % Sample BFP positions (Gaussian)
    r0 = sqrt(-0.5*w_beam^2*log(rand(N_photons, 1)));
    theta0 = 2 * pi * rand(N_photons, 1);
    x0 = r0 .* cos(theta0);
    y0 = r0 .* sin(theta0);

    parfor i = 1:N_photons
        % Prop. BFP to lens (z = f)
        r_bfp = [x0(i), y0(i), 0];
        u1 = [0, 0, 1];
        r_lens = r_bfp + f_obj * u1;

        % Get ray direction after lens to focus at z = 2f
        focus = [0, 0, 2 * f_obj - 0e-6];
        delta = focus - r_lens;
        u = delta / norm(delta);

        % Enter tissue at z = 2f - z0
        z_entry = 2 * f_obj - z0;
        s_entry = (z_entry - r_lens(3)) / u(3);
        r = r_lens + s_entry * u;

        % Get travel distance before first scattering
        s = -ls * log(rand());
        r_new = r + s * u;

        % Get unscattered photons
        if r_new(3) >= 2 * f_obj
            % For unscattered photon backtrack to focal plane
            s_back = (2 * f_obj - r(3)) / u(3);
            r_focus = r + s_back * u;
            photon_hits(i, :) = r_focus(1:2);
            photon_count(i) = true;
        end
    end
    Focal_percent = sum(photon_count) / N_photons;

    % Bin all valid photons
    xvals = photon_hits(photon_count, 1);
    yvals = photon_hits(photon_count, 2);
    E_focus = histcounts2(yvals, xvals, ybins, xbins);

    % Plot
    imagesc(xbins*1e6, ybins*1e6, E_focus);
    axis image
    xlabel('x (μm)');
    ylabel('y (μm)');
    %title(sprintf('Photon density at focus (z = %.1f µm)', 2*f_obj*1e6));
    colorbar
end

