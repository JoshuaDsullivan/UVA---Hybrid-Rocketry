%% tank_optimizer_tube.m
% Optimize a hollow tube (annular cylinder) tank:
% - uses mean diameter for thin-wall hoop stress
% - computes annular shell volume (outer minus inner)
% - performs bolt checks (bearing, tear, shear)
% - optimizes over stocked OD & wall options, bolt diameters, bolt count, with length constrained to 0-1 m
%
% Units: SI (m, Pa, kg). Change inputs below as needed.

clear; close all; clc;

%% ---------------- User inputs ----------------
% Given values
sigma_y_Al = 3.1*10^8;       % Pa (Al-6061 T6 yield)
tank_pressure = 6894760;      % Pa (internal pressure)
required_volume = 0.001819536;% m^3 (total tank volume required)
MinFS = 2;                    % minimum factor of safety required (hoop + bolts)
rho_al = 2700;                % kg/m^3 (density of aluminum)

% Bolt variables (user-changeable)
n_list = [2 3 4 5 6 8 10 11 12 14];  % bolts per end (vector to search)
dis = 0.003;  % m -- distance (from bolt row to tank top). set as desired.

% Bolt diameters to search (inches -> meters)
bolt_diam_in = [1/8, 3/16, 1/4];    % common small bolt sizes
bolt_diam_m = bolt_diam_in * 0.0254;

% Standard stocked tube OD (mm) and wall (mm) lists -- replace with vendor list if desired
OD_in = [3.75];
%OD_in = [2.75, 3, 3.25 ];   % outer diameters (mm)
wall_in = [.25];
%[.065, .125, .25, .375];                  % common walls (mm)

OD = OD_in*0.0254;  % m
wall = wall_in*0.0254;  % m
% mean diameter for thin-wall hoop stress (now replaced with thick-wall)

%% ---------- material properties for bolts (approx) ----------
% For shear capacity of steel bolts (used for FS_shear)
sigma_y_steel = 1.1721*10^9;   % Pa (adjust for actual bolt grade)
tau_st = 0.577 * sigma_y_steel;  % conservative conversion

%% ---------- results storage ----------
res = struct('OD',[],'t',[],'ID',[],'Dmean',[],'L',[],'bolt_d',[],'n',[],'mass_shell',[],'sigma_hoop',[],'FS_hoop',[],'FS_bearing',[],'FS_tear',[],'FS_shear',[]);
idx = 0;

%% ---------- Search loop ----------
for iOD = 1:length(OD)
    for it = 1:length(wall)
        OD_i = OD(iOD);
        t_i  = wall(it);
        ID_i = OD_i - 2*t_i;        % inner diameter
        r_i = ID_i / 2;
        r_o = OD_i / 2;
        if ID_i <= 0
            continue;               % invalid geometry
        end

        % compute required length for the annular tube (use cross-sectional area of the inner fluid)
        % NOTE: required_volume is total internal volume (fluid ullage included).
        % If the user provided "Total Tank Volume" as external volume, adjust accordingly.
        L_i = required_volume / (pi*(ID_i/2)^2);
        if L_i <= 0 || L_i > 1.0
            continue;  % enforce length range 0-1 m
        end

        % mean diameter for thin-wall hoop stress
        D_mean = (OD_i + ID_i) / 2;

        % thin-wall hoop stress: sigma_theta = p * D_mean / (2 * t)
        sigma_hoop = tank_pressure * (r_o^2 + r_i^2) / (r_o^2 - r_i^2);
        FS_hoop = sigma_y_Al / sigma_hoop;

        % require hoop FS
        if FS_hoop < MinFS
            continue;
        end

        % annular shell volume and mass
        vol_shell = pi*(OD_i^2 - ID_i^2)/4 * L_i;   % m^3
        mass_shell = vol_shell * rho_al;

        % Internal pressure axial force on endcap using inner diameter area
        endcap_area = pi*(ID_i/2)^2;
        F_total = tank_pressure * endcap_area;   % N, net axial separating force

        % iterate bolt sizes and counts
        for ib = 1:length(bolt_diam_m)
            b_d = bolt_diam_m(ib);
            A_bolt = pi*(b_d/2)^2;   % bolt shear area (gross)

            for in_idx = 1:length(n_list)
                n_b = n_list(in_idx);
                if n_b <= 0
                    continue;
                end

                % Force distribution per bolt (axial)
                F_per_bolt = F_total / n_b;

                % Bearing stress on tube wall at bolt hole (approx):
                % sigma_bearing = F_per_bolt / (b_d * t)
                % where b_d is bolt diameter (bearing width) and t is tube wall thickness.
                sigma_bearing = F_per_bolt / (b_d * t_i);
                FS_bearing = sigma_y_Al / sigma_bearing;

                % Tear (tensile) stress approximation (your original formulation):
                % sigma_t = F_total / (2 * n * dis * t)
                % Note: dis is radial lever distance from bolt row to axis (m). Ensure dis > 0.
                sigma_t = F_total / (2 * n_b * dis * t_i);
                tau_conv = 0.577 * sigma_y_Al;
                FS_tear = tau_conv / sigma_t;

                % Shear on bolt (single shear assumption)
                sigma_s = F_per_bolt / A_bolt;
                FS_shear = tau_st / sigma_s;

                % All checks must meet MinFS
                if FS_bearing >= MinFS && FS_tear >= MinFS && FS_shear >= MinFS
                    idx = idx + 1;
                    res(idx).OD = OD_i;
                    res(idx).t = t_i;
                    res(idx).ID = ID_i;
                    res(idx).Dmean = D_mean;
                    res(idx).L = L_i;
                    res(idx).bolt_d = b_d;
                    res(idx).n = n_b;
                    res(idx).mass_shell = mass_shell;
                    res(idx).sigma_hoop = sigma_hoop;
                    res(idx).FS_hoop = FS_hoop;
                    res(idx).FS_bearing = FS_bearing;
                    res(idx).FS_tear = FS_tear;
                    res(idx).FS_shear = FS_shear;
                end
            end
        end
    end
end

%% ---------- Report ----------
if idx == 0
    fprintf('No feasible tube designs found. Try larger wall thicknesses, larger OD, higher n, or adjust dis.\n');
    return;
end

T = struct2table(res);
T_sorted = sortrows(T,'L');

n_print = min(20,height(T_sorted));
fprintf('Top %d feasible tube designs (sorted by shell mass):\n', n_print);
fprintf('OD(in)  t(in)  ID(in)  Dmean(mm)  L(mm)  bolt_d(in)  n  mass_shell(kg)  FS_hoop  FS_bear  FS_tear  FS_shear\n');
for k = 1:n_print
    row = T_sorted(k,:);
    fprintf('%6.2f  %5.2f  %6.2f  %8.2f  %6.1f   %6.3f      %2d  %10.4f      %5.2f    %5.2f    %5.2f    %5.2f\n', ...
        row.OD/0.0254, row.t/0.0254, row.ID/0.0254, row.Dmean*1000, row.L*1000, row.bolt_d/0.0254, row.n, row.mass_shell, row.FS_hoop, row.FS_bearing, row.FS_tear, row.FS_shear);
end

tank_options = T_sorted;  % full table returned to workspace

% Quick plot
figure;

scatter(T_sorted.Dmean*1000, T_sorted.mass_shell, 40, 'filled');
xlabel('Mean Diameter (mm)'); ylabel('Shell mass (kg)');
title('Feasible tube designs: shell mass vs mean diameter');
grid on;
