%% Dust Cloud Particle Simulation
% n particles, n satellites
clear all
close all
clc

%% Constants

obj = [];
obj.G = 6.67430e-11;                       % [m^3/kg/s^2]
obj.mass_Sun = 1.988400e30;                % [kg]
obj.mass_Earth = 5.9722e24;                % [kg]
obj.distance_Sun_Earth = 1.495978707e11;   % [m]
obj.distance_Sun_L1 = 1.480018707e11;      % [m]
obj.speed_light = 2.99792458e8;            % [m/s]
obj.solar_luminosity = 3.846e26;           % [W]
obj.omega_Earth = 2*pi/(365.25*86400);     % [rad/s]

% Sun fixed at origin
Sun_pos = [0; 0; 0];

% Earth fixed at +x axis
Earth_pos = [obj.distance_Sun_Earth; 0; 0];   % in km for the plot

%% Spacecraft parameters
obj.epsilon_0 = 8.854e-12;      % [F/m]
obj.Voltage_SC = 100000;        % [V]
obj.L_SC = 2 * 10e3;            % tether length [m]
obj.a_SC = 0.001;               % [m]
obj.mass_tether = obj.L_SC * (2*pi* obj.a_SC^2) * 8960;   % [kg]

obj.Capacitance_SC = 2*pi* obj.epsilon_0* obj.L_SC / log(2* obj.L_SC/ obj.a_SC);
obj.Charge_SC = obj.Voltage_SC * obj.Capacitance_SC;     % [C]

% Dust particle properties
obj.radius_particle = 3e-7; % Particle radius [m]
obj.density_particle = 220;     % [kg/m^3]
obj.cross_sectional_area = pi * obj.radius_particle^2; % Cross-sectional area [m^2]

obj.volume_particle = (4/3)*pi*(obj.radius_particle^3);
obj.mass_particle = obj.volume_particle * obj.density_particle;

obj.Charge_to_mass_ratio = 0.784; % [C kg−1]
% 90% hollow silica particles of 0.3micron radius is 13.41
% Fille silica particles of 0.3micron radius is 0.784
obj.particle_charge = obj.Charge_to_mass_ratio * obj.mass_particle; % Particle charge [C]

% Q_sca Computed by Mark [Particle Radius [μm]	Spectrally averaged Q_sca]

obj.g_sca_asym = 0.6; % Set to 0.6 for scattering, and 0 for absorbing

given_particle_radius_Qsca_array = [0.1     0.007176244150242173;
    0.2     1.4934342934524942;
    0.3     2.4466003650566077;
    1.0     2.2854053039113142;
    1.25    1.9253989635168711;
    2.5     0.85394165707732];

obj.Q_sca = interp1(given_particle_radius_Qsca_array(:,1),given_particle_radius_Qsca_array(:,2), round(obj.radius_particle * 1e6,2), 'linear');
obj.Q_rad_pressure_coefficient = obj.Q_sca * (1 - obj.g_sca_asym); % Radiation efficiency factor

%% SPICE Initialization

flag_file_location = 3;

switch flag_file_location

    case 3
        path_to_MuSCAT_Supporting_Files = '../MuSCAT_Supporting_Files/';
        path_to_MuSCAT_v2 = '../MuSCAT_Matlab_v2/';

    otherwise
        error('Should not reach here!')
end

addpath(genpath(path_to_MuSCAT_Supporting_Files));
addpath(genpath(path_to_MuSCAT_v2));
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/de440s.bsp']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/naif0012.tls']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/pck00011.tpc']);

%% Hexagonal honeycomb in y-z plane

hex_rings = 59;           % rings around center -> 8267 sat. = 52 rings
Columb_interSC_distance = 39000; % spacing [m]
d = Columb_interSC_distance;     % radial spacing

points = [];
for q = -hex_rings:hex_rings
    r1 = max(-hex_rings, -q - hex_rings);
    r2 = min(hex_rings, -q + hex_rings);
    for r = r1:r2
        y = d * (sqrt(3) * q + sqrt(3)/2 * r);
        z = d * (3/2 * r);
        points = [points; y z];
    end
end

honeycomb_radius = max( sqrt( points(:,1).^2 + points(:,2).^2 ) );

x_L1 = obj.distance_Sun_L1;
obj.num_Columb_SC = size(points,1);
Columb_Position_SC_array = [ x_L1 * ones(obj.num_Columb_SC,1), points(:,1), points(:,2) ]; % Nx3: [x y z]

edge_point_idx = logical(vecnorm(points,2,2) > honeycomb_radius*sind(60));


% Honeycomb y-z figure
plot_handle = figure();
clf
set(plot_handle, 'Name','Hexagonal Honeycomb (y-z plane)')
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','outerposition',[0 0 0.5 0.5])
set(plot_handle,'PaperPositionMode','auto');
% figure('Name','Hexagonal Honeycomb (y-z plane)','Color',[1 1 1]);
hold on
scatter(points(:,1), points(:,2), 60, 'filled');
scatter(points(edge_point_idx,1), points(edge_point_idx,2), 'r', 'LineWidth',2);
axis equal; grid on;
xlabel('y [m]'); ylabel('z [m]');
title(sprintf('Hexagonal Honeycomb at x = %.3e m, (no. of SC = %d)', x_L1, obj.num_Columb_SC));
set(gca,'FontSize',14, 'FontName','Times New Roman')

% Remove edge points
points = points(~edge_point_idx,:);
obj.num_Columb_SC = size(points,1);
Columb_Position_SC_array = [ x_L1 * ones(obj.num_Columb_SC,1), points(:,1), points(:,2) ]; % Nx3: [x y z]

inscribed_circle_radius = max(vecnorm(points, 2, 2));
edge_point_idx = logical(vecnorm(points,2,2) >= inscribed_circle_radius - 2*Columb_interSC_distance);

obj.num_Columb_SC_base = size(Columb_Position_SC_array,1);

radius_dust_cloud = 2.6706e+06; % [m]
if radius_dust_cloud > 0.8*inscribed_circle_radius
    warning('radius_dust_cloud > 0.8*inscribed_circle_radius!')
    radius_dust_cloud
    0.8*inscribed_circle_radius
end

% Honeycomb y-z figure
plot_handle = figure();
clf
set(plot_handle, 'Name','Circle with Hexagonal Honeycomb (y-z plane)')
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','outerposition',[0 0.5 0.5 0.5])
set(plot_handle,'PaperPositionMode','auto');
% figure('Name','Hexagonal Honeycomb (y-z plane)','Color',[1 1 1]);
hold on
scatter(points(:,1), points(:,2), 60, 'filled');
scatter(points(edge_point_idx,1), points(edge_point_idx,2), 'r', 'LineWidth',2);
axis equal; grid on;
xlabel('y [m]'); ylabel('z [m]');
title(sprintf('Circle with Hexagonal Honeycomb at x = %.3e m, (no. of SC = %d)', x_L1, obj.num_Columb_SC));
set(gca,'FontSize',14, 'FontName','Times New Roman')


% Step Forward
num_forward_steps = 20;
factor_forward_step = 0.02;
old_Columb_Position_SC_array = Columb_Position_SC_array;

for i_fs = 1:1:num_forward_steps

    new_Columb_Position_SC_array = old_Columb_Position_SC_array(edge_point_idx,:) - [factor_forward_step*i_fs*Columb_interSC_distance, 0, 0]; % [m]

    new_Columb_Position_SC_array(:,2:3) = (1 + 0*i_fs/10)*new_Columb_Position_SC_array(:,2:3);

    Columb_Position_SC_array = [Columb_Position_SC_array; new_Columb_Position_SC_array];
end

% Step Backward
num_bands = 10;
for i_b = 1:1:num_bands

    min_dist = (i_b-1)*(inscribed_circle_radius/num_bands);
    max_dist = (i_b)*(inscribed_circle_radius/num_bands);

    band_point_idx = logical(vecnorm(Columb_Position_SC_array(:,2:3),2,2) >= min_dist) & logical(vecnorm(Columb_Position_SC_array(:,2:3),2,2) < max_dist);

    Columb_Position_SC_array(band_point_idx,1) = Columb_Position_SC_array(band_point_idx,1) + (num_bands-i_b)*Columb_interSC_distance;

end

% % Add opposing SC
% old_Columb_Position_SC_array = old_Columb_Position_SC_array - [factor_forward_step*num_forward_steps*interSC_distance + 500000, 0, 0]; % [m]
% new_opposing_Columb_Position_SC_array = [];
%
% num_bands = 10;
% for i_b = 1:1:num_bands
%
%     min_dist = (i_b-1)*(inscribed_circle_radius/num_bands);
%     max_dist = (i_b)*(inscribed_circle_radius/num_bands);
%
%     band_point_idx = logical(vecnorm(old_Columb_Position_SC_array(:,2:3),2,2) >= min_dist) & logical(vecnorm(old_Columb_Position_SC_array(:,2:3),2,2) < max_dist);
%
%     old_Columb_Position_SC_array(band_point_idx,1) = old_Columb_Position_SC_array(band_point_idx,1) - 10*(num_bands-i_b)*interSC_distance;
%
%     old_Columb_Position_SC_array(band_point_idx,2:3) = 1*old_Columb_Position_SC_array(band_point_idx,2:3);
%
%     if mod(i_b,4) == 0
%         new_opposing_Columb_Position_SC_array = [new_opposing_Columb_Position_SC_array; old_Columb_Position_SC_array(band_point_idx,:)];
%     end
%
% end
%
% % Columb_Position_SC_array = [Columb_Position_SC_array; old_Columb_Position_SC_array];
% Columb_Position_SC_array = [Columb_Position_SC_array; new_opposing_Columb_Position_SC_array];


obj.num_Columb_SC = size(Columb_Position_SC_array,1);

% Mirror Spacecraft
Mirror_interSC_distance = 500; % spacing [m]
num_forward_steps = 50;
obj.num_Mirror_SC_per_step = ceil(2*pi*inscribed_circle_radius/Mirror_interSC_distance);
delta_theta = 360/obj.num_Mirror_SC_per_step; % [deg]
theta_array = [delta_theta:delta_theta:360]';

obj.num_Mirror_SC = 0;
Mirrior_Position_SC_array = [];
Mirrior_Position_SC_orientation_unit_vector = [];
% for i_steps = 1:1:num_forward_steps
%
%     this_Mirrior_Position_SC_array = [ (x_L1 - 5*(i_steps-1)*Mirror_interSC_distance ) * ones(obj.num_Mirror_SC_per_step,1), inscribed_circle_radius*sind(theta_array), inscribed_circle_radius*cosd(theta_array) ];
%
%     Mirrior_Position_SC_array = [Mirrior_Position_SC_array; this_Mirrior_Position_SC_array];
% end

% obj.num_Mirror_SC = size(Mirrior_Position_SC_array,1);
%
% Mirrior_Position_SC_orientation_unit_vector = [Mirrior_Position_SC_array(:,1), zeros(obj.num_Mirror_SC,2)] - Mirrior_Position_SC_array;
% Mirrior_Position_SC_orientation_unit_vector = Mirrior_Position_SC_orientation_unit_vector./vecnorm(Mirrior_Position_SC_orientation_unit_vector,2,2);

% See equations from here: https://en.wikipedia.org/wiki/Solid_angle
obj.Mirror_cone_half_apex_angle = 1; % [deg]
obj.Mirror_cone_solid_angle = 4*pi*(sind(obj.Mirror_cone_half_apex_angle))^2; % [steradian]
obj.Mirror_SC_area = 100*100; % [m^2]

plot_handle = figure();
clf
set(plot_handle, 'Name','All SC Position')
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','outerposition',[0.5 0 0.5 1])
set(plot_handle,'PaperPositionMode','auto');

hold on; grid on;

plot3(Columb_Position_SC_array(:,1)/1e3, Columb_Position_SC_array(:,2)/1e3, Columb_Position_SC_array(:,3)/1e3, 'or', 'MarkerSize',5 );

if ~isempty(Mirrior_Position_SC_array)
    plot3(Mirrior_Position_SC_array(:,1)/1e3, Mirrior_Position_SC_array(:,2)/1e3, Mirrior_Position_SC_array(:,3)/1e3, 'ob', 'MarkerSize',5 );
end

title(sprintf('Circle with Hexagonal Honeycomb and Step Forward at x = %.3e m, (Columb SC = %d, Mirror SC = %d)', x_L1, obj.num_Columb_SC, obj.num_Mirror_SC));
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
view(3)
set(gca,'FontSize',14, 'FontName','Times New Roman')

drawnow limitrate


%% Generate Random Dust Particle Initial Positions (in y-z plane)

obj.threshold_particle_base_distance = 10000; % m

obj.first_dynamic_step = 1*obj.threshold_particle_base_distance;
obj.second_dynamic_step = 1*obj.threshold_particle_base_distance;


obj.num_particles = 1;

theta = 2*pi*rand(obj.num_particles,1);
r_rand = 0.8*inscribed_circle_radius * rand(obj.num_particles,1);

y_rand = r_rand .* cos(theta);
z_rand = r_rand .* sin(theta);

x_rand = (obj.distance_Sun_L1 - obj.threshold_particle_base_distance) * ones(obj.num_particles,1);

positions_particles = [x_rand, y_rand, z_rand];
velocities_particles = zeros(obj.num_particles, 3);

Particle_Pos_Vel_init_array = [positions_particles, velocities_particles];

% save('Particle_Pos_Vel_init_array.mat','Particle_Pos_Vel_init_array')
% load('Particle_Pos_Vel_init_array.mat')

% y0 = zeros(6 * obj.num_particles, 1);
%
% for p = 1:obj.num_particles
%     idx = (p-1)*6 + (1:6);
%     y0(idx) = [ positions_particles(p,:)'; velocities_particles(p,:)' ];
% end

%% Time span
time_step = 10*60; % [sec]
tspan = [0:time_step:10*24*60*60]; % sec

Particle_Pos_Vel_full_array = zeros(length(tspan),6,obj.num_particles);

k = 1;
for p = 1:obj.num_particles
    Particle_Pos_Vel_full_array(k,:,p) = Particle_Pos_Vel_init_array(p,:);
end


%% Integration for all particles

for k = 2:1:length(tspan)

    disp(['Time ',num2str(tspan(k)),' sec'])

    y0 = [];
    for p = 1:obj.num_particles
        y0 = [y0; Particle_Pos_Vel_full_array(k-1,:,p)'];
    end

    [t, Y] = ode45(@(t,Y) dustDynamics_all_particles(t, Y, obj, Columb_Position_SC_array, Mirrior_Position_SC_array, Mirrior_Position_SC_orientation_unit_vector), [0:1:time_step], y0);

    y_end = Y(end,:);

    for p = 1:obj.num_particles
        Particle_Pos_Vel_full_array(k,:,p) = y_end((p-1)*6+1:(p-1)*6+6);
    end



end

% %% Extract trajectories (m)
% positions = zeros(length(t), obj.num_particles, 3);
%
% for p = 1:obj.num_particles
%     idx = (p-1)*6 + (1:3);
%     positions(:, p, :) = Y(:, idx);
% end

%% ALL PARTICLE TRAJECTORIES

disp('Plotting now!')

plot_handle = figure();
clf
set(plot_handle, 'Name','Particle Motion')
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','outerposition',[0 0 0.5 1])
set(plot_handle,'PaperPositionMode','auto');

hold on; grid on;

plot3(Columb_Position_SC_array(:,1)/1e3, Columb_Position_SC_array(:,2)/1e3, Columb_Position_SC_array(:,3)/1e3, 'or', 'MarkerSize',5 );

if ~isempty(Mirrior_Position_SC_array)
    plot3(Mirrior_Position_SC_array(:,1)/1e3, Mirrior_Position_SC_array(:,2)/1e3, Mirrior_Position_SC_array(:,3)/1e3, 'ob', 'MarkerSize',5 );
end


xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');

for p = 1:obj.num_particles

    y_km = Particle_Pos_Vel_full_array(:,1:3,p)/1e3;
    % Sun_pos_km = [0 0 0];
    % Earth_pos_km = [distance_Sun_Earth/1000, 0, 0];

    plot3(y_km(:,1), y_km(:,2), y_km(:,3), '-k', 'LineWidth', 2);
    hold on

    % Plot particle start and end positions with large markers
    plot3(y_km(1,1), y_km(1,2), y_km(1,3), 'og', 'MarkerSize', 10, 'MarkerFaceColor','g');  % start
    plot3(y_km(end,1), y_km(end,2), y_km(end,3), 'ob', 'MarkerSize', 10, 'MarkerFaceColor','b'); % end

end

% x_min = -1e7; x_max = 2e8;   % km
% y_min = -5e7; y_max = 5e7;   % km
% z_min = -5e7; z_max = 5e7;   % km
% xlim([x_min x_max]);
% ylim([y_min y_max]);
% zlim([z_min z_max]);
%
% light('Position',[1 0 0],'Style','infinite');
% material dull;
set(gca,'FontSize',14, 'FontName','Times New Roman')
% legend({'Sun','Earth','Dust Particle'}, 'Location','best');
% view(45,20);

view(3)

drawnow limitrate
% saveas(gcf, 'Sun_Earth_Dust_Trajectory.png');

%% X-Position vs Time

plot_handle = figure();
clf
set(plot_handle, 'Name','Particle Motion')
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','outerposition',[0.5 0 0.5 1])
set(plot_handle,'PaperPositionMode','auto');

hold on; grid on;

for p = 1:obj.num_particles

    plot(tspan/(60*60),Particle_Pos_Vel_full_array(:,1,p)/1e3,'-')

end

plot([tspan(1) tspan(end)]/(60*60),1e-3*(obj.distance_Sun_L1)*[1 1],'-k','LineWidth',2)
plot([tspan(1) tspan(end)]/(60*60),1e-3*(obj.distance_Sun_L1 - obj.first_dynamic_step)*[1 1],'-r','LineWidth',2)
plot([tspan(1) tspan(end)]/(60*60),1e-3*(obj.distance_Sun_L1 - obj.first_dynamic_step - obj.second_dynamic_step)*[1 1],'--r','LineWidth',2)

xlabel('Time [hours]');
ylabel('Position in X axis [km]');
set(gca,'FontSize',14, 'FontName','Times New Roman')




%% Function for Equations of Motion
function dYdt = dustDynamics_all_particles(t, Y, obj, Columb_Position_SC_array, Mirrior_Position_SC_array, Mirrior_Position_SC_orientation_unit_vector)

% persistent lastPrintTime
% if isempty(lastPrintTime)
%     lastPrintTime = -1e9;
% end

dYdt = zeros(6*obj.num_particles,1);

r1_max = max(Y(1:6:end));
flag_dynamic_case = 0;
% Dynamic check
if (obj.distance_Sun_L1 - r1_max <= obj.first_dynamic_step)
    % Keep Base on
    flag_dynamic_case = 1;

elseif (obj.distance_Sun_L1 - r1_max <= obj.first_dynamic_step + obj.second_dynamic_step)
    % switch off base proportionately
    flag_dynamic_case = 2;
    this_gain = 1 - ((obj.distance_Sun_L1 - r1_max - obj.first_dynamic_step)/obj.second_dynamic_step);

else % (obj.distance_Sun_L1 - r(1) > 3*obj.threshold_particle_base_distance)
    % switch off base completely
    flag_dynamic_case = 3;

end
% flag_dynamic_case


for p=1:1:obj.num_particles

    % Loop over each particle (n particles), compute forces including Coulomb
    r = Y((p-1)*6+1:(p-1)*6+3);
    v = Y((p-1)*6+4:(p-1)*6+6);
    r_mag = norm(r);

    % Solar Radiation Pressure
    F_SRP_mag = (obj.solar_luminosity * obj.cross_sectional_area * obj.Q_rad_pressure_coefficient) / (4*pi* obj.speed_light * r_mag^2);
    F_SRP_vec = F_SRP_mag * (r / r_mag);

    % Gravitational Forces
    F_grav_Sun   = -obj.G * obj.mass_Sun  * obj.mass_particle * (r / r_mag^3);
    r_rel_earth = r - [obj.distance_Sun_Earth; 0; 0];
    r_rel_earth_norm = norm(r_rel_earth);
    F_grav_Earth = -obj.G * obj.mass_Earth * obj.mass_particle * (r_rel_earth / r_rel_earth_norm^3);

    % Centripetal
    F_centripetal = obj.mass_particle * obj.omega_Earth^2 * r;

    % Coulomb: sum contributions from all satellites in Columb_Position_SC_array
    F_Coulomb = [0;0;0];
    % for k = 1:size(Columb_Position_SC_array,1)
    %     R = r - Columb_Position_SC_array(k,:)';
    %     d = norm(R);
    %
    %     % if d == 0
    %     %     continue;
    %     % end
    %     % Coulomb force: (1/(4*pi*eps0)) * q1*q2 / d^2 * unit_vector = (const) * (q1*q2/d^3)*R
    %     F_Coulomb = F_Coulomb + (1/(4*pi*epsilon_0)) * (particle_charge * Charge_SC / d^3) * R;
    % end
    Particle_Columb_Position_SC_vector = r' - Columb_Position_SC_array;
    distance_Particle_Columb_Position_SC = vecnorm(Particle_Columb_Position_SC_vector,2,2);
    F_Coulomb_SC_scalar = (1/(4*pi* obj.epsilon_0)) * (obj.particle_charge * obj.Charge_SC) * (1./distance_Particle_Columb_Position_SC).^3;


    % Dynamic check
    if flag_dynamic_case == 1 % (obj.distance_Sun_L1 - r(1) <= obj.first_dynamic_step)
        % Keep Base on

    elseif flag_dynamic_case == 2 % (obj.distance_Sun_L1 - r(1) <= obj.first_dynamic_step + obj.second_dynamic_step)
        % switch off base proportionately
        % this_gain = 1 - ((obj.distance_Sun_L1 - r(1) - obj.first_dynamic_step)/obj.second_dynamic_step);
        F_Coulomb_SC_scalar(1:obj.num_Columb_SC_base) = this_gain*F_Coulomb_SC_scalar(1:obj.num_Columb_SC_base);

    else % (obj.distance_Sun_L1 - r(1) > 3*obj.threshold_particle_base_distance)
        % switch off base completely
        F_Coulomb_SC_scalar(1:obj.num_Columb_SC_base) = 0;
    end

    F_Coulomb_SC = F_Coulomb_SC_scalar .* Particle_Columb_Position_SC_vector;
    F_Coulomb = (sum(F_Coulomb_SC))';


    % Mirror: sum contributions from all satellites in Mirror_Position_SC_array
    F_Mirror = [0;0;0];

    if ~isempty(Mirrior_Position_SC_array)
        Particle_Mirror_Position_SC_unit_vector = r' - Mirrior_Position_SC_array;
        distance_Particle_Mirror_Position_SC = vecnorm(Particle_Mirror_Position_SC_unit_vector,2,2);
        Particle_Mirror_Position_SC_unit_vector = Particle_Mirror_Position_SC_unit_vector./vecnorm(Particle_Mirror_Position_SC_unit_vector,2,2);
        Particle_Mirror_Angles = acosd(dot(Particle_Mirror_Position_SC_unit_vector', Mirrior_Position_SC_orientation_unit_vector')'); % [deg]

        Mirror_SC_idx = find(Particle_Mirror_Angles <= obj.Mirror_cone_half_apex_angle);

        for i_m = 1:1:length(Mirror_SC_idx)

            this_idx = Mirror_SC_idx(i_m);
            this_SC_total_power = (obj.solar_luminosity/ (norm(Mirrior_Position_SC_array(this_idx,:)))^2 ) * obj.Mirror_SC_area; % [W]
            this_SC_flux_at_particle = this_SC_total_power/(obj.Mirror_cone_solid_angle * (distance_Particle_Mirror_Position_SC(this_idx))^2 ); % [W/m^2]

            this_F_Mirror_mag = this_SC_flux_at_particle * (obj.cross_sectional_area * obj.Q_rad_pressure_coefficient) / (4*pi* obj.speed_light);

            F_Mirror = F_Mirror + (this_F_Mirror_mag * Particle_Mirror_Position_SC_unit_vector(this_idx,:))';

        end
    end

    % Total acceleration
    a_total = (F_grav_Sun + F_grav_Earth + F_centripetal + F_SRP_vec + F_Coulomb + F_Mirror) / obj.mass_particle;

    % Derivative
    this_dYdt = [v; a_total];

    dYdt((p-1)*6+1:(p-1)*6+6,1) = this_dYdt;

end


% % Periodic debug print
% if t - lastPrintTime > 50000
%     fprintf("\nTime = %.2f s\n", t);
%     % print magnitudes for particle 1 (as representative)
%     fprintf("  |F_sun grav|     = %.3e N\n", norm(F_grav_Sun));
%     fprintf("  |F_earth grav|   = %.3e N\n", norm(F_grav_Earth));
%     fprintf("  |F_SRP_vec|      = %.3e N\n", norm(F_SRP_vec));
%     fprintf("  |F_Coulomb|      = %.3e N\n", norm(F_Coulomb));
%     fprintf("  |F_Mirror|       = %.3e N\n", norm(F_Mirror));
%     fprintf("  |F_centripetal|  = %.3e N\n", norm(F_centripetal));
%     fprintf("  |a_total|        = %.3e m/s^2\n", norm(a_total));
%     lastPrintTime = t;
% end

end