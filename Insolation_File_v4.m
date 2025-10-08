% Solar Incidence on Earth with and without Dust Cloud

close all, 
clear all, clc

tic

%% Pre-Set Constants

radius_Earth = 6371; % [km]

radius_Sun = 696340; % [km]

mass_Earth = 5.9722e24; % [kg]

mass_Sun = 1.988400e30; % [kg]

omega_Earth = 2*pi/(365.25*86400); % [rad/sec]

G = 6.67430e-11; % [N⋅m2⋅kg−2]

solar_luminosity = 3.846e26; % [W]

speed_light = 2.99792458e8; % [m/s]

%% SPICE Initialization

flag_file_location = 1;

path_to_MuSCAT_Supporting_Files = '../MuSCAT_Supporting_Files/'; 
path_to_MuSCAT_v2 = '../muscat/'; 

addpath(genpath(path_to_MuSCAT_Supporting_Files));
addpath(genpath(path_to_MuSCAT_v2));
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/de440s.bsp']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/naif0012.tls']);
cspice_furnsh([path_to_MuSCAT_Supporting_Files,'SPICE/pck00011.tpc']);

%% Initialize Simulation

time_utc = '2025-12-21T06:00:00'; % this_time_utc

num_rays = 1e6; % this_num_rays

% save(['all_data_',num2str(num_rays),'_',time_utc,'.mat'])
save(replace(['all_data_',num2str(num_rays),'_',time_utc,'.mat'],':','_'))

all_rays_array = zeros(num_rays,14); % [Boolean hits Dust Cloud, Boolean hits Earth, [x1 y1 z1], [x2 y2 z2], [x3 y3 z3], [x4 y4 z4] ]

%% Earth Position and Orientation

% Convert UTC to ET
et_t = cspice_str2et(time_utc);

% Earth position at the defined time
state = cspice_spkezr('EARTH', et_t, 'J2000', 'NONE', 'SUN');
Earth_pos = state(1:3)'; % [km]

distance_Sun_Earth = norm(Earth_pos); % [km]

% Solve for Sun-Earth L1 point

syms x

eqn = (G*mass_Sun/((x^2) * (1e6))) - (G*mass_Earth/( (distance_Sun_Earth - x)^2 * (1e6) )) - ((omega_Earth^2) * x * 1e3) == 0;

S = solve(eqn,x);

distance_Sun_L1 = double(S(1)); % [km]

distance_Earth_L1 = (distance_Sun_Earth - distance_Sun_L1); % [km]

ratio_L1_Earth = distance_Sun_L1/distance_Sun_Earth;

% Earth orientation at the defined time
rotation_matrix_Earth = cspice_pxform('IAU_EARTH','J2000',et_t);

solar_constant_Earth = solar_luminosity/(4*pi*(distance_Sun_Earth^2)* 1e6); % [W/m2]

%% Outer radius at Dust Cloud location

x_location = distance_Sun_L1; % [km] Location of Dust Cloud

temp_x = (distance_Sun_Earth * radius_Earth)/(radius_Sun - radius_Earth); % [km]

radius_location = (radius_Earth / temp_x ) * (temp_x + distance_Sun_Earth - x_location); % [km]

radius_dust_cloud = 2.2089e+03; % [km] From Dust_Cloud_Mass_Marks_Equations.m, line 502

%% Rotation Matrix for Ray

% Aligning the Primary Vector
primary_vector = func_normalize_vec([distance_Sun_Earth 0 0]);
desired_primary_vector = func_normalize_vec(Earth_pos);

% Compute rotation matrix
if 1 - dot(primary_vector, desired_primary_vector) <= 1e-10
    Rotation_matrix_ray = eye(3);
    warning('Rotation_matrix_ray is Identity!')

elseif 1 + dot(primary_vector, desired_primary_vector) <= 1e-10

    third_vector = func_normalize_vec([primary_vector(2:3), primary_vector(1)]);
    v = func_normalize_vec(cross(primary_vector, third_vector));

    Rotation_matrix_ray = -eye(3) + 2*(v')*v;
    warning('Rotation_matrix_ray is Badly computed!')

else
    v = cross(primary_vector, desired_primary_vector);
    Rotation_matrix_ray = eye(3) + skew(v) + skew(v)^2*(1/(1+dot(primary_vector, desired_primary_vector)));

    % This is Rodrigue's formula; it gives a proper rotation matrix
    % Rot_primary = eye(3) + skew(v) + skew(v)^2*(1-dot(primary_vector, desired_primary_vector))/(norm(v)^2);
end



%% Ray Tracing

flag_Sun_ray_generator = 'uniformly random';

for i=1:1:num_rays

    if mod(i, 100) == 0
        disp(i)
    end

    switch flag_Sun_ray_generator

        case 'bands'

            % Select random location on Sun for the start of ray
            x1 = 0; % [km]
            r1 = rand*radius_Sun; % [km]
            th1 = rand*2*pi; % [rad]
            y1 = r1*cos(th1); % [km]
            z1 = r1*sin(th1); % [km]

            % Select random location on SRP-balanced Sun-Earth L1 for ray
            x2 = x_location; % [km]
            r2 = rand*radius_location; % [km]
            th2 = rand*2*pi; % [rad]
            y2 = r2*cos(th2); % [km]
            z2 = r2*sin(th2); % [km]

        case 'uniformly random'

            % Select random location on Sun for the start of ray
            x1 = 0; % [km]

            flag_rand_variables = 0;
            while flag_rand_variables == 0
                yr1 = 2*(rand-0.5);
                zr1 = 2*(rand-0.5);
                if norm([yr1 zr1]) < 1
                    flag_rand_variables = 1;
                end
            end

            y1 = radius_Sun*yr1; % [km]
            z1 = radius_Sun*zr1; % [km]

            % Select random location on SRP-balanced Sun-Earth L1 for ray
            x2 = x_location; % [km]

            flag_rand_variables = 0;
            while flag_rand_variables == 0
                yr2 = 2*(rand-0.5);
                zr2 = 2*(rand-0.5);
                if norm([yr2 zr2]) < 1
                    flag_rand_variables = 1;
                end
            end

            y2 = radius_location*yr2; % [km]
            z2 = radius_location*zr2; % [km]

        otherwise

            error('Ray generator is not defined!')

    end

    % Does the ray pass through Dust Cloud
    if norm([y2 z2]) <= radius_dust_cloud
        % Yes
        all_rays_array(i,1) = 1;
    end

    % Does the ray hit Earth?
    x3 = distance_Sun_Earth; % [km]

    y3 = (y2 - y1)*(x3 - x1)/(x2 - x1) + y1; % [km]
    z3 = (z2 - z1)*(x3 - x1)/(x2 - x1) + z1; % [km]

    all_rays_array(i,3:11) = [x1 y1 z1 x2 y2 z2 x3 y3 z3]; 

    if (norm([y3 z3])/radius_Earth) <= 1
        % Yes
        all_rays_array(i,2) = 1;

        % % Intersection with Earth's Sphere (using Matlab solve)
        % syms t
        % eqn2 = norm(t*[x2 y2 z2] + (1-t)*[x3 y3 z3] - [distance_Sun_Earth 0 0]) == radius_Earth;
        % S2 = solve(eqn2,t);
        %
        % if double(S2(1)) > 0
        %     t = double(S2(1));
        % elseif double(S2(2)) > 0
        %     t = double(S2(2));
        % else
        %     error('t is wrong!')
        % end
       
        % Intersection with Earth's Sphere (using quadratic equation)

        A = (x2 - x3)^2 + (y2 - y3)^2 + (z2 - z3)^2;
        B = 2*(x2 - x3)*(x3 - distance_Sun_Earth) + 2*(y2 - y3)*y3 + 2*(z2 - z3)*z3;
        C = (x3 - distance_Sun_Earth)^2 + y3^2 + z3^2 - radius_Earth^2;

        t = (-B + sqrt(B^2 - 4*A*C))/(2*A);
        if t < 0
            t = (-B - sqrt(B^2 - 4*A*C))/(2*A);
        end

        Sphere_Intersect_point_ray = t*[x2 y2 z2] + (1-t)*[x3 y3 z3];

        all_rays_array(i,12:14) = Sphere_Intersect_point_ray;   

        % % Earth Intersect Point in J2000 frame
        % Sphere_Intersect_point_J2000 = (Rotation_matrix_ray * Sphere_Intersect_point_ray')'; % [km]
        %
        % Sphere_Intersect_point_relative_Earth = Sphere_Intersect_point_J2000 - Earth_pos; % [km]
        %
        % [radius, lon, lat] = cspice_reclat(rotation_matrix_Earth' * Sphere_Intersect_point_relative_Earth'); % [radius, longitude [rad], latitude [rad] ]
        %
        % longitude = rad2deg(lon); % [deg]
        % latitude = rad2deg(lat); % [deg]
        %
        % all_rays_array(i,3) = latitude; % [deg]
        % all_rays_array(i,4) = longitude; % [deg]

    end


end

disp(['Rays hitting Earth = ',num2str(100*sum(all_rays_array(:,2) == 1)/num_rays),' % ']) % Approx 46%

disp(['Rays hitting Earth and Dust Cloud = ',num2str(100*sum( (all_rays_array(:,2) == 1) & (all_rays_array(:,1) == 1) )/num_rays),' % ']) % Approx 17%

%% Find Latitude and Longitude

all_rays_lat_long = zeros(num_rays,2); % [Latitude Longitude]

for i = 1:1:num_rays

    if all_rays_array(i,2) == 1

        Sphere_Intersect_point_ray = all_rays_array(i,12:14);

        % Earth Intersect Point in J2000 frame
        Sphere_Intersect_point_J2000 = (Rotation_matrix_ray * Sphere_Intersect_point_ray')'; % [km]

        Sphere_Intersect_point_relative_Earth = Sphere_Intersect_point_J2000 - Earth_pos; % [km]

        [radius, lon, lat] = cspice_reclat(rotation_matrix_Earth' * Sphere_Intersect_point_relative_Earth'); % [radius, longitude [rad], latitude [rad] ]

        longitude = rad2deg(lon); % [deg]
        latitude = rad2deg(lat); % [deg]

        all_rays_lat_long(i,1) = latitude; % [deg]
        all_rays_lat_long(i,2) = longitude; % [deg]

    end

end

toc

% save(['all_data_',num2str(num_rays),'_',time_utc,'.mat'])
save(replace(['all_data_',num2str(num_rays),'_',time_utc,'.mat'],':','_'))

%% Plot World Map

plot_handle = figure(1);
clf
set(plot_handle, 'Name','Mission Dashboard')
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','outerposition',[0 0 0.5 0.5])
set(plot_handle,'PaperPositionMode','auto');

hold on

idx_Earth = logical( (all_rays_array(:,2) == 1) & (all_rays_array(:,1) == 0) );
plot(all_rays_lat_long(idx_Earth,2), all_rays_lat_long(idx_Earth,1), '.r')

idx_Earth = logical( (all_rays_array(:,2) == 1) & (all_rays_array(:,1) == 1) );
plot(all_rays_lat_long(idx_Earth,2), all_rays_lat_long(idx_Earth,1), '.b')

load coastlines.mat
plot(coastlon, coastlat, 'k','LineWidth',2)

grid on
axis equal

xlim([-180 180])
ylim([-90 90])

xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
title(['Date = ',time_utc])
set(gca,'FontSize',15, 'FontName','Times New Roman')

% saveas(plot_handle,['Rays_',num2str(num_rays),'_',time_utc,'.png'])
saveas(plot_handle,replace(['Rays_',num2str(num_rays),'_',time_utc,'.png'],':','_'))

