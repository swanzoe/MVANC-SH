function Noise_canbine_field3_Map(kgrid,medium,PML_size,source,sensor,slit_mask,Nx,Ny,Index,dx,dy,slit_x_pos,FileName)
%% ------------------------------------------------------------------------
% Noise_canbine_field3_Map : run a k-Wave simulation and render / save the
%                            sound-pressure-level (SPL) map.
%
% This variant of the simulation helper records the full-field RMS
% pressure, converts it to dB (referenced to 20 microPa), plots the SPL
% map and a contour plot over the receiver region, reports the average
% and error-microphone SPL, and saves the result to a .mat file.
%
% Inputs:
%   kgrid      : kWaveGrid object.
%   medium     : k-Wave medium structure.
%   PML_size   : PML size [grid points].
%   source     : k-Wave source structure.
%   sensor     : k-Wave sensor structure (must record 'p_rms').
%   slit_mask  : barrier-with-slit mask (2-D array).
%   Nx, Ny     : grid dimensions [grid points].
%   Index      : mask marking the physical error microphones (painted
%                a fixed colour on the SPL map).
%   dx, dy     : grid spacing [m].
%   slit_x_pos : x-position of the barrier [grid points].
%   FileName   : path of the .mat file in which 'map' and 'RMS' are saved.
% ------------------------------------------------------------------------

% Set the k-Wave simulation options.
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, ...
    'DisplayMask', slit_mask,'DataCast', 'single'};

% Run the simulation.
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Reshape the RMS pressure back to the 2-D grid.
figure;
sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Ny]);

% Compute and report the average SPL over the receiver region and the
% SPL at the two physical error microphones.
a_record          = sensor_data.p_rms(1:slit_x_pos-2,:);
b_a               = sensor_data.p_rms(1:slit_x_pos-21,:);
b_record          = 20*log10(a_record *10^6/20)    ;
E_p_record        = 20*log10(sensor_data.p_rms(slit_x_pos-21, Ny/2-8)*10^6/20);
s_n               = numel(b_a);
s_p_record        = sqrt(sum(b_a.^2,'all')/s_n);
s_p_record        = 20*log10(s_p_record*10^6/20)    ;
charter = ['Average SPL:',num2str(s_p_record),' Erro SPL:',num2str(E_p_record)];
disp(charter);

% Render the full-field SPL map (dB). The barrier and the physical mics
% are painted a fixed value so they are not color-scaled.
RMS               = sensor_data;
sensor_data.p_rms = 20*log10(sensor_data.p_rms*10^6/20);
mx = max(abs(sensor_data.p_rms(:)));
mi = min(abs(sensor_data.p_rms(:)));
sensor_data.p_rms(slit_mask == 1) = 256;
sensor_data.p_rms(Index==1)= 256;
sensor_data.p_rms(Index==1)= 256;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, sensor_data.p_rms, [65, 100]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Average SPL');
hold on ;
scatter(-8*dx*1000,(Nx/4-11)*dy*1000,[],'r');
hold on ;
scatter(8*dx*1000,(Nx/4-11)*dy*1000,[],'r');
hold on

% Overlay the virtual microphone circle.
M_vm = 10;          % number of virtual microphones
R = 0.9;            % target-region radius [m]
theta = 2*pi*(0:M_vm-1)/M_vm;
x_vm = 15 +round(R*cos(theta)/dx);
y_vm = round(R*sin(theta)/dy);

scatter(y_vm*dy*1000,x_vm*dx*1000,[],'go');

% Contour plot of the receiver-region SPL, then save the map.
figure ;
[b_y,b_x]=size(b_record);
y_m   = (1:b_y)*dy;
x_m   = (1:b_x)*dx;
x_m   =  x_m-mean(x_m);
[X,Y] = meshgrid(x_m,y_m);
contour(X,Y,flipud(b_record),'ShowText','on');
grid on ;
hold on ;
scatter(-8*dx,19*dy,[],'r');
hold on ;
scatter(8*dx,19*dy,[],'r');
map = flipud(b_record);
save(FileName,'map','RMS');
end
