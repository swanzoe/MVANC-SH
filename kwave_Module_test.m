clc
clear
close all
%% Default simulation parameters
%
% This script is the main entry point for a 2-D modal active noise control
% (ANC) simulation built on the k-Wave acoustics toolbox. 
%
% Pipeline:
%   1. Build the k-Wave grid, barrier-with-slit medium and time array.
%   2. Generate the primary (tuning / control) disturbance signals.
%   3. Record the primary disturbances at the physical and virtual mics.
%   4. Identify the secondary paths from each secondary source to every
%      physical and virtual microphone using normalized LMS.
%   5. Tune the control filters on the virtual-microphone circle
%      (Stage 1: ModuleANC_MV_2D).
%   6. Fit auxiliary filters that map virtual residuals back to the
%      physical microphones (Stage 2: AuxiliaryLMS).
%   7. Run the online control stage (ContrFxLMS) and render the final
%      anti-noise field with k-Wave.
scale    = 1 ;
PML_size = 10;
Nx       = scale * 256-2*PML_size;
Ny       = scale * 128-2*PML_size;
Ly       = 2;
dx       = Ly/Ny; % grid point spacing in the x direction [m].
dy       = dx  ; % grid point spacing in the y direction [m].
c0       = 343 ; % the sound speed [m/s].
rho0     = 1.29; % the density of the air [kg/m^3].
kgrid    = kWaveGrid(Nx, dx, Ny, dy);

% Create a mask describing a rigid barrier that contains a slit.
slit_thickness = scale * 2;             % [grid points]
slit_width     = floor(0.6*Ny/Ly);
slit_x_pos     = Nx - Nx/4+10 ;         % [grid points]
slit_offset    = Ny/2 - slit_width/2 - 1;% [grid points]
slit_mask      = zeros(Nx, Ny);
slit_mask(slit_x_pos:slit_x_pos + slit_thickness, 1:1 + slit_offset) = 1;
slit_mask(slit_x_pos:slit_x_pos + slit_thickness, end - slit_offset:end) = 1;

% Assign the barrier material properties (a scaled acoustic impedance) to
% the mask; the remaining domain keeps the air properties.
barrier_scale      = 20;
medium.sound_speed = c0 * ones(Nx, Ny);
medium.density     = rho0 * ones(Nx, Ny);
medium.sound_speed(slit_mask == 1) = barrier_scale * c0  ;
medium.density(slit_mask == 1)     = barrier_scale * rho0;

% Assign the reference sound speed used by the k-Wave CFL condition.
medium.sound_speed_ref = c0;

% Find the time step at the stability limit.
c_ref = medium.sound_speed_ref;
c_max = barrier_scale * c0;
k_max = max(kgrid.k(:));
dt_limit = 2 / (c_ref * k_max) * asin(c_ref / c_max);

% Create the time array, with the time step just below the stability limit.
dt = 0.95 * dt_limit;  % [s]
t_end = 10e-2;         % [s]
kgrid.setTime(round(t_end / dt) + 1, dt);

% Build the time-varying (chirped) sinusoidal primary source.
% Pri_control / Pri_tuning differ only by the chirp rate k, so that the
% tuning run sweeps frequency faster than the control run.
source_mag  = 30;
source_freq = 500;
k = 500/t_end;
Pri_control = source_mag * sin(2 * pi * source_freq * kgrid.t_array + pi *k/2 * kgrid.t_array.^2);
% Pri_control = filterTimeSeries(kgrid, medium, Pri_control);

Pri_tuning = source_mag * sin(2 * pi * source_freq * kgrid.t_array + pi *k * kgrid.t_array.^2);
% Pri_tuning = filterTimeSeries(kgrid, medium, Pri_tuning);

%% Build the disturbance source and the physical / virtual sensor layout
% Primary noise source (speaker) placed at the bottom edge of the domain.
source.p_mask = zeros(Nx, Ny);
source.p_mask(end, Ny/2)  = 1;
source.p = Pri_tuning;

sensor_mic.mask = zeros(Nx, Ny);% all sensors (physical + virtual)
sensor_pm.mask  = zeros(Nx, Ny);% physical error microphones
sensor_vm.mask  = zeros(Nx, Ny);% virtual microphones

% Place M_vm virtual microphones uniformly on a circle of radius R around
% the centre of the slit.
M_vm = 10;          % number of virtual microphones
R = 0.9;            % target-region radius [m]
theta = 2*pi*(0:M_vm-1)/M_vm;
x_vm = Nx/2+15 +round(R*cos(theta)/dx);
y_vm = Ny/2 + round(R*sin(theta)/dy);
sensor_vm.mask(sub2ind([Nx,Ny],x_vm, y_vm)) = 1;

% Reference / physical microphone positions on the receiver side.
x_ref = slit_x_pos;
y_ref = [Ny/2-10,Ny/2+10];
% sensor_pm.mask(x_ref, y_ref) = 1;% ref microphone (commented out)
x_phy = (slit_x_pos-21)*ones(1,2);
y_phy = [Ny/2-10,Ny/2+10];
sensor_pm.mask(sub2ind([Nx,Ny],x_phy, y_phy)) = 1;% physical error mics

% Combined sensor mask used for plotting.
sensor_mic.mask = sensor_pm.mask | sensor_vm.mask;

figure;
spy(slit_mask,'k',8);hold on;
scatter(y_ref,x_ref,  'bo','filled');hold on;
scatter(y_phy,x_phy, 'yo', 'filled');hold on;
scatter(y_vm,x_vm, 'go');
legend('','Ref Mic','','Phy Mic','Vir Mic')

% Modal-ANC options structure (2-D case only).
SHOpts = struct('Fs', {1/dt}, 'c', {c0}, 'R', {R});

%% Record the tuning-stage disturbances
% Primary source active, ANC off.
sensor_vm.record = {'p', 'p_final'};
[~,Sig_tuning_vm]= Noise_canbine_field(kgrid,medium,PML_size,source,sensor_vm,slit_mask);

sensor_pm.record = {'p', 'p_final'};
[~,Sig_tuning]= Noise_canbine_field(kgrid,medium,PML_size,source,sensor_pm,slit_mask);

%% Identify all secondary paths
% ----------------------- Secondary Speaker 1 ---------------------------
source.p_mask = zeros(Nx, Ny);
source.p_mask(slit_x_pos-2,slit_offset-5)=1;
source.p = Pri_tuning;

% Run the forward simulation and fit a tap-delay-line secondary path with
% a normalized LMS filter (dsp.LMSFilter).
[Xin,d]= Noise_canbine_field(kgrid,medium,PML_size,source,sensor_pm,slit_mask);
M   = 1024*4;
muS = 0.05;
Secon1Matrix = zeros(M,nnz(sensor_pm.mask));
for nn= 1:nnz(sensor_pm.mask)
    PathEstimator = dsp.LMSFilter('Method','Normalized LMS','StepSize', muS, ...
    'Length', M);
    [yS,eS,Secon1Matrix(:,nn)] = PathEstimator(Xin,d(:,nn));
    figure
    subplot(2,1,1);
    plot(eS)      ;
    grid on       ;
    title('Residual error');
    subplot(2,1,2);
    plot(Secon1Matrix(:,nn));
    grid on       ;
    title('Secondary Path S1');
end

% Same source, but now record at the virtual microphones.
[Xin,d]= Noise_canbine_field(kgrid,medium,PML_size,source,sensor_vm,slit_mask);
Secon1Matrix_Vm = zeros(M,nnz(sensor_vm.mask));
for nn= 1:nnz(sensor_vm.mask)
    PathEstimator = dsp.LMSFilter('Method','Normalized LMS','StepSize', muS, ...
    'Length', M);
    [yS,eS,Secon1Matrix_Vm(:,nn)] = PathEstimator(Xin,d(:,nn));
    figure
    subplot(2,1,1);
    plot(eS)      ;
    grid on       ;
    title('Residual error');
    subplot(2,1,2);
    plot(Secon1Matrix_Vm(:,nn));
    grid on       ;
    title('SecondaryVic Path S1');
end


% ----------------------- Secondary Speaker 2 ---------------------------
source.p_mask = zeros(Nx, Ny);
source.p_mask(slit_x_pos-2,end - slit_offset+5)=1;
source.p = Pri_tuning;
[Xin,d] = Noise_canbine_field(kgrid,medium,PML_size,source,sensor_pm,slit_mask);

M = 1024*4 ;
muS = 0.05 ;
Secon2Matrix = zeros(M,nnz(sensor_pm.mask));
for nn= 1:nnz(sensor_pm.mask)
    PathEstimator = dsp.LMSFilter('Method','Normalized LMS','StepSize', muS, ...
    'Length', M);
    [yS,eS,Secon2Matrix(:,nn)] = PathEstimator(Xin,d(:,nn));
    figure
    subplot(2,1,1);
    plot(eS)      ;
    grid on       ;
    title('Residual error');
    subplot(2,1,2);
    plot(Secon2Matrix(:,nn));
    grid on       ;
    title('Secondary Path S2');
end

[Xin,d]= Noise_canbine_field(kgrid,medium,PML_size,source,sensor_vm,slit_mask);
Secon2Matrix_Vm = zeros(M,nnz(sensor_vm.mask));
for nn= 1:nnz(sensor_vm.mask)
    PathEstimator = dsp.LMSFilter('Method','Normalized LMS','StepSize', muS, ...
    'Length', M);
    [yS,eS,Secon2Matrix_Vm(:,nn)] = PathEstimator(Xin,d(:,nn));
    figure
    subplot(2,1,1);
    plot(eS)      ;
    grid on       ;
    title('Residual error');
    subplot(2,1,2);
    plot(Secon2Matrix_Vm(:,nn));
    grid on       ;
    title('SecondaryVic Path S2');
end

%% Multichannel virtual-microphone ANC stage parameters
close all
% Switch back to the actual primary source and the slower "control" chirp.
source.p_mask = zeros(Nx, Ny);
source.p_mask(end, Ny/2)  = 1;
source.p = Pri_control;

[~,Sig_control]= Noise_canbine_field(kgrid,medium,PML_size,source,sensor_pm,slit_mask);
[Xin,Sig_control_vm]= Noise_canbine_field(kgrid,medium,PML_size,source,sensor_vm,slit_mask);

Dv_tuning = Sig_tuning_vm;
Dp_tuning = Sig_tuning;
Sv = cat(3, Secon1Matrix_Vm, Secon2Matrix_Vm);
Sp = cat(3, Secon1Matrix, Secon2Matrix);

nfft = length(Pri_tuning);
L_path = M;
S_num = 2;      % number of secondary sources
V_num = M_vm;   % number of virtual microphones
E_num = 2;      % number of physical (error) microphones

pri_noise = Pri_control;
% Build the filtered-reference signals used during tuning.
[Fx_v_tuning, Fx_p_tuning] = Ref_Sig_Gene(Sv,Sp,Pri_tuning,nfft,L_path,S_num,V_num,E_num);
Dv_control = Sig_control_vm;
Dp_control = Sig_control;
% Build the filtered-reference signals used during online control.
[Fx_v_control, Fx_p_control] = Ref_Sig_Gene(Sv,Sp,Pri_control,nfft,L_path,S_num,V_num,E_num);

%% Create the left and right control filters
close all
L_fil = 1024*2;

% Stage 1: tune the control filters on the virtual microphone circle in
% the circular-harmonic (modal) domain.
u1 = 0.000003; % Step size of the modal FxLMS algorithm.
[W_tuning,Er_tuning,Info] = ModuleANC_MV_2D(L_fil,S_num,V_num,nfft,Fx_v_tuning,Dv_tuning,u1,SHOpts);
% [W_tuning,Er_tuning,~] = ModuleANC_MV(L_fil,S_num,V_num,nfft,Fx_v_tuning,Dv_tuning,u1,SHOpts);
% figure;plot(Info.FrameIdx,10*log10(Info.Cost));
figure;
for i = 1:V_num
    subplot(V_num/2,2,i);
    plot(Dv_tuning(:,i)');hold on
    plot(Er_tuning(i,:));
end
sgtitle('TS1:Error for VM', 'FontSize', 14);

% Stage 2: fit auxiliary filters that translate the virtual residuals back
% to the physical error microphones.
u2 = 0.00000001; % Step size of the auxiliary LMS algorithm.
[H_Vm,Er_Vm_au] = AuxiliaryLMS(L_fil,S_num,E_num,nfft,Fx_p_tuning,Dp_tuning,Pri_tuning',W_tuning,u2);
figure;
for i = 1:E_num
    subplot(2,1,i)
    plot(Er_Vm_au(i,:));
end
sgtitle('TS2:Error for VM ', 'FontSize', 14);

%% Online control stage
u3 = 0.0000000006; % Step size of the online FxLMS algorithm.
[WC_control,Er_Ph_control,Ev_Vm_control] = ContrFxLMS(L_fil,S_num,V_num,E_num,nfft,Fx_p_control,Fx_v_control,Dp_control,Dv_control,pri_noise',H_Vm,u3);
set(groot,'defaultAxesTickLabelInterpreter','latex');
figure;
for i = 1:V_num
    subplot(V_num/2,2,i);
    plot(Dv_control(:,i));hold on
    plot(Ev_Vm_control(i,:));
end
sgtitle('CS:Error for VM', 'FontSize', 14);

figure
fs = 1/kgrid.dt;
subplot(1,2,1)
plot(1:L_fil,W_tuning(1:L_fil),1:L_fil,WC_control(1:L_fil))
grid on
xlabel('Taps','Interpreter','latex')
ylabel('${{\bf{w}}_{11}}$','Interpreter','latex')
legend('Control filter 1 in the tuning stage','Control filter 1 in the control stage','Interpreter','latex')
xlim([0,L_fil])
subplot(1,2,2)
[H, F] = freqz(W_tuning(1:L_fil), 1, L_fil, fs);
plot(F, abs(H));
hold on
[H, F] = freqz(WC_control(1:L_fil), 1, L_fil, fs);
plot(F, abs(H));
grid on
xlim([0,3000])
xlabel('Frequency (Hz)','Interpreter','latex')
ylabel('$|{W_{11}}(f)|$','Interpreter','latex')

%% Render the final anti-noise field with k-Wave
% Split the converged control filter into left / right branches.
left_W = WC_control(1:L_fil);
right_W = WC_control(1+L_fil:end);

% Create the anti-noise: primary noise at the bottom edge, plus the two
% secondary speakers driven by the left/right control filters.
source.p_mask = zeros(Nx, Ny);
source.p_mask(end, Ny/2)   = 1;
source.p_mask(slit_x_pos-2,slit_offset-5)=1;
source.p_mask(slit_x_pos-2,end - slit_offset+5)=1;

y1 = filter(left_W,1,Pri_control);
y2 = filter(right_W,1,Pri_control);
source.p = [y2;Pri_control;y1] ;

% Mark the two physical error microphones on the receiver side.
sensor.mask = zeros(Nx, Ny)       ;
sensor.mask(slit_x_pos-21, Ny/2-8) = 1;
sensor.mask(slit_x_pos-21,Ny/2+8)  = 1;
Index = sensor.mask;

% Record the full-field RMS pressure for plotting.
sensor.mask = ones(Nx, Ny);
sensor.record = {'p_rms','u_final', 'p_final'};
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, ...
    'DisplayMask', slit_mask,'DataCast', 'single'};

% Run the final simulation.
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

figure;
sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Ny]);

% Compute and report the SPL over the receiver region.
a_record          = sensor_data.p_rms(1:slit_x_pos-2,:);
b_a               = sensor_data.p_rms(1:slit_x_pos-21,:);
b_record          = 20*log10(a_record *10^6/20)    ;
E_p_record        = 20*log10(sensor_data.p_rms(slit_x_pos-21, Ny/2-8)*10^6/20);
s_n               = numel(b_a);
s_p_record        = sqrt(sum(b_a.^2,'all')/s_n);
s_p_record        = 20*log10(s_p_record*10^6/20)    ;
charter = ['Average SPL:',num2str(s_p_record),' Erro SPL:',num2str(E_p_record)];
disp(charter);

% Plot the average SPL map (dB). The barrier and the physical microphones
% are painted a fixed value so they are not color-scaled.
RMS               = sensor_data;
sensor_data.p_rms = 20*log10(sensor_data.p_rms*10^6/20);
mx = max(abs(sensor_data.p_rms(:)));
mi = min(abs(sensor_data.p_rms(:)));
sensor_data.p_rms(slit_mask == 1) = 256;
sensor_data.p_rms(Index==1)= 256;
sensor_data.p_rms(Index==1)= 256;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, sensor_data.p_rms, [mi, mx]);
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
