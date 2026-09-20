function [Xin,d]=Noise_canbine_field(kgrid,medium,PML_size,source,sensor,slit_mask)
%% ------------------------------------------------------------------------
% Noise_canbine_field : run one k-Wave 2-D forward simulation and return
%                       the recorded sensor signals.
%
% This is the workhorse used throughout the main script: it launches a
% kspaceFirstOrder2D simulation with the given source / sensor masks,
% plots the final pressure field and the input/output time waveforms, and
% returns the source signal and the sensor recordings.
%
% Inputs:
%   kgrid     : kWaveGrid object (spatial and temporal grid).
%   medium    : k-Wave medium structure (sound_speed, density, ...).
%   PML_size  : size of the perfectly matched layer [grid points].
%   source    : k-Wave source structure (p_mask, p, ...).
%   sensor    : k-Wave sensor structure (mask, record, ...).
%   slit_mask : 2-D logical array marking the barrier-with-slit cells.
%
% Outputs:
%   Xin       : source pressure signal transposed (row vector).
%   d         : sensor pressure recordings transposed, cast to double.
% ------------------------------------------------------------------------

% Set the k-Wave simulation options.
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, ...
    'DisplayMask', slit_mask+sensor.mask,'DataCast', 'single'};

% Run the simulation.
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% Plot the final pressure field overlaid on the masks.
figure;
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, ...
    sensor_data.p_final + source.p_mask + sensor.mask+slit_mask, [-1, 1]);
colormap(getColorMap);
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;

% Plot the input and sensor time signals with auto-scaled units.
figure;
[t_sc, scale, prefix] = scaleSI(max(kgrid.t_array(:)));

subplot(2, 1, 1);
plot(kgrid.t_array * scale, source.p, 'k-');
xlabel(['Time [' prefix 's]']);
ylabel('Signal Amplitude');
axis tight;
title('Input Pressure Signal');

subplot(2, 1, 2);
plot(kgrid.t_array * scale, sensor_data.p, 'r-');
xlabel(['Time [' prefix 's]']);
ylabel('Signal Amplitude');
axis tight;
title('Sensor Pressure Signal');

% Return the source and sensor signals (transposed to row vectors).
Xin = source.p'      ;
d   = double(sensor_data.p') ;
end
