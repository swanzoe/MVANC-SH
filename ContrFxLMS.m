function [WC,ErPhysic,ErVirt] = ContrFxLMS(LenFilter,NumSource,NumVM,NumPM,N,Fx_p,Fx_v,DisturPhysic,DisturVirt,PriNoise,H,StepSize)
%% ------------------------------------------------------------------------
% ContrFxLMS : online (control-stage) multichannel filtered-x LMS.
%
% This is the third and final stage of the virtual-microphone ANC
% pipeline. It runs the converged controller on the actual primary noise:
%   * the auxiliary filters H (fitted by AuxiliaryLMS) are used to monitor
%     the residual at the physical microphones;
%   * the control filters WC are updated on those physical residuals using
%     the filtered references Fx_p;
%   * the virtual residuals are also reported (via Fx_v) for diagnosis.
%
% Inputs:
%   LenFilter    : L, length of each control filter.
%   NumSource    : K, number of secondary sources.
%   NumVM        : M, number of virtual microphones.
%   NumPM        : J, number of physical (error) microphones.
%   N            : number of simulation samples.
%   Fx_p         : filtered reference at the physical microphones,
%                  size (N+L-1, J, K).
%   Fx_v         : filtered reference at the virtual microphones,
%                  size (N+L-1, M, K).
%   DisturPhysic : primary disturbances at the J physical mics, size (N,J).
%   DisturVirt   : primary disturbances at the M virtual mics, size (N,M).
%   PriNoise     : primary reference noise, length N.
%   H            : auxiliary filters from AuxiliaryLMS, size (J*L, 1).
%   StepSize     : adaptation step size.
%
% Outputs:
%   WC           : online control filters, size (K*L, 1).
%   ErPhysic     : residual errors at the physical microphones, (J, N).
%   ErVirt       : residual errors at the virtual microphones, (M, N).
% ------------------------------------------------------------------------

FX = zeros(NumPM,NumSource*LenFilter); % Filtered-reference matrix (physical).
FV = zeros(NumVM,NumSource*LenFilter); % Filtered-reference matrix (virtual).
WC = zeros(NumSource*LenFilter,1)    ; % Control filter at the control stage.
ErPhysic = zeros(NumPM,N)  ;
ErVirt   = zeros(NumVM,N)  ;
XH       = zeros(N,NumPM)  ; % Auxiliary-filter outputs per physical mic.

% Precompute the auxiliary-filter outputs XH by filtering the primary
% noise with each branch of H.
for i = 1:NumPM
    XH(:,i) = filter(H((i-1)*LenFilter+1:i*LenFilter),1,PriNoise);
end

for i=1:N
    % Build the physical and virtual filtered-reference rows.
    for j = 1:NumPM
        for kk = 1:NumSource
            FX(j,(kk-1)*LenFilter+1:kk*LenFilter) = Fx_p(i+LenFilter-1:-1:i,j,kk)';
        end
    end
    for j = 1:NumVM
        for kk = 1:NumSource
            FV(j,(kk-1)*LenFilter+1:kk*LenFilter) = Fx_v(i+LenFilter-1:-1:i,j,kk)';
        end
    end

    % Residuals at the physical and virtual microphones after applying
    % the current control filter WC.
    Ep      = DisturPhysic(i,:)' + FX*WC ;
    ErVirt(:,i) = DisturVirt(i,:)' + FV*WC ;

    % Physical residual with the auxiliary output subtracted, and the
    % filtered-x LMS update of WC.
    Eh      = Ep-XH(i,:)'      ;
    WC      = WC - (StepSize*Eh'*FX)';
    ErPhysic(:,i) = Ep               ;
end

end
