function [Fx_v,Fx_p] = Ref_Sig_Gene(Sv,Sp,PriNoise,N,L,K,M,J)
%% ------------------------------------------------------------------------
% Ref_Sig_Gene : generate filtered-reference signals from the secondary
%                path estimates.
%
% Convolves the primary noise with each estimated secondary path and
% pads the result with L-1 leading zeros so that the output can be
% indexed as x'(n+L-1) by the adaptive filters.
%
% Inputs:
%   Sv       : secondary paths to the virtual microphones,
%              size (L_path, M, K).
%   Sp       : secondary paths to the physical microphones,
%              size (L_path, J, K).
%   PriNoise : primary reference noise, length N.
%   N        : number of simulation samples.
%   L        : FIR length of the control filters (padding length is
%              L-1 samples).
%   K        : number of secondary sources.
%   M        : number of virtual microphones.
%   J        : number of physical microphones.
%
% Outputs:
%   Fx_v     : filtered references at the virtual mics,
%              size (N+L-1, M, K).
%   Fx_p     : filtered references at the physical mics,
%              size (N+L-1, J, K).
% ------------------------------------------------------------------------

% Virtual-microphone filtered references.
Fx_v = zeros(N+L-1,M,K)     ;
for i = 1:K
    for j = 1:M
        Fx_v(:,j,i) = [zeros(1,L-1),filter(Sv(:,j,i),1,PriNoise)];
    end
end

% Physical-microphone filtered references.
Fx_p = zeros(N+L-1,J,K)     ;
for i = 1:K
    for j = 1:J
        Fx_p(:,j,i) = [zeros(1,L-1),filter(Sp(:,j,i),1,PriNoise)];
    end
end
end
