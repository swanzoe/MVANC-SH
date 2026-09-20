function [Dv,Dp,Fx_v,Fx_p] = CreatReferenceSignal(Pv,Pp,Sv,Sp,PriNoise,N,L,K,M,J)
%% ------------------------------------------------------------------------
% CreatReferenceSignal : build primary disturbances and filtered
%                        references from identified paths.
%
% This is a convenience wrapper that convolves the primary noise with the
% estimated primary and secondary paths to produce the disturbance and
% filtered-reference signals used by the adaptive ANC stages.
%
% Inputs:
%   Pv       : primary path from the primary source to each virtual
%              microphone, size (L_path, M).
%   Pp       : primary path from the primary source to each physical
%              microphone, size (L_path, J).
%   Sv       : secondary paths to the virtual microphones,
%              size (L_path, M, K).
%   Sp       : secondary paths to the physical microphones,
%              size (L_path, J, K).
%   PriNoise : primary reference noise, length N.
%   N        : number of simulation samples.
%   L        : FIR length used to pad the filtered references.
%   K        : number of secondary sources.
%   M        : number of virtual microphones.
%   J        : number of physical microphones.
%
% Outputs:
%   Dv       : primary disturbances at the virtual microphones, (N, M).
%   Dp       : primary disturbances at the physical microphones, (N, J).
%   Fx_v     : filtered references at the virtual mics,
%              size (N+L-1, M, K) (padded with L-1 leading zeros).
%   Fx_p     : filtered references at the physical mics,
%              size (N+L-1, J, K) (padded with L-1 leading zeros).
% ------------------------------------------------------------------------

% Primary disturbances at the virtual microphones.
Dv = zeros(N,M)             ;
for i=1:M
    Dv(:,i) = filter(Pv(:,i),1,PriNoise);
end

% Primary disturbances at the physical microphones.
Dp = zeros(N,J)             ;
for i=1:J
    Dp(:,i) = filter(Pp(:,i),1,PriNoise);
end

% Filtered references at the virtual microphones.
Fx_v = zeros(N+L-1,M,K)     ;
for i = 1:K
    for j = 1:M
        Fx_v(:,j,i) = [zeros(L-1,1);filter(Sv(:,j,i),1,PriNoise)];
    end
end

% Filtered references at the physical microphones.
Fx_p = zeros(N+L-1,J,K)     ;
for i = 1:K
    for j = 1:J
        Fx_p(:,j,i) = [zeros(L-1,1);filter(Sp(:,j,i),1,PriNoise)];
    end
end

% XH = zeros(N,J)     ;
% for i = 1:J
%     XH(:,i) = filter(H((i-1)*L+1:i*L),1,PriNoise);
% end
end
