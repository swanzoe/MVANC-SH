function [H,Er]=AuxiliaryLMS(LenFilter,NumSource,NumPM,N,Fx_p,DisturPhysic,PriNoise,W,StepSize)
%% ------------------------------------------------------------------------
% AuxiliaryLMS : auxiliary (virtual-sensing) LMS stage.
%
% After the modal control filters W have been tuned on the virtual
% microphone circle, this stage fits an auxiliary filter bank H that
% predicts the residual at each physical error microphone from the
% primary noise. It bridges the virtual-microphone controller (which
% drives W) to the physical microphones used during online control.
%
% Update equations (sample-by-sample):
%   Ep(i)  = DisturPhysic(i) - FX(i) * W        (residual predicted at
%                                                the physical mics)
%   YH(i)  = X(i) * H                            (auxiliary output)
%   Eh(i)  = Ep(i) - YH(i)                       (auxiliary error)
%   H      = H + StepSize * Eh(i) * X(i)
%
% Inputs:
%   LenFilter      : L, length of each control / auxiliary filter.
%   NumSource      : K, number of secondary sources.
%   NumPM          : J, number of physical (error) microphones.
%   N              : number of simulation samples.
%   Fx_p           : filtered reference signals at the physical
%                    microphones, size (N+L-1, J, K).
%   DisturPhysic   : primary disturbances at the J physical microphones,
%                    size (N, J).
%   PriNoise       : primary reference noise, length N.
%   W              : modal control filters from Stage 1, size (K*L, 1)
%                    (one branch per secondary source).
%   StepSize       : adaptation step size.
%
% Outputs:
%   H              : auxiliary filters, size (J*L, 1). Each microphone's
%                    filter occupies a block of length L.
%   Er             : auxiliary error signals, size (J, N).
% ------------------------------------------------------------------------

H  = zeros(NumPM*LenFilter,1);   % Auxiliary filters.
SH = zeros(NumPM*LenFilter,1);   % Increment buffer.
FX = zeros(NumPM,NumSource*LenFilter); % Filtered-reference matrix.
YH = zeros(NumPM,1)  ;            % Auxiliary filter output per mic.
Er = zeros(NumPM,N)  ;
X  = [zeros(LenFilter-1,1); PriNoise];  % Reference delay line.

for i = 1:N
    % Build the filtered-reference rows for every physical mic / source.
    for j = 1:NumPM
        for kk = 1:NumSource
            FX(j,(kk-1)*LenFilter+1:kk*LenFilter) = Fx_p(i+LenFilter-1:-1:i,j,kk)';
        end
        % Auxiliary output for microphone j.
        YH(j) =  X(i+LenFilter-1:-1:i)'*H((j-1)*LenFilter+1:j*LenFilter);
    end

    % Residual predicted at the physical mics after applying W, then the
    % auxiliary error.
    Ep = DisturPhysic(i,:)'-FX*W ;
    Eh = Ep - YH ;

    % LMS update of the auxiliary filters.
    for j = 1:NumPM
        SH((j-1)*LenFilter+1:j*LenFilter) = StepSize*Eh(j)*X(i+LenFilter-1:-1:i);
    end
    H       = H + SH;
    Er(:,i) = Eh    ;
end
end
