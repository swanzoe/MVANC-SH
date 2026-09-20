function [W,Er,Info] = ModuleANC_MV_2D(LenFilter,NumSource,NumMic,N,Fx,Distur,StepSize,CHOpts)
%% ------------------------------------------------------------------------
% ModuleANC_MV_2D : 2-D circular-harmonic-domain multichannel FxLMS.
%
% This function implements ONLY the first-stage modal ANC update:
%
%   boundary residuals
%        -> short-time spectra
%        -> circular-harmonic coefficients
%        -> regional modal cost
%        -> control-filter update
%
% No virtual-sensing auxiliary-filter stage is included here.
%
% Inputs:
%   LenFilter : L, length of each control filter.
%   NumSource : K, number of secondary sources.
%   NumMic    : M, number of error microphones uniformly distributed on
%               the boundary circle.
%   N         : number of simulation samples.
%   Fx        : filtered reference signals.
%               Size:
%                 (N+L-1, M, K)       for J = 1
%                 (N+L-1, M, K, J)    for J > 1
%               with
%                 Fx(n+L-1,m,k,j) = x'_{jkm}(n).
%   Distur    : primary disturbances at the M boundary microphones,
%               size (N,M).
%   StepSize  : adaptation step size.
%   CHOpts    : optional structure:
%       .Fs          sampling frequency [16000]
%       .c           sound speed [343]
%       .R           radius of target circle [0.30]
%       .Angle       microphone azimuths in rad [uniform on 0,2*pi)
%       .Weights     circular quadrature weights [2*pi/M]
%       .NDFT        spectral-analysis length [2^nextpow2(max(L,64))]
%       .Hop         samples between updates [NDFT]
%       .Band        [Fmin Fmax] in Hz [0 Fs/2]
%       .MaxOrder    maximum circular-harmonic order
%                    [floor((M-1)/2)]
%       .Window      'rect' | 'hann' ['rect']
%       .Weighting   'none' | 'radial' ['none']
%       .Epsilon     Bessel-denominator regularization [1e-6]
%       .GradScale   additional numerical gradient scaling [1]
%       .Norm        'none' | 'power' ['none']
%
% Outputs:
%   W     : final control filters, size (K*L*J,1).
%           Block w_{kj} starts at
%           ((j-1)*K+(k-1))*L+1.
%   Er    : actual online residual pressures, size (M,N).
%   Info  : circular-harmonic geometry, retained orders and convergence
%           diagnostics.
%
% -------------------------------------------------------------------------
% 2-D circular-harmonic model
%
% Real orthonormal circular-harmonic basis:
%
%   C_0(theta)     = 1/sqrt(2*pi)
%   C_qc(theta)    = cos(q theta)/sqrt(pi)
%   C_qs(theta)    = sin(q theta)/sqrt(pi), q >= 1
%
% Boundary coefficient:
%
%   a_q(i) ~= sum_m omega_m E_m(i) C_q(theta_m)
%
% Regional cost:
%
%   J = 1/2 sum_i sum_q lambda_q(i) |a_q(i)|^2
%
% where, for radial weighting,
%
%   lambda_q(i) = int_0^R |rho_q(k_i,r)|^2 r dr,
%   rho_q(k_i,r) = J_q(k_i r)/J_q(k_i R).
%
% The function uses an EXPLICIT tap-spectrum gradient.  This is slower than
% the FFT-shift acceleration used in the original 3-D prototype, but it is
% preferable for validating the derivation because it avoids periodic-shift
% assumptions in recovering the FIR-tap gradient.
%% ------------------------------------------------------------------------

if nargin < 8
    CHOpts = struct();
end

L = LenFilter;
K = NumSource;
M = NumMic;

%% ---- 0. Input dimensions ------------------------------------------------
if ndims(Fx) < 4
    J = 1;
    Fx = reshape(Fx,size(Fx,1),M,K,1);
else
    J = size(Fx,4);
end

if size(Distur,1) ~= N || size(Distur,2) ~= M
    error('ModuleANC_MV_2D: Distur must have size (N,NumMic).');
end

if size(Fx,1) < N+L-1 || size(Fx,2) ~= M || size(Fx,3) ~= K
    error('ModuleANC_MV_2D: Fx must have size at least (N+L-1,M,K,J).');
end

%% ---- 1. Options ----------------------------------------------------------
Opt = CHOpts;
Opt = SetDefault(Opt,'Fs',16000);
Opt = SetDefault(Opt,'c',343);
Opt = SetDefault(Opt,'R',0.30);
Opt = SetDefault(Opt,'NDFT',2^nextpow2(max(L,64)));
Opt = SetDefault(Opt,'Hop',Opt.NDFT);  % strict block-gradient default
Opt = SetDefault(Opt,'Band',[0 Opt.Fs/2]);
Opt = SetDefault(Opt,'MaxOrder',floor((M-1)/2));
Opt = SetDefault(Opt,'Window','rect');
Opt = SetDefault(Opt,'Weighting','none');
Opt = SetDefault(Opt,'Epsilon',1e-6);
Opt = SetDefault(Opt,'GradScale',1);
Opt = SetDefault(Opt,'Norm','none');

Ndft = round(Opt.NDFT);
Hop  = max(1,round(Opt.Hop));

if Ndft < 1 || Hop < 1
    error('ModuleANC_MV_2D: NDFT and Hop must be positive.');
end

if Ndft > N
    error('ModuleANC_MV_2D: NDFT must not exceed N.');
end

% Unlike the accelerated 3-D prototype, the explicit tap-gradient below
% does not require NDFT >= L.  Keeping NDFT >= L is nevertheless convenient
% for spectral resolution and is recommended.
if Ndft < L
    warning('ModuleANC_MV_2D: NDFT < LenFilter. This is allowed here, but NDFT >= L is recommended.');
end

%% ---- 2. Circular microphone geometry ------------------------------------
% In 2-D, the microphones lie on r = R and are uniformly distributed in
% azimuth.  This is the circular analogue of equal-area spherical sampling.
if ~isfield(Opt,'Angle') || isempty(Opt.Angle)
    Opt.Angle = 2*pi*(0:M-1)'/M;
end

Ang = Opt.Angle(:);

if numel(Ang) ~= M
    error('ModuleANC_MV_2D: Angle must contain NumMic entries.');
end

% Circular quadrature weights.  The weights are normalized so that
% sum_m omega_m = 2*pi.
if ~isfield(Opt,'Weights') || isempty(Opt.Weights)
    Wq = (2*pi/M)*ones(M,1);
else
    Wq = Opt.Weights(:);
    if numel(Wq) ~= M
        error('ModuleANC_MV_2D: Weights must contain NumMic entries.');
    end
    Wq = 2*pi*Wq/sum(Wq);
end

% For real circular harmonics up to Q, the number of basis functions is
% 2*Q+1.  Hence the sampling condition is 2*Q+1 <= M.
Qcap = max(0,min(round(Opt.MaxOrder),floor((M-1)/2)));

%% ---- 3. Frequency bins and retained order -------------------------------
Freq = (0:floor(Ndft/2))*Opt.Fs/Ndft;

% DC is excluded from the adaptive modal cost.  Nyquist is retained when it
% falls inside the requested band.
Sel  = Freq > 0 & Freq >= Opt.Band(1) & Freq <= Opt.Band(2);
Bins = find(Sel)-1;                       % zero-based DFT bin index

if isempty(Bins)
    error('ModuleANC_MV_2D: the retained frequency band contains no positive-frequency DFT bin.');
end

nB = numel(Bins);

Kw = 2*pi*Freq(Bins+1)/Opt.c;            % acoustic wavenumber k_i
Kw = Kw(:).';

% Practical spatial-bandwidth rule for a disk:
% Q_i ~ ceil(k_i R), capped by microphone sampling.
Ord = min(Qcap,max(0,ceil(Kw*Opt.R)));    % 1 x nB

%% ---- 4. Real circular-harmonic basis ------------------------------------
[Ych,HarmOrder,ModeLabel] = RealCHMatrix(Qcap,Ang);
NumMode = size(Ych,2);                    % = 2*Qcap+1

% Frequency-dependent mode truncation.
Mask = double(HarmOrder <= Ord);          % NumMode x nB

% Optional radial weighting corresponding to the area-integrated acoustic
% energy over the target disk.
if strcmpi(Opt.Weighting,'radial')
    Beta = RadialOrderWeight2D(Qcap,Kw,Opt.R,Opt.Epsilon);
    Bta  = Mask.*Beta;
elseif strcmpi(Opt.Weighting,'none')
    Beta = ones(NumMode,nB);
    Bta  = Mask;
else
    error('ModuleANC_MV_2D: Weighting must be ''none'' or ''radial''.');
end

%% ---- 5. Spectral window --------------------------------------------------
switch lower(Opt.Window)
    case 'hann'
        Win = hann(Ndft,'periodic');
    case 'rect'
        Win = ones(Ndft,1);
    otherwise
        error('ModuleANC_MV_2D: Window must be ''rect'' or ''hann''.');
end

%% ---- 6. Adaptive process -------------------------------------------------
W  = zeros(K*L*J,1);
FX = zeros(M,K*L*J);
Er = zeros(M,N);

NumFrame  = floor((N-Ndft)/Hop)+1;
CostHist  = zeros(max(NumFrame,0),1);
ModalNorm = zeros(max(NumFrame,0),1);
FrameIdx  = zeros(max(NumFrame,0),1);
nFrame    = 0;

for n = 1:N

    % 6.1 Actual online residual at the current sample.
    FX = BuildFXMatrix(Fx,n,L,M,K,J);
    Ev = Distur(n,:).' - FX*W;
    Er(:,n) = Ev;

    % 6.2 Modal-domain update.
    if n >= Ndft && mod(n-Ndft,Hop) == 0

        % Recompute the residual frame using the CURRENT W throughout the
        % whole frame.  This keeps the spectral gradient consistent even
        % when Hop < NDFT and avoids mixing residuals produced by different
        % historical controller values.
        Ewin = zeros(Ndft,M);
        for ell = 0:Ndft-1
            s = n-ell;
            FXs = BuildFXMatrix(Fx,s,L,M,K,J);
            Ewin(ell+1,:) = (Distur(s,:).' - FXs*W).';
        end

        % Normalized short-time spectra E_m(i,n).
        EfFull = fft(Ewin.*Win,Ndft,1)/Ndft;
        Ef     = EfFull(Bins+1,:).';       % M x nB

        % 6.3 Circular-harmonic coefficients.
        A = Ych.'*(Wq.*Ef);                % NumMode x nB
        A = A.*Mask;

        % 6.4 Back-project the weighted conjugate modal residual to the
        % microphone domain:
        %
        % B_m(i) = sum_q lambda_q(i) a_q^*(i) C_q(theta_m).
        B = Ych*(Bta.*conj(A));            % M x nB

        % 6.5 Explicit FIR-tap gradient.
        %
        % For tap l:
        % X'_{jkm,l}(i,n) =
        %   sum_{ell=0}^{NDFT-1} x'_{jkm}(n-ell-l)
        %   exp(-j 2*pi*i*ell/NDFT).
        %
        % The explicit calculation is intentionally used for correctness
        % testing and avoids the circular-shift approximation.
        Grad = zeros(L,K*J);
        Px   = 0;

        for j = 1:J
            for k = 1:K
                col = (j-1)*K+k;

                for l = 0:L-1
                    Xseq = zeros(Ndft,M);

                    for ell = 0:Ndft-1
                        s = n-ell-l;
                        row = s+L-1;       % Fx(row,...) stores x'(s)

                        if row >= 1 && row <= size(Fx,1)
                            Xseq(ell+1,:) = reshape(Fx(row,:,k,j),1,M);
                        end
                    end

                    Xfull = fft(Xseq.*Win,Ndft,1)/Ndft;
                    Xsel  = Xfull(Bins+1,:).';     % M x nB

                    % g_l = sum_i Re{ sum_m omega_m B_m(i) X'_m(i,l) }.
                    Grad(l+1,col) = real(sum(sum((Wq.*B).*Xsel)));

                    if strcmpi(Opt.Norm,'power')
                        Px = Px + sum(sum((Wq.*ones(1,nB)).*abs(Xsel).^2));
                    end
                end
            end
        end

        Mu = StepSize*Opt.GradScale;

        if strcmpi(Opt.Norm,'power')
            Px = Px/(max(1,L*K*J*nB));
            Mu = Mu/(Px+eps);
        elseif ~strcmpi(Opt.Norm,'none')
            error('ModuleANC_MV_2D: Norm must be ''none'' or ''power''.');
        end

        W = W + Mu*reshape(Grad,L*K*J,1);

        % 6.6 Convergence diagnostics.
        nFrame = nFrame+1;
        ModalNorm(nFrame) = sum(abs(A(:)).^2);
        CostHist(nFrame)  = 0.5*sum(sum(Bta.*abs(A).^2));
        FrameIdx(nFrame)  = n;
    end
end

%% ---- 7. Information for post-processing ---------------------------------
Info.Wmat        = reshape(W,L,K,J);
Info.Basis       = Ych;
Info.ModeOrder   = HarmOrder;
Info.ModeLabel   = ModeLabel;
Info.Weights     = Wq;
Info.Angle       = Ang;
Info.Xmic        = Opt.R*cos(Ang);
Info.Ymic        = Opt.R*sin(Ang);
Info.Bins        = Bins;
Info.Freq        = Freq(Bins+1);
Info.Wavenum     = Kw;
Info.Order       = Ord;
Info.MaxOrder    = Qcap;
Info.Mask        = Mask;
Info.RadialWeight = Beta;
Info.OrderWeight  = Bta;
Info.ModalNorm   = ModalNorm(1:nFrame);
Info.Cost        = CostHist(1:nFrame);
Info.FrameIdx    = FrameIdx(1:nFrame);
Info.Options     = Opt;

% For a sufficiently sampled uniform circular array, this should be close
% to zero for the retained basis.
Info.OrthError = norm(Ych.'*(Wq.*Ych)-eye(NumMode),'fro');

end

%% =========================================================================
function Opt = SetDefault(Opt,Field,Value)
if ~isfield(Opt,Field) || isempty(Opt.(Field))
    Opt.(Field) = Value;
end
end

%% =========================================================================
function FX = BuildFXMatrix(Fx,n,L,M,K,J)
% Build the M x (K*L*J) filtered-reference matrix at sample n:
%
% [ x'_{11m}(n) ... x'_{11m}(n-L+1) | ... ].
FX = zeros(M,K*L*J);

for m = 1:M
    for k = 1:K
        for j = 1:J
            Off = ((j-1)*K+(k-1))*L;
            FX(m,Off+1:Off+L) = reshape(Fx(n+L-1:-1:n,m,k,j),1,L);
        end
    end
end
end

%% =========================================================================
function [Y,Order,Label] = RealCHMatrix(Qmax,Ang)
% Real orthonormal circular-harmonic matrix.
%
% Column order:
%   q = 0,
%   cos(theta), sin(theta),
%   cos(2 theta), sin(2 theta), ...
%
% Normalization:
%   integral_0^{2*pi} C_a(theta) C_b(theta) dtheta = delta_ab.

M = numel(Ang);
NumMode = 2*Qmax+1;

Y     = zeros(M,NumMode);
Order = zeros(NumMode,1);
Label = cell(NumMode,1);

Y(:,1) = ones(M,1)/sqrt(2*pi);
Order(1) = 0;
Label{1} = 'q0';

idx = 2;

for q = 1:Qmax
    Y(:,idx) = cos(q*Ang)/sqrt(pi);
    Order(idx) = q;
    Label{idx} = sprintf('cos%d',q);

    Y(:,idx+1) = sin(q*Ang)/sqrt(pi);
    Order(idx+1) = q;
    Label{idx+1} = sprintf('sin%d',q);

    idx = idx+2;
end
end

%% =========================================================================
function Beta = RadialOrderWeight2D(Qmax,Kw,R,Eps)
% Area-energy weighting for a 2-D circular target region.
%
% lambda_q(k_i) = int_0^R |rho_q(k_i,r)|^2 r dr,
%
% rho_q(k_i,r) = J_q(k_i r)/J_q(k_i R).
%
% A regularized magnitude-squared ratio is used:
%
% |rho_q|^2 ~= |J_q(k_i r)|^2 / (|J_q(k_i R)|^2 + Eps).
%
% The integral is normalized by R^2/2 only to keep the numerical step-size
% scale convenient.  Remove the normalization if the absolute physical
% volume/area energy scaling is required.

nB = numel(Kw);
NumMode = 2*Qmax+1;

Nr = 128;
r  = linspace(0,R,Nr).';

BetaOrder = zeros(Qmax+1,nB);

for q = 0:Qmax
    for i = 1:nB
        jr = besselj(q,Kw(i)*r);
        jR = besselj(q,Kw(i)*R);

        rho2 = abs(jr).^2/(abs(jR).^2 + Eps);
        BetaOrder(q+1,i) = trapz(r,rho2.*r)/(R^2/2);
    end
end

% Map order weights to the real basis:
% q=0 -> one mode; q>=1 -> cosine and sine modes share the same weight.
Beta = zeros(NumMode,nB);
Beta(1,:) = BetaOrder(1,:);

idx = 2;
for q = 1:Qmax
    Beta(idx,:)   = BetaOrder(q+1,:);
    Beta(idx+1,:) = BetaOrder(q+1,:);
    idx = idx+2;
end
end
