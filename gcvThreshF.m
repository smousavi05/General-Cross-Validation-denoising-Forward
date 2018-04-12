%**************************************************************************
%   GCV denoising 
%**************************************************************************
%
%   [USAGE] 
%   [dn] = gcvThresh(data,opt)
% 
%   [INPUTS] 
%   data:    a structure including the information of input data
%   opt:     a structure including parameters needed for the CWT and
%            denoising
% 
%   [OUTPUTS]
%   dn:      a structure including denoised signal and associated wavelet 
%   coefficients. 
% 
%   [References] 
%   [1]  Mousavi S. M., and C. A. Langston (2016). Hybrid Seismic Denoising 
%   Using Wavelet Block Thresholding and Higher Order Statistics, 
%   Bulletin of Seismological Society of America, 106 (4), 1380-1393, 
%   DOI:10.1785/0120150345. 
%
%   [2]  Mousavi, S. M., and C. A. Langston (2017). Automatic Noise-Removal/
%   Signal-Removal Based on the General-Cross-Validation Thresholding in 
%   Synchrosqueezed domains, and its application on earthquake data, 
%   Geophysics.82(4), V211-V227 doi: 10.1190/geo2016-0433.1                                                                                 
%
%   [3]  Weyrich, N. & Warhola, G.T. (1995). De-noising using wavelets and cross validation.
%   In "Approximation Theory, Wavelets and Applications", (S.P. Sing Ed.), NATO ASI, 
%   Series C, 454, pp. 523--532.
%
%   [4]  Jansen, M., Malfait, M. & Bultheel, A. (1997). Generalized cross validation for
%   wavelet thresholding. Signal Processing, 56, 33--44.
%
%   [5]  Mousavi, S. M., and C. A. Langston, (2016). Fast and novel microseismic 
%   detection using time-frequency analysis. SEG Technical Program Expanded 
%   Abstracts 2016: pp. 2632-2636. doi: 10.1190/segam2016-13262278.1 
%
%-------------------------------------------------------------------------- 
%   S. Mostafa Mousavi 
%   mmousavi@stanford.edu
%   Last time modified: November, 22, 2017
%---------------------------------------------------------------------------

%%%   Copyright (c) 2015 S. Mostafa Mousavi
%%%   All rights reserved.
%%%   This software is provided for non-commercial research purposes only. 
%%%   No warranty is implied by this distribution. Permission is hereby granted
%%%   without written agreement and without license or royalty fees, to use, 
%%%   copy, modify, and distribute this code and its documentation, provided 
%%%   that the copyright notice in its entirety appear in all copies of this 
%%%   code, and the original source of this code, 

function [dn] = gcvThreshF(data,opt)
lensec = length(data.noisy)*data.dt;
 nw = ceil(lensec./opt.wsiz); % number of moving ; 
 width = opt.wsiz*(1/data.dt);  
 
   ed = 0; denoised =[];
   h = waitbar(0.30,'Denoising...')
   for k = 1:width:length(data.noisy);
       h = waitbar((k)./(length(data.noisy)))

       bg = ed +1; 
       ed = width + ed;
       
       if ed > length(data.noisy);
          dd = data.noisy(bg:length(data.noisy));
       else
          dd = data.noisy(bg:ed); 
       end
        
normC = max(dd);
if normC > 0
dd = (dd)/normC;
end

[wtcof1,as] = cwt_fw(dd,opt.type,opt.nv,data.dt);
[nr nc] = size(wtcof1);
kk = zeros(nr,nc);
wk = wtcof1;
dn.wnoisy = wtcof1;
dn.as = as;

% %%% pre-processing
for i = 1:nr
  
    w = real(wtcof1(i,:));
    Vk = 24/length(wtcof1(i,:)); % eq(12) in [1]
    lk = sqrt(Vk)/sqrt(1-0.9);
    kr = (sum((w-mean(w)).^4)./length(w)) ./ (var(w,1).^2)-3; % calculating Kurtosis eq (11) in [1]
    kr(isnan(kr)) = 0;
    if abs(kr) <= lk*opt.gc;
    wk(i,:) = 0;   
    end     
    
end

wk(isnan(wk))=0; wpre=wk;
dn.wpre = wk;
% dn.xpre = cwt_iw(wk, opt.type, opt.nv);

%%% main thresholding
for q = 1:nr;
    if any(wk(q,:)) 
    i = wk(q,:);
    
% Fibonacci
 F0 = 0;
 F1 = 1;
 F2 = 1;
 num = 1;
 
 m = 100;
 b = max(abs(i));
 a =  b / m; 
 tol = 0.0005;
  
 while (b - a) / F2 > tol
    F0 = F1;
    F1 = F2;
    F2 = F0 + F1;
    num = num + 1;
 end
  
 v = a + (F1/F2) * (b - a);
 fv = gcv(i,v);
 u = a + (F0/F2) * (b - a);
 fu = gcv(i,u);
     
 for k = 1:num - 1
 if (fu > fv)
    a = u;
    u = v;
    fu = fv;
    v = b-u+a;
    fv = gcv(i,v);
    firstmove = 1;
 else
    b = v;
    v = u;
    fv = fu;
    u = a+b-v;
    fu = gcv(i,u);
 end      
    a = b/m;  
  end  

thr = b;

wk(q,:)= SoftThresh(i,thr);
    end
end
dn.wgcv = wk;
ddn = cwt_iw(wk, opt.type, opt.nv); % recunstructing the denoised signal.
if normC > 0
ddn = ddn.*normC; % scaling the denoised signal to the size of original one.
end

denoised = [denoised;ddn'];
   end

dn.xgcv = denoised;
close(h)
% for qq = 1:nr;
%  dn.wnoise(qq,:) = abs(wtcof1(qq,:)) - abs(wk(qq,:));
% end
% dn.xnoise = cwt_iw(dn.wnoise, opt.type, opt.nv);


function  gcv = gcv(j,t)

% gcv:    Applies soft-threshold t to all elements of x and computes the generalized 
%         cross-validation threshold needed by the gcvThreh . 
% Inputs
%   x		 input wavelet coefficients, length= 2^J
%   t     soft threshold
% Outputs
%   gcv	 generalized cross-validation threshold

  N = length(j);
  res = abs(j) - t;
  N0 = sum(res <= 0);
  y = sign(j).*res.*(res >= 0);

  if (N0 == 0)
     gcv = 0;
  else
     gcv = N * (norm(y-j) / N0 )^2 ;
  end

  function x = SoftThresh(y,t)
%  Inputs 
%    y     Noisy Data 
%    t     Threshold
%  Outputs 
%    x     sign(y)(|y|-t)_+
%
	res = (abs(y) - t);
	res = (res + abs(res))/2;
	x   = sign(y).*res;
    

  
    
  function [Wx,as] = cwt_fw(x, type, nv, dt, opt)
% Forward continuous wavelet transform, discretized, as described
% in Mallat, S., Wavelet Tour of Signal Processing 3rd ed.Sec. 4.3.3.
%
% [INPUTS]
%     x: input signal vector.
%  type: wavelet type, string
%    nv: number of voices 
%    dt: sampling period 
%   opt: options structure
%
% [OUTPUTS]
%    Wx: [na x n] size matrix (rows = scales, cols = times)
%    as: na length vector containing the associated scales
%
%---------------------------------------------------------------------------------
%    Modified after a wavelet transform by Eugene Brevdo
%---------------------------------------------------------------------------------
    opt = struct();
    opt.rpadded = 0;

    x = x(:); % Turn into column vector
    n = length(x);

    % Padding the signal 
    N = 2^(1+round(log2(length(x)+eps)));
    n1 = floor((N-n)/2); 
    n2 = n1;if (mod(2*n1+n,2)==1), n2 = n1 + 1; end

    xl = padarray(x(:), n1, 'pre');
    xr = padarray(x(:), n2, 'post');
 
    x = [xl(1:n1); x(:); xr(end-n2+1:end)];

    % Choosing more than this means the wavelet window becomes too short
    noct = log2(N)-1;
    assert(noct > 0 && mod(noct,1) == 0);
    assert(nv>0 && mod(nv,1)==0); 
    assert(dt>0);
    assert(~any(isnan(x)));
    
    na = noct*nv;
    as = 2^(1/nv) .^ (1:1:na);
    
    Wx = zeros(na, N);
    
    x = x(:).';
    xh = fft(x);
    
    % for each octave
    for ai = 1:na
        a = as(ai);
        psih = wfilth(type, N, a, opt);
        xcpsi = ifftshift(ifft(psih .* xh));
        Wx(ai, :) = xcpsi;
    end

    % Shorten W to proper size (remove padding)
    if (~opt.rpadded)
        Wx = Wx(:, n1+1:n1+n);
    end

    % Output a for graphing purposes, scale by dt
    as = as * dt;

function x = cwt_iw(Wx, type, nv)
% The inverse wavelet transform
%
% Implements Eq. (4.67) of Mallat, S., Wavelet Tour of Signal Processing 3rd ed.
%
% Inputs:
%  Wx: wavelet transform of a signal, see help cwt_fw
%  type: wavelet used to take the wavelet transform,
%        see help cwt_fw and help wfiltfn
%  opt: options structure used for forward wavelet transform.
%
% Output:
%  x: the signal, as reconstructed from Wx
%
%---------------------------------------------------------------------------------
%    Modified after a wavelet transform written by Eugene Brevdo
%---------------------------------------------------------------------------------
     opt = struct()
    [na, n] = size(Wx);

    % Padding the signal 
    N = 2^(1+round(log2(n+eps)));
    n1 = floor((N-n)/2); 
    n2 = n1;if (mod(2*n1+n,2)==1), n2 = n1 + 1; end
    Wxp = zeros(na, N);

    Wxp(:, n1+1:n1+n) = Wx;
    Wx = Wxp; clear Wxp;
disp(nv)
    noct = log2(N)-1;
    as = 2^(1/nv) .^ (1:1:na);
    
    assert(mod(noct,1) == 0);
    assert(nv>0 && mod(nv,1)==0); 

    % the admissibility coefficient Cpsi
    switch type
      case 'shannon',
        Cpsi = log(2);
      otherwise
        psihfn = wfiltfn(type, opt);
        Cpsi = quadgk(@(x) (conj(psihfn(x)).*psihfn(x))./x, 0, Inf);
    end
    
    % Normalize
    Cpsi = Cpsi / (4*pi);
         
    x = zeros(1, N);
    for ai=1:na
        a = as(ai);
        Wxa = Wx(ai, :);
        psih = wfilth(type, N, a, opt);

        % Convolution theorem 
        Wxah = fft(Wxa);
        xah = Wxah .* psih;
        xa = ifftshift(ifft(xah));
        x = x + xa/a;
    end

     % Take real part and normalize by log_e(a)/Cpsi
     x = log(2^(1/nv))/Cpsi * real(x);

     % Keep the unpadded part
     x = x(n1+1: n1+n);


function [psih] = wfilth(type, N, a, opt)
% Outputs the FFT of the wavelet of family 'type' with parameters
% in 'opt', of length N at scale a: (psi(-t/a))^.
%
% [Inputs]
%   type: wavelet type 
%   N: number of samples to calculate
%   a: wavelet scale parameter 
%   opt: wavelet options 
%   opt.dt: delta t 
%
% [Outputs]
%   psih: wavelet sampling in frequency domain 
%---------------------------------------------------------------------------------
    opt = struct(); 

    k = 0:(N-1);
    xi = zeros(1, N);

    xi(1:N/2+1) = 2*pi/N*[0:N/2];
    xi(N/2+2:end) = 2*pi/N*[-N/2+1:-1];

    psihfn = wfiltfn(type, opt);
    psih = psihfn(a*xi);

    % Normalizing 
    psih = psih * sqrt(a) / sqrt(2*pi);

    % Center around zero in the time domain
    psih = psih .* (-1).^k;


function [psihfn] = wfiltfn(type, opt)
% Wavelet transform function of the wavelet filter in question,
% fourier domain.
%
% [Input]
%   type: string (see below)
%   opt: options structure, e.g. struct('s',1/6,'mu',2)
%
% [Output]
%   psihfn: mother wavelet function ( mexican hat, morlet, shannon, or hermitian)

% Example:
%  psihfn = wfiltfn('bump', struct('mu',1,'s',.5));
%  plot(psihfn(-5:.01:5));
%---------------------------------------------------------------------------------
    switch type
        case 'mhat', % mexican hat
        if ~isfield(opt,'s'), s = 1; else s = opt.s; end
	psihfn = @(w) -sqrt(8)*s^(5/2)*pi^(1/4)/sqrt(3)*w.^2.*exp(-s^2*w.^2/2);
      case 'morlet',
        % can be used with synsq for large enough s (e.g. >5)
        if ~isfield(opt,'mu'), mu = 2*pi; else mu = opt.mu; end
        cs = (1+exp(-mu^2)-2*exp(-3/4*mu^2)).^(-1/2);
        ks = exp(-1/2*mu^2);
        psihfn = @(w)cs*pi^(-1/4)*(exp(-1/2*(mu-w).^2)-ks*exp(-1/2*w.^2));
      case 'shannon',
        psihfn = @(w)exp(-i*w/2).*(abs(w)>=pi & abs(w)<=2*pi);
      case 'hhat', % hermitian hat
        psihfn = @(w)2/sqrt(5)*pi^(-1/4)*w.*(1+w).*exp(-1/2*w.^2);
%       case 'mostafa',
%         load ss
%         if ~isfield(opt,'mu'), mu = 5; else mu = opt.mu; end
%         if ~isfield(opt,'s'), s = 1; else s = opt.s; end
%         psihfnorig = @(w)(0.0720*w.^8)+(0.2746*w.^7)+(0.2225*w.^6)+(-0.2781*w.^5)+(-0.3884*w.^4)+(0.0735*w.^3)+(-0.3354*w.^2)+(-0.0043*w)+(0.3675);
%         psihfn = @(w) psihfnorig((w-mu)/s);
     otherwise
        error('Unknown wavelet type: %s', type);
    end 

      

 
