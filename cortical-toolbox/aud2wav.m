function [x0, xmin, errv] = aud2wav(v5, x0, paras)
% AUD2WAV fast inverse auditory spectrum (for band 180 -7246 Hz)
%	[x0, xmin, errv] = aud2wav(v5, x0, [L_frm, tc, fac, shft, ...
%			iter, DISP, SND]);	
%	v5		: auditory spectrogram (N-by-M)
%	x0		: the projected (guessed) acoustic output (input).
%	xmin	: the sequence with minimum error
%	errv	: error vector.
%
%	COCHBA  = (global) [cochead; cochfil]; 
%		cochead : 1-by-M, M-channel header
%		p  = real(cochead); filter order
%		CF = imag(cochead); characteristic frequency
%	cochfil : (L-1)-by-M, M-channel filterbank
%		B = real(cochfil); MA (Moving Average) coefficients.
%		A = imag(cochfil); AR (AutoRegressive) coefficients.
%
%	PARAS	= [L_frm, tc, fac, shft, iter, DISP, SND];
%	L_frm	: frame length, typically, 16 ms or 2^[natural #] ms.
%	tc		: time const, typically, 64 ms = 1024 pts for 16 kHz.
%			if tc == 0, the leaky integration turns to short-term average
%	fac		: nonlinear factor (critical level ratio), typically, .01.
%			The less the value, the more the compression
%			fac = 0: y = (x > 0), full compression
%			fac = -1, y = max(x, 0), half-wave rectifier
%			fac = -2, y = x, linear function
%	shft	: shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
%			etc. SF = 16K * 2^[shft].
%	iter	: # of iterartions
%	DISP	: display the new spectrogram (1) or not (0).
%	SND		: play the sound (1) or not (0).
%
%	AUD2WAV inverts auditory spectrogram V5 back to acoustic input.
%	The COCHBA (in AUD24.MAT) should have been loaded and set to 
%	global beforehand.
%	See also: WAV2AUD	

% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 20-Jan-98

% Revision: Taishih Chi (tschi@isr.umd.edu), NSL, UMD
% v1.10: 24-Aug-98, bug in channel selection [N, M] = size(v5); --> [N, MM] = size(v5);
% v1.11: 02-Oct-98, fix auditory spectrogram. Thanks to David Klein
% v2.00: 24-Jul-01, modify reconstruction algorithm for nonlinear hair cell

% get filter bank,
%	L: max. # of order + 2;
%	M: no. of channels
global COCHBA;
[L, M] = size(COCHBA);	% p_max = L - 2;

% options
INIT	= 1;		% initial channel (highest freq.)
DIFF	= 0;		% difference matching
DEOVR   = 0;		% de-overlapping

% extract parameters 
shft	= paras(4);				% octave shift
fac	= paras(3);				% nonlinear factor
L_frm	= paras(1) * 2^(4+shft);		% frame length (in points)

if paras(2),
	alph = exp(-1/(paras(2)*2^(4+shft)));	% leaky integ. coeff.
	alp1 = exp(-paras(1)/paras(2));		% for one frame
else,
	alph = 0;				% short-term average
	alp1 = 0;
end;

haircell_tc=0.5;				% hair cell time constant
beta= exp(-1/(haircell_tc*2^(4+shft)));

iter	= paras(5);				% # of iter
DISP	= paras(6);				% image/plot
SND	= paras(7);				% play sound for each iter.

% fix the gcf;
if DISP,
	fig_aud  = gcf;
end;

% fix auditory spectrum
if ~isreal(v5), v5 = aud_fix(v5); end;

% get data, allocate memory for ouput 
[N, MM]	= size(v5);
v5_new	= v5;
v5_mean = mean(v5(:));
v5_sum2 = sum(v5(:).^2);
L_x	= N * L_frm;

% de-overlapping vector
vt = v5;	% target vectors
if DEOVR,
	for n = 2:N,
	   		vt(n, :) = max(vt(n, :)-alp1*vt(n-1, :), 0);
	end;
end;

% initial guess
L_x0 = length(x0);
x0   = x0(:);
if L_x0 == 0,
	x0	= rand(L_x, 1)-.5;	% uniform random sequence
	x0	= unitseq(x0);		% N(0, 1)
elseif L_x0 < L_x,
	x0(L_x) = 0;			% zero-padding	
else,
	x0	= x0(1:L_x);		% truncation
end;

% iteration
xmin = x0;	% initial sequence with minimum-error
emin = Inf;	% initial minimum error
errv = [];	% error vector

for idx = 1:iter,
	pb;
	% normalization (Norm(0, 1))
	if fac == 0,	% default: No
		x0 = unitseq(x0);
	end;

	% projected v5
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% last channel
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if INIT,	% default: No
		p		= real(COCHBA(1, M)); 
		NORM	= imag(COCHBA(1, M));
		B		= real(COCHBA((0:p)+2, M));
		A		= imag(COCHBA((0:p)+2, M));
		y1_h	= filter(B, A, x0); 		% forward filtering
		y2_h	= sigmoid(y1_h, fac);		% nonlinear oper.
		if (fac ~= -2), y2_h = filter(1, [1 -beta], y2_h); end;
		%y_cum	= filter(B, A, flipud(y1_h));	% reverse filtering
		%y_cum	= flipud(y_cum)/NORM;
		y_cum	= 0;
	else,
		y1_h	= 0;
		y2_h	= 0;
		y_cum	= 0;
	end;


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% All other channels
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for ch = (M-1):-1:1,

		p	= real(COCHBA(1, ch));		% order of ARMA filter
		NORM	= imag(COCHBA(1, M));
		B	= real(COCHBA((0:p)+2, ch));	% moving average coeff.
		A	= imag(COCHBA((0:p)+2, ch));	% auto-regressive coeff.
		
		% forwarding
		y1 = filter(B, A, x0);		% filter bank
		y2 = sigmoid(y1, fac);		% nonlinear op.
		if (fac ~= -2), y2 = filter(1, [1 -beta], y2); end;
		y3 = y2 - y2_h;			% difference (L-H)
		y4 = max(y3, 0);		% half-wave rect.
		if alph,
			y5 = filter(1, [1 -alph], y4); % leaky integ. 
			vx = y5(L_frm*(1:N));	% new aud. spec.
		else,
			vx = mean(reshape(y4, L_frm, N))';
		end;
		
		v5_new(:, ch) = vx;

		% matching
		s = ones(N, 1);
		for n = 1:N,
			if DEOVR & alp1 & (n>1),
				vx(n) = max(vx(n)-alp1*vx(n-1), 0);
			end;

			% scaling vector
			if vx(n),
				s(n) = vt(n, ch) ./ vx(n);
			elseif vt(n, ch),
				s(n) = 2;	% double it for now
			else,
				s(n) = 1;
			end;
		end;
		%if fac == 0,
		 	%%s = (s-1)/2 + 1;
			%s = sqrt(s);
		%end;
		%dy = y1 - y1_h;
		%?? hard scaling
		s	= (s * ones(1, L_frm))';
		s	= s(:);

		if (fac == -2)			% linear hair cell
			dy = y3;
			y1 = dy.*s;
		else				% nonlinear hair cell
			ind = find(y3 >= 0);
			y1(ind) = y3(ind).*s(ind);
			maxy1p = max(y1(ind));
			ind = find(y3 < 0);
			y1(ind) = y3(ind)/abs(min(y3(ind)))*maxy1p;
		end	
			
		y1_h = y1;
		y2_h = y2;



		% inverse wavelet transform
		y1	= filter(B, A, flipud(y1));	% reverse filtering
		y_cum	= y_cum + flipud(y1)/NORM;	% accumulation
		
	end;

	% previous performance
	v5_r = v5_new / mean(mean(v5_new)) * v5_mean;	% relative v5
	err = sum(sum((v5_r - v5).^2)) ./ v5_sum2;	% relative error
	err = round(err * 10000) / 100;
	era = sum(sum((v5_new - v5).^2)) ./ v5_sum2;	% absolute error
	era = round(era * 10000) / 100;	
	errv = [errv; err era];

	if err < emin,			% minimum error found
		emin = err; xmin = x0; 
	elseif (err-100) > emin,	% blow up !
		y_cum = unitseq(sign(y_cum)+rand(size(y_cum)));
	end;

	% inverse filtering/normalization
	x0 = y_cum*1.01;		% pseudo-normalization

	% play sound (if necessary)
	if SND,
		if shft ~= -1,
			R1 = interpft(x0, round(length(x0)*2^(shft+1)));
		else,
			R1 = x0;
		end; 
		R1 = R1/max(R1)*.9;
		sound(R1, 8000);
	end;

	% display performance message
	errstr = sprintf(['#%i, Err.: %5.2f%% (rel.); %5.2f%% (abs.);' ...
		' Energy: %3.1e'], idx, err, era, sum(x0.^2));
	disp(errstr);

	% plot the auditory spectrogram (if necessary)
	if DISP,
		%subplot(6, 1, 5);
		figure(fig_aud);
		%image(real_col(v5_new)'); axis xy;
		aud_plot(v5_new, paras);
		title(errstr);
		drawnow;
	end;
	pd;

end;

disp(['Minimum Error: ' num2str(emin) '%']);

