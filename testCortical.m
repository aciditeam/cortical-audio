% Include the toolboxes
addpath(genpath('.'));
% Audio file test
audioFile = 'soundfiles/avoice00a.wav';
audioFileSignal = audioread(audioFile);

% Load cortical toolbox environment (colormap, filterbank and parameters)
loadload;

% -------------------
% Moving Ripple Generator
% -------------------

% ---
% mvripfft - moving ripple generator using FFT
% Am: modulation depth, 0 < Am < 1, DEFAULT = .9;
Am = .9;
% Rt: rate (Hz), integer preferred, typically, 1 .. 128, DEFAULT = 6;
Rt = 6;
% Om: scale (cyc/oct), any real number, typically, .25 .. 4, DEFAULT = 2; 
Om = 2;
% Ph: (optional) symmetry (Pi) at f0, -1 < Ph < 1, DEFAULT = 0.
Ph = 0;
% T0: duration (sec), DEFAULT = 1.
T0 = 1;
% f0: center freq. (Hz), DEFAULT = 1000.
f0 = 440;
% SF: sample freq. (Hz), must be power of 2, DEFAULT = 16384
SF = 44100;
% BW: excitation band width (oct), DEFAULT = 5.8.
BW = 5.8;
% RO: roll-off (dB/oct), 0 means log-spacing, DEFAULT = 0;
RO = 0;
% df: freq. spacing, in oct (RO=0) or in Hz (RO>0), DEFAULT = 1/16.
df = 1/16;
% ph_c: component phase
ph_c = [];
% Call to generate a moving ripple stimulus
[s, ph_c, fdx] = mvripfft([Am, Rt, Om, Ph], [T0, f0, SF, BW, RO, df], ph_c);
% s is the time waveform
figure; plot(s);

% ---
% multimvripfft       - generalized version of mvripfft (cf. previous)
% [s, profile] = multimvripfft(rippleList, cond);
% rippleList = [Am1, w1, Om1, Ph1
%               Am2, w2, Om2, Ph2
%               ....
%               AmN, wN, OmN, PhN];
% cond = [T0, f0, SF, BW, RO, df]

% -------------------
% Auditory Processing.
% -------------------

% ---
% cochfil - cochlear filter coefficient reader
%	[CF, B, A] = cochfil(n, shft);
%	HH = cochfil(CHAR);
%	n: channel indices, e.g., 60 or 11:100
%	shft: octave shift
%	CF: characteristic frequencies of the selected channels
%	B, A: IIR coefficients of the selected channel (only one)
%	HH: CHAR = 'o', overall response;
%		CHAR = 'n', normalized overall response.
HH = cochfil('n');

% ---
% wav2aud - fast auditory spectrogramm (for band 180 - 7246 Hz)
% frmlen : frame length, typically, 8, 16 or 2^[natural #] ms.
frmlen = 16;
%	tc	: time const., typically, 4, 16, or 64 ms, etc.
%		  if tc == 0, the leaky integration turns to short-term avg.
tc = 16;
%	fac	: nonlinear factor (critical level ratio), typically, .1 for
%		  a unit sequence, e.g., X -- N(0, 1);
%		  The less the value, the more the compression.
%		  fac = 0,  y = (x > 0),   full compression, booleaner.
%		  fac = -1, y = max(x, 0), half-wave rectifier
%		  fac = -2, y = x,         linear function
fac = -1;
% shft : shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
shft = 0;
% filt : filter type ('p' = Powen's IIR, 'p_o' = steeper group delay)
filt = 'p';
% v: verbose mode
v = 1;
% x: signal waveform
x = audioFileSignal;
% Call to the actual auditory spectrogram transform
v5 = wav2aud(x, [frmlen, tc, fac, shft], filt, v);
% v5 provides the auditory spectrogram
figure; imagesc(v5);

%% ---
% wav2y1 - cochlear filter output
% x: the acoustic input.
x = audioFileSignal;
% filt: filter type ('p' = Powen's IIR, 'p_o' = steeper group delay)
filt = 'p';
% ch_sel: channel selection
ch_sel = -1;
% Call the cochlear filter function
y1 = wav2y1(x, ch_sel, filt);
% y1 : the filter bank output, Lx-by-M 

%% ---
% aud2wav - fast inverse auditory spectrum (for band 180 - 7246 Hz)
% v5: auditory spectrogram (N-by-M)
v5 = v5;
% x0: the projected (guessed) acoustic output (input).
x0 = [];
% L_frm: frame length, typically, 16 ms or 2^[natural #] ms.
L_frm = 16;
% tc: time const, typically, 64 ms = 1024 pts for 16 kHz.
tc = 64;
% fac : nonlinear factor (critical level ratio), typically, .01.
fac = -1;
% shft: shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
shft = 0;
% iter: # of iterartions
iter = 10;
% disp: display the new spectrogram (1) or not (0).
disp = 1;
% snd: play the sound (1) or not (0).
snd = 0;
[x0, xmin, err] = aud2wav(v5, x0, [L_frm, tc, fac, shft, iter, disp, snd]);
% x0: resulting waveform
sound(x0, 44100);
% xmin: the sequence with minimum error
% errv: error vector.

% -- Unprocessed function from auditory aspects --
%	aud2rst     - rate-scale matrix directly from auditory spectrogram
%	aud2sf      - auditory to scale-frequency energy
%	halfregu    - half duration regulator 
%	aud_fix     - fix auditory spectrogram (employed to fix the quality of cor2aud)
%	aud_post    - auditory post windowing (remove the mean of auditory spectrum apply window (to remove the edge effect in the cortical representation)
%	stereaus    - stereausis
%	s = stereaus(x, [frmlen, tc, fac, shft]);
%	x	: the acoustic input, L_x-by-2 matrix.
%	STEREAUS computes the auditory spectrogram for an acoustic waveform.
%	This function takes the advantage of IIR filter's fast performance
%	which not only reduces the computaion but also saves memory space


% A1 Cortical Processing.
% ---
% aud2cors - static cortical representation
% y: auditory spectrum (M-by-1)
y = v5(:, (end/2))';
% sv: char. ripple freq's (K-by-1), e.g., sv = 2.^(-2:.5:3)
sv = 2 .^ (-2:.5:3);
% ser: sample ripple frequency, default = 24 ch/oct
srf = 24;
% full: non-truncation factor
full = 0;
% BP: all bandpass filter bank
BP = 0;
z = aud2cors(y, sv, srf, full, BP);
% z: cortical representation (M-by-K)

% ---
% aud2cor - dynamic cortical representation
% y: auditory spectrogram, N-by-M (samples x channels)
y = v5;
% fullT (fullX): fullness of temporal (spectral) margin in [0, 1].
fullT = 0;
fullX = 0;
% BP: pure bandpass indicator
BP = 0;
% rv: rate vector in Hz, e.g., 2.^(1:.5:5).
rv = 2 .^ (1:.5:5);
% sv: scale vector in cyc/oct, e.g., 2.^(-2:.5:3).
sv = 2 .^ (-2:.5:3);
% fname: cortical representation file (fname='tmpxxx': no data saved)
fname = 'corticalTest.txt';
% disp: saturation level for color panels. If disp < 0 max-normalized.
disp = -1;
% Call the cortical representation function
cr = aud2cor(y, [paras FULLT FULLX BP], rv, sv, fname, disp);
% cr: cortical representation (4D, scale-rate(up-down)-time-freq)

% ---
% gen_cort - generate (bandpass) cortical temporal filter transfer function
% fc: characteristic frequency
fc = 440;
% L: length of the filter, power of 2 is preferable.
L = 2048;
% stf: sample rate.
stf = 22100;
% pass: (vector) [idx K] idx = 1, lowpass; 1<idx<K, bandpass; idx = K, highpass.
pass = [2 3];
%	GEN_CORT generate (bandpass) cortical temporal filter for various
%	length and sampling rate. The primary purpose is to generate 2, 4,
%	8, 16, 32 Hz bandpass filter at sample rate ranges from, roughly
%	speaking, 65 -- 1000 Hz. Note the filter is complex and non-causal.
h = gen_cort(fc, L, stf, pass);

% ---
% gen_corf - generate cortical spectral filter
% fc: characteristic frequency
fc = 440;
% L: length of the filter, power of 2 is preferable.
L = 2048;
% srf: sample rate.
srf = 22100;
% kind: (1 = Gabor function; 2 = Gaussian Function)
kind = 2;
%	GEN_CORF generate (bandpass) cortical filter for various length and
%	sampling rate. The primary purpose is to generate 2, 4, 8, 16, 32 Hz
%	bandpass filter at sample rate 1000, 500, 250, 125 Hz. This can also
%	be used to generate bandpass spatial filter .25, .5, 1, 2, 4 cyc/oct
%	at sample ripple 20 or 24 ch/oct.
h = gen_corf(fc, L, srf, KIND);

%	cor2rst		  - rate-scale matrix along time axis

%	cor2auds	  - inverse static cortical representation
%	cor2aud		  - inverse dynamic cortical representation
%	scalspec	  - scale spectrum

%% ---------------------------
% Generating STRF functions
% ---------------------------
% SCHEMATC  plot the schematic
%	xh = schematc(x);
%	xh = schematc(x, octshft, FCORNAME);
%	x: input at 8 kHz
%
%	This function present an overview of the NSLtools
%	for 1 second, 8000 Hz sample rate, 8 ms frame
%	rv = 2.^(1:5); sv = 2.^(-2:3);
% Auther: Powen Ru (powen@isr.umd.edu), NSL, UMD
% v1.00: 01-Aug-97
% v1.01: 23-Sep-97, added octshft;
octshft = -1;
% time, frequency
L = length(x);
paras(4) = octshft;							% octave shift
paras(1:2) = paras(1:2) / 2 ^(octshft+1);	% frame size and step
paras(3) = -2;								% linear
SF = 16000 * 2^octshft;						% sample rate
T = L / SF * 1000;
% Rate-scale types
hilo_r = {'Slow';'Fast'};
hilo_s = {'Coarse';'Fine'};
nepo = '-+';
% Parameters
sgn_sel = [-1 1];
rv = [4 8 16 32] * 2^(octshft+1);
sv = [.25 .5 1 2 4 8];
r_sel = [2 4];	% 8, 32 Hz
s_sel = [3 5];  % 1, 4 cyc/oct
r_sel1 = [2 4];  % 8, 32 Hz
s_sel1 = [3 5];  % 1, 4 cyc/oct

WD = .2; 
HG = .075;

trng = 192;	 % ms
frng = 1.0;	 % oct
tres = 8;	 % ms
tlen = 64;	 % pts
ch_oct = 24; % ch/oct
flen = 32;	 % ch

tdx = fix(trng/tres); tdx = (1:tdx);
tdx0 = tdx * tres;
fdx = fix(frng*ch_oct); fdx = (-fdx:fdx);
fdx0 = fdx/ch_oct; fdx = fdx + flen;
for sgn_idx = 1:2,
for rdx = 1:2,
for sdx = 1:2,
	LF = .4 + sgn_sel(sgn_idx)*((rdx-1)*.22 + .11);
	BT = .35 + (sdx-1) * .175;
	% rate response
	R1 = gen_cort(rv(r_sel1(rdx)), tlen, 1000/tres);
	r1 = ifft(R1, tlen*2);
	r1 = r1(1:tlen);
	% scale response
	R2 = gen_corf(sv(s_sel1(sdx)), flen, ch_oct);
	r2 = ifft(R2, flen*2);
	r2 = r2([(1:flen)+flen 1:flen]);
	% total response
	if sgn_idx == 1,
		r1 = conj(r1);
	end;
	r = r2 * r1.';
	r = r(fdx, tdx);
    % Plot corresponding STRF
    figure; imagesc(real(r));
    % 
	sub_c = axes('position', [LF, BT, WD, HG]);
    % CHECK
    % CHECK
    % If needed we can pick the corresponding response to the STRF
    % CHECK
    % CHECK
	%z = cor_pick(FCORNAME, 3-sgn_idx, r_sel(rdx), s_sel(sdx));
	%image(cplx_col(z, mean(y(:))/4)');  axis xy
	set(gca, 'xtick', [], 'ytick', [], 'box', 'on');
	
end;
end;
end;


% Manipulations.
%	cor_sel		  - cortical channel selection (weighting)
%	cor_mor		  - cortical morphing
%	cor_map		  - full cortical "mapping"
%	cor_maps	  - static cortical "mapping"
%	aud_maps	  - direct static cortical "mapping"
%	aud_dil		  - dilate auditory spectrogram
%	aud_patch	  - mix timbre-pitch from two auditory spectrograms
%	wav_mor		  - morph two sounds
%
% Utilities.
%	cplxmean	  - re-defined mean for complex number
%	cor_pick	  - pick one channel
%	cor_enhs   	  - enhancement for static cortical representation
%	cor_dist	  - static cortical "distance"
%	cor_info	  - display information of .cor file
%	cochfil		  - cochlear filter bank reader
%	auddist2	  - auditory "distance"
%	rs_sugg		  - rate, scale vectors suggestion
%	rst_view	  - view dynamic cortical representation
%	shiftmat	  - shift a matrix	
%
% Graphics / Performance.
%	schematc	  - plot schematic (try this first !)
%	cor_plot	  - plot selected panels
%	cor_plts  	  - plot static cortical representation
%	aud_plot	  - plot auditory spectrogram
%	subplot1	  - custom subplot design
%	cplx_col	  - indices of [a1map] colormap for complex number
%	real_col	  - indices of [a1map] colormap for real number
%	image_c		  - mono-image for complex matrix using arrows
%	image_q		  - mono-image for complex matrix using quivers
%	image_pm	  - mono-image for matrix with positive-negative numbers
%	figname		  - change the name of current figure
%	figsize		  - set current figure size
%	a1fig		  - new figure with A1 colormap and appropriate name
%	isa1map		  - true for the A1 colormap
%	pb, pd		  - computer performance measure
%	time_est	  - estimate elapse time