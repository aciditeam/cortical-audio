% Include the toolboxes
addpath(genpath('.'));
% Load the cortical toolbox
loadload
% Obtain full path informations
[path, file, ext] = fileparts(audioFile);
% Read the corresponding signal
[sig, fs] = audioread(audioFile);
% First perform audio descriptor extraction
descriptors = processDescriptor(audioFile);
% Parameters to compute the FFT
resampleTo = 22050;
% First resample the signal (similar to ircamDescriptor)
if (fs ~= resampleTo)
    % Be a bit clever and limit the amount of ops
    gcdVal = gcd(fs, resampleTo);
    % Resample the signal
    sig = resample(sig, resampleTo / gcdVal, fs / gcdVal);
    fs = resampleTo;
end
winSize = round(0.060000 * fs);
hopSize = round(0.010000 * fs);
% Finally compute the FFT
fftMat = spectrogram(sig, blackman(winSize), (winSize - hopSize));
% Check that the number of windows are consistent
if (size(fftMat, 2) ~= size(descriptors.SpectralCentroid.value, 1))
    error('Incoherent number of windows');
end
% Final number of time points
nWin = size(fftMat, 2);
% Process STRF (cochlear filter output) representation of sound
frmlen = 16;
% tc : time constant, typically, 4, 16, or 64 ms, etc.
tc = 16;
% fac : nonlinear factor = 0 (full compression), -1 (half-wave rectifier), -2 (linear function)
fac = -1;
% shft : shifted by # of octave, e.g., 0 for 16k, -1 for 8k,
shft = 0;
% filt : filter type ('p' = Powen's IIR, 'p_o' = steeper group delay)
filt = 'p';
% Call to the actual auditory spectrogram transform
y1 = wav2aud(sig, [frmlen, tc, fac, shft], filt, 0);
y1 = y1';
% Find common divisor with number of windows
gcdVal = gcd(size(y1, 2), nWin);
% Final matrix of cochlear data
finalY1 = zeros(size(y1, 1), nWin);
% Resample to the number of time points
for t = 1:size(y1, 1)
    fprintf('Resampling cochlear bank %d.\n', t);
    finalY1(t, :) = resample(y1(t, :), nWin / gcdVal, size(y1, 2) / gcdVal);
end
% Re-organize descriptors as matrix
descNames = fieldnames(descriptors);
finalDescNames = {};
% Find fitting descriptors
for d = 1:length(descNames)
    if size(descriptors.(descNames{d}).value, 1) == nWin
        finalDescNames = [finalDescNames descNames{d}];
    end
end
descMatrix = zeros(length(finalDescNames), nWin);
for d = 1:length(finalDescNames)
    descMatrix(d, :) = descriptors.(finalDescNames{d}).value;
end
fID_desc = fopen([path '/' file '.descriptors.txt'], 'w');
fID_fft = fopen([path '/' file '.fft.txt'], 'w');
fID_aud = fopen([path '/' file '.auditory.txt'], 'w');
% Now finally output all representations
for t = 1:size(fftMat, 2)
    if (sum(abs(fftMat(:, t))) < 0.1)
        fprintf('Skipping window %d (Empty FFT)\n', t);
        continue;
    end
    % First output all descriptors
    for d = 1:size(descMatrix, 1)
        fprintf(fID_desc, '%f ', descMatrix(d, t));
    end
    fprintf(fID_desc, '\n');
    % Then output FFT representation
    for d = 1:size(fftMat, 1)
        fprintf(fID_fft, '%f ', fftMat(d, t));
    end
    fprintf(fID_fft, '\n');
    % First output all descriptors
    for d = 1:size(finalY1, 1)
        fprintf(fID_aud, '%f ', finalY1(d, t));
    end
    fprintf(fID_aud, '\n');
end
fclose(fID_desc);
fclose(fID_fft);
fclose(fID_aud);