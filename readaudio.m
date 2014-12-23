function []= readaudio(filename)

% Read an audio waveform

[d,sr] = wavread(filename);

% Save to mp3 (default settings)

mp3write(d,sr,'test.mp3');
[d2,sr] = mp3read(filename);

% mp3 encoding involves some extra padding at each end; we attempt
% to cut it off at the start, but can't do that at the end, because
% mp3read doesn't know how long the original was.  But we do, so..
% Chop it down to be the same length as the original

d2 = d2(1:length(d),:);

% What is the SNR (distortion)?

ddiff = d - d2;
disp(['SNR is ',num2str(10*log10(sum(d(:).^2)/sum(ddiff(:).^2))),' dB']);

subplot(211)
specgram(d(:,1),1024,sr);
subplot(212)
plot(1:5000,d(10000+(1:5000),1),1:5000,d2(10000+(1:5000)));

% Yes, pretty close
%
% NB: lame followed by mpg123 causes a little attenuation; you
% can get a better match by scaling up the read-back waveform:

ddiff = d - 1.052*d2;
disp(['SNR is ',num2str(10*log10(sum(d(:).^2)/sum(ddiff(:).^2))),' dB']);