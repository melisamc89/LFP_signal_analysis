BW = 1;
BW2 = 2;
WINDOWSIZE=0.2;
SRATE=n_fs;

filt_eeg= COMODOfilter(eeg',n_fs,75,BW2,'eegfilt');
filt_eeg= COMODOfilter(eeg',n_fs,8,BW2,'eegfilt');


[~,peaks] = findpeaks(real(filt_eeg));
peaks=downsample(peaks,37); 
[aaa,bbb] = meshgrid(-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs),peaks);
ppp = aaa+bbb;
ppp(ppp<=0 | ppp>length(eeg)) = length(eeg+1);
vect2 = mean(eeg(ppp));

plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,vect2,'b','linewidth',2);
xlim([-WINDOWSIZE +WINDOWSIZE]);
xlabel('Time (s)');