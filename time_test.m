clear all
%% Synthetic stationary signal with PAC
%  Synthesize a stationary signal with PAC with the following parameters
amplitude_frequency = 7;   % Amplitude-giving frequency in Hz
phase_frequency = 3;       % Phase-giving frequency in Hz
modulated_amp = .3;         % Amplitude of the amplitude-giving signal 
                            % relative to the phase-giving signal
nonmodulated_amp = .3;      % Amplitude of the non-modulated portion of the  amplitude-given signal relative to the amplitude of the modulation part
noise_amp = 0.1;           % Amplitude of noise relative to the phase-giving signal
modulation_center = 5.0265; % Phase with preffered modulation
modulation_width = 0.75;    % Width of modulation
regularity = .5;            % Regularity of modulation
duration = 100;              % Signal total duration, in seconds
sampling_rate = 500;       % Sampling rate, in Hz
jitter = [15.0 1.0];      

amplitude_frequencies = 1:2:30;
phase_frequencies = 2:1:10;

duration=100;
noise_amp=0.1;
fakesignal = COMODOfakesignal(duration,sampling_rate,amplitude_frequency,...
            phase_frequency,modulated_amp,nonmodulated_amp,noise_amp,...
            modulation_center,modulation_width,regularity);
comodulogram(fakesignal,amplitude_frequencies,phase_frequencies,sampling_rate,'FilterWidth',[2 1]);
mi=comodulogram(fakesignal,amplitude_frequencies,phase_frequencies,sampling_rate,'FilterWidth',[2 1]);


amplitude_frequency = 9;   % Amplitude-giving frequency in Hz
phase_frequency = 3;       % Phase-giving frequency in Hz
j=1;
interval=[0.1:0.1:1,1.5:0.5:5,6:10];
for i=interval
    %subplot(5,1,j)
    como=comodulogram(fakesignal,amplitude_frequency,phase_frequency,sampling_rate,'WindowType',...
        'Hanning','WindowOverlap',[round(sampling_rate*i) round(sampling_rate*i)-1]);
    MI(j)=mean(como);
    MI_std(j)=std(como);
    t=0:1/sampling_rate:length(como)/sampling_rate-1/sampling_rate;
    %plot(t,como)
    %hold on
    j=j+1;
end

errorbar(interval(1:end-1),MI,MI_std,'r')
hold on
plot(interval(1:end-1),MI,'Linewidth',2)
xlabel('Time Window [s]','Fontsize',15)
ylabel('Modularity index')
legend('stdMI','meanMI')

%% Surrogate test

%there is a problema with the in time analysis and the sorrogate

i=0.5;
amplitude_frequency = 80;   % Amplitude-giving frequency in Hz
phase_frequency = 10;       % Phase-giving frequency in Hz
como=comodulogram(fakesignal,amplitude_frequency,phase_frequency,sampling_rate,'Sorrogate',10,...
     'WindowType','Hanning','WindowOverlap',[round(sampling_rate*i) round(sampling_rate*i)-1],'FilterWidth',[2.5 1]);
 

%% look for a part of coupled activity in the trial

duration=[2,0.5,2];
noise_amp=[1,0.1,1];
signal=[];
for i=1:3
    fakesignal = COMODOfakesignal(duration(i),sampling_rate,amplitude_frequency,...
            phase_frequency,modulated_amp,nonmodulated_amp,noise_amp(i),...
            modulation_center,modulation_width,regularity,jitter);
        signal=[signal,fakesignal];
end                
COMODOinspect(signal,sampling_rate,amplitude_frequency,phase_frequency,18);
comodulogram(signal,amplitude_frequencies,phase_frequencies,sampling_rate);

j=1;
for i=[0.75,1,2]
    subplot(3,1,j)
    j=j+1;
    como=comodulogram(signal,amplitude_frequency,phase_frequency,sampling_rate,'WindowType',...
        'Hanning','WindowOverlap',[sampling_rate*i sampling_rate*i-1],'FilterWidth',[2.5 1]);
    t=0:1/sampling_rate:length(como)/sampling_rate-1/sampling_rate;
    plot(t,como)
    %hold on
end

amplitude_frequency = 80;   % Amplitude-giving frequency in Hz
phase_frequency = 10;  
duration=[2,0.5,2];
noise_amp=[1,0.1,1];
for j=1:100
    signal=[];
    for i=1:3
        fakesignal = COMODOfakesignal(duration(i),sampling_rate,amplitude_frequency,...
                phase_frequency,modulated_amp,nonmodulated_amp,noise_amp(i),...
                modulation_center,modulation_width,regularity,jitter);
            signal=[signal,fakesignal];
    end
    como=comodulogram(signal,amplitude_frequency,phase_frequency,sampling_rate,'WindowType',...
        'Hanning','WindowOverlap',[sampling_rate sampling_rate-1],'FilterWidth',[2.5 1]);
    matrix(j,:)=como;
    como2=comodulogram(signal,40,15,sampling_rate,'WindowType',...
        'Hanning','WindowOverlap',[sampling_rate sampling_rate-1],'FilterWidth',[2.5 1]);
    matrix2(j,:)=como2;
end
t=0:1/sampling_rate:length(como)/sampling_rate-1/sampling_rate;
errorbar(t,mean(matrix),std(matrix))
hold on
errorbar(t,mean(matrix),var(matrix),'r')
plot(t,mean(matrix),'k','Linewidth',2)
axis([0 2.5 0 0.1])

%% Dependance on the overlaping window




%% longer signal with many coupling on it.

amplitude_frequency = 80;   % Amplitude-giving frequency in Hz
phase_frequency = 10;       % Phase-giving frequency in Hz
modulated_amp = .3;         % Amplitude of the amplitude-giving signal 
                            % relative to the phase-giving signal
nonmodulated_amp = .3;      % Amplitude of the non-modulated portion of the  amplitude-given signal relative to the amplitude of the modulation part
noise_amp = 0.1;           % Amplitude of noise relative to the phase-giving signal
modulation_center = 5.0265; % Phase with preffered modulation
modulation_width = 0.75;    % Width of modulation
regularity = .5;            % Regularity of modulation
duration = 100;              % Signal total duration, in seconds
sampling_rate = 500;       % Sampling rate, in Hz
jitter = [15.0 1.0];      

%comodulogram(fakesignal,amplitude_frequencies,phase_frequencies,sampling_rate,'FilterWidth',[2.5 1])

position=round(rand(10,1)*length(fakesignal));

for j=1:15
    amplitude_frequency = 80;   % Amplitude-giving frequency in Hz
    phase_frequency = 10;       % Phase-giving frequency in Hz
    noise_amp = 0.5;           % Amplitude of noise relative to the phase-giving signal
    fakesignal = COMODOfakesignal(duration,sampling_rate,amplitude_frequency,...
                phase_frequency,modulated_amp,nonmodulated_amp,noise_amp,...
                modulation_center,modulation_width,regularity,jitter);
    new_fakesignal=fakesignal;
    amplitude_frequency = 80;   % Amplitude-giving frequency in Hz
    phase_frequency = 10;       % Phase-giving frequency in Hz
    for i=1:length(position)
        period = COMODOfakesignal(2,sampling_rate,amplitude_frequency,...
                phase_frequency,modulated_amp,nonmodulated_amp,0.1,...
                modulation_center,modulation_width,regularity,jitter);
        new_fakesignal=[new_fakesignal(1:position(i)),period,new_fakesignal(position(i)+1:end)];
    end
    como=comodulogram(new_fakesignal,amplitude_frequency,phase_frequency,sampling_rate,'WindowType',...
        'Hanning','WindowOverlap',[2*sampling_rate sampling_rate],'FilterWidth',[2.5 1]);
    matrix3(j,:)=como;
end
time=0:1/sampling_rate:length(new_fakesignal)/sampling_rate-1/sampling_rate;
nf=round(sampling_rate);
t=0:1:length(como)-1;
t=t(1:length(como));
