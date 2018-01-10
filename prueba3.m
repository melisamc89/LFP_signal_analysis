%% Before starting, load everything

clear all
close all
files=Data_Listing();

load('1CA1','-mat');
load('1EC1','-mat');
load('1EC2','-mat');
rat=1;
area={CA1,EC1,EC2};
lfp={Days(CA1),Days(EC1),Days(EC2)};
name=['1CA1';'1EC1';'1EC2'];

load('2CA1','-mat');
load('2EC1','-mat');
load('2EC2','-mat');
rat=2;
area={CA1_2,EC1_2,EC2_2};
lfp={Days(CA1_2),Days(EC1_2),Days(EC2_2)};
name=['2CA1';'2EC1';'2EC2'];


load('3CA1','-mat');
load('3EC1','-mat');
load('3EC2','-mat');
rat=3;
area={CA1_3,EC1_3,EC2_3};
lfp={Days(CA1_3),Days(EC1_3),Days(EC2_3)};
name=['3CA1';'3EC1';'3EC2'];

n_fs=400;

%% prueba

a=1;
i=8;
prot='C';
session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
se=1;
order=Wished_Register_Order(area{a},lfp{a}(i,:));
register=Rat_Register(files,area{a},rat,order);
[eeg,fs]=ReadEEG(register,session_number(se));
n_fs=floor(fs/floor(fs/n_fs));
eeg=resample(eeg,1,floor(fs/n_fs)); 

time=0:1/n_fs:length(eeg)/n_fs-1/n_fs;

amplitude_frequency = 9;   % Amplitude-giving frequency in Hz
phase_frequency = 3;  
filter_width=[3 2];

for i=[6,12,24,50] %100,200,500,1000]
    nf=round(0.75*n_fs*i);
    como=comodulogram(eeg',amplitude_frequency,phase_frequency,n_fs,'WindowType','Hanning','WindowOverlap',[n_fs*i round(0.5*n_fs*i)],'FilterWidth',filter_width);
    t=0:1/nf:time(end)/nf-1/nf;
    t=t(1:length(como));
    plot(t,como)
    hold on
end

%% In time average analysis

[n m]=size(area);
prot='C';
totaltime=20;
non_mov_time=8;
amplitude_frequency = 9;   % Amplitude-giving frequency in Hz
phase_frequency = 3;  
filter_width=[3 1.5];
window_size=2*n_fs;
time=0:1:2*(totaltime+non_mov_time);


for a=1:m
    count=1;
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs)); 
                pos=Loading_Pos(register,session_number(se));
                tf = isfield(pos.data, 'log');
                if ~(prot=='C') tf=1; end
                if tf==1
                    [sr sl]=ReadStartTime(pos);
                    s1=round((sr-non_mov_time)*n_fs);
                    s2=round((sl-non_mov_time)*n_fs);                    
                    er=round((non_mov_time+totaltime)*n_fs)+s1;
                    el=round((non_mov_time+totaltime)*n_fs)+s2;
                    count=1;
                    for trial=1:length(s1)-1
                        como1=comodulogram(eeg(s1:er)',amplitude_frequency,phase_frequency,n_fs,'WindowType','Hanning','WindowOverlap',[window_size window_size-10],'FilterWidth',filter_width);
                        como2=comodulogram(eeg(s2:el)',amplitude_frequency,phase_frequency,n_fs,'WindowType','Hanning','WindowOverlap',[window_size window_size-10],'FilterWidth',filter_width);
                        data_2(count,:)=[como1;como2];
                        count=count+1;
                    end
                end
            end
        end
    end
    time=0:1:2*(totaltime+non_mov_time);

    
end


%% LEaring in coupling throw days


for a=1:m
    count=1;
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs)); 
                pos=Loading_Pos(register,session_number(se));
                tf = isfield(pos.data, 'log');
                if ~(prot=='C') tf=1; end
                if tf==1
                   como1=comodulogram(eeg',9,3,n_fs,'FilterWidth',[3 1.5]);
                   data5(count)=como1;
                   count=count+1;
                end
            end
        end
    end
end