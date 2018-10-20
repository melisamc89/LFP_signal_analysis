clear all
directory='/home/melisa/Escritorio/Melisa/Doctorado/';
%Rat for the experiment
%setting parameters for filters
low1=1.5;
low2=6;
high1=5;
high2=12;
n_fs=400;

%for EC1
load('3EC1','-mat');
areaname='EC3';
area=EC1_3;
list=area;
%Loading the avaible files
files=Data_Listing();
totaltime(3)=40;
rat=3;

%for EC1
load('1EC1','-mat');
areaname='EC1';
area=EC1;
list=area;
%Loading the avaible files
files=Data_Listing();
totaltime(1)=24;
rat=1;

%for EC1
load('2EC1','-mat');
areaname='EC2';
area=EC1_2;
list=area;
%Loading the avaible files
files=Data_Listing();
totaltime(2)=24;
rat=2;

%creating the list of different LFP of the recordings. It is not the same
%as the number of recodings and cells.
lfp=Days(list);
n_fs=400;
counter=1;
for index= 1:length(lfp) %goes for every recorded session
        order = Wished_Register_Order(area,lfp(index,:)); %choose the corresponding day of register
        register=Rat_Register(files,area,rat,order); %load the day register
        cells=SimultaneousRecordings(area,lfp(index,1:7)); %create a list of cells that where measured the same day
        [spikes_times session_number]=SpikesAssembly(files,rat,area,cells,'C'); % load spike times of all cells measured
        %in protocol 'C', measured the same day.
        [ncells nse]=size(session_number); %number of cells and repetitions of protocol C
        if nse
            for se=1:nse
                pos=Loading_Pos(register,session_number(1,se));
                tf = isfield(pos.data, 'log');
                if tf==1
                    [eeg fs]=ReadEEG(register,session_number(1,se));
                    eeg_vector=resample(eeg,1,floor(fs/n_fs));
                    n_fs=fs/floor(fs/n_fs);
                    [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
                    if qual>0.5
                        [sr sl]=ReadStartTime(pos);
                        sr=floor(sr*n_fs);sl=floor(sl*n_fs);
                        s1=floor(sr-7*n_fs); s1(1)=s1(1)+1;
                        s2=floor(sl-7*n_fs); s2(1)=s2(1)+1;
                        er=floor(sr+totaltime(rat)*n_fs);er(1)=er(1)+1;
                        el=floor(sl+totaltime(rat)*n_fs);el(1)=el(1)+1;
                        for trial=1:length(sl)-1
                            Y1=eeg_vector(s1(trial):er(trial));
                            Y2=eeg_vector(s2(trial):el(trial));
                            
                            signal{counter}{trial}=[Y1;Y2];
                        end
                        counter=counter+1;
                    end
                end
            end
        end
end