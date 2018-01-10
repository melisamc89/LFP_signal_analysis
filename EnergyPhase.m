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

%% Mean conmodulogram 

lowf1=1;
lowf2=6;
highf1=20;
highf2=120;
fstep=2;
filter_width=[2 1];
n_fs=400;

protocol={'A','B','C'};
[n m]=size(area);
for a=1:m
    for prot=1:length(protocol)
        count=1;
        clear matrix
        for i=1:length(lfp{a})
            session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),protocol{prot});
            if length(session_number)
                for se=1:length(session_number)
                    order=Wished_Register_Order(area{a},lfp{a}(i,:));
                    register=Rat_Register(files,area{a},rat,order);
                    [eeg,fs]=ReadEEG(register,session_number(se));
                    [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
                    if qual>0.5
                        n_fs=floor(fs/floor(fs/n_fs));
                        eeg=resample(eeg,1,floor(fs/n_fs)); 
                        pos=Loading_Pos(register,session_number(se));
                        tf = isfield(pos.data, 'log');
                        if ~(prot=='C') tf=1; end
                        if tf==1
                            figure(i)
                            amplitude_frequencies = 6:1:120;
                            window_lenght = .5;
                            [ENERGY_PHASE, BIN_CENTERS]=COMODOenergyphase(eeg',n_fs,amplitude_frequencies,9,window_lenght);
                            COMODOenergyphase(eeg',n_fs,amplitude_frequencies,8,window_lenght);
                            saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC\',int2str(rat),'_',name(a,:),'_',lfp{a}(i,:),'_phase'),'png')
                            close(i)
                            matrix(count,:,:)=ENERGY_PHASE;
                            count=count+1;
                        end
                    end
                end
            end
        end
        figure(a)
        [n m l]=size(matrix);
        imagesc(reshape(mean(matrix),[m l]));
        nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC\',int2str(rat),'_COMODO_mean_phase',name(a,:));
        save(nombre,'matrix','-mat')
        saveas(a,nombre,'png')
        close(a)
    end
end

%%