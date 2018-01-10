%% Before starting load everything

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

lowf1=log(-0.01);
lowf2=4;
highf1=log(0.8);
highf2=20;
fstep=0.50;
n_fs=400;
number=20;
freq_data.phasef1=lowf1;
freq_data.phasef2=highf1;
freq_data.ampf1=lowf2;
freq_data.ampf2=highf2;
freq_data.sf=n_fs;
freq_data.ampstep=fstep;
freq_data.phasestep=logspace(-0.01,0.8,number);

protocol={'A','B','C'};
[n m]=size(area);
for a=2:2
    for prot=[1,2,3]
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
                            %como=comodulogram(eeg',lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterWidth',filter_width);
                            como=comodulogram(eeg',lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como);
                            colormap(jet)
                            %saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC\',int2str(prot),'_',name(a,:),'_',lfp{a}(i,:)),'png')
                            caxis('auto')
                            saveas(i,strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/10-05/',int2str(prot),'_',name(a,:),'_',lfp{a}(i,:)),'png')
                            close(i)
                            matrix(count,:,:)=como;
                            count=count+1;
                        end
                    end
                end
            end
        end
        figure(a)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix),[number length(lowf2:fstep:highf2)]));
        colormap(jet)
        %nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC\',int2str(prot),'_COMODO_mean_',name(a,:));
        nombre=strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/10-05/',int2str(prot),'_COMODO_mean_',name(a,:));
        save(nombre,'matrix','freq_data','-mat')
        saveas(a,nombre,'png')
        close(a)
    end
end
