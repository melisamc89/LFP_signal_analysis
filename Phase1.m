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

protocol={'A','B','C'};
[n m]=size(area);
a=2;
clear matrix
prot=1;
count=1;
n_fs=400;

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
                figure(i)
                BW2 = 2;
                WINDOWSIZE=1.5;
                SRATE=n_fs;
                filt_eeg= COMODOfilter(eeg',n_fs,3,BW2,'eegfilt');
                [~,peaks] = findpeaks(real(filt_eeg));
                %peaks=downsample(peaks,140); 
                [aaa,bbb] = meshgrid(-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs),peaks);
                ppp = aaa+bbb;
                ppp(ppp<=0 | ppp>length(eeg)) = length(eeg+1);
                vect2 = mean(eeg(ppp));
                plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,vect2,'b','linewidth',2);
                xlim([-WINDOWSIZE +WINDOWSIZE]);
                xlabel('Time (s)');
                saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Phase\Delta\',name(a,:),'_',lfp{a}(i,:)),'png')
                close(i)
                matrix(count,:,:)=vect2;
                count=count+1;
              end
          end
       end
end

nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Phase\',int2str(prot),'_',name(a,:),'_COMODO_mean_speeds');
save(nombre,'matrix','-mat')

