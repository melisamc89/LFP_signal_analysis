%% Before starting, load everything

clear all
close all
files=Data_Listing();

lowf1=1;
lowf2=6;
highf1=4.5;
highf2=20;
fstep=0.5;
filter_width=[0.5 0.25];
n_fs=400;
sf=50;
n_osc=3;
wind=3;

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
clear matrix1 matrix2 matrix3
prot=1;
count=1;
clear st_1 st_2 et_1 et_2;
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
                 if length(eeg)>3*wind*60*n_fs
                 Y1=eeg(1:wind*60*n_fs)';% init
                 Y2=eeg(round(length(eeg)/2)-1.5*wind*60*n_fs:round(length(eeg)/2)+1.5*wind*60*n_fs)';% middle
                 Y3=eeg(end-wind*60*n_fs:end)';% last
                 figure(i)
                 subplot(1,3,1)
                 como1=comodulogram(Y1,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                 COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como1);
                 subplot(1,3,2)
                 como2=comodulogram(Y2,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                 COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como2);
                 subplot(1,3,3)
                 como3=comodulogram(Y3,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                 COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                 como=[como1,como2,como3];
                 color=[max(max(como))];
                 for j=1:3
                      subplot(1,3,j)
                      caxis([0 color])
                 end
                 saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC_DELTA\Open\parts\',name(a,:),'_',lfp{a}(i,:)),'png')
                 close(i)
                 matrix1(count,:,:)=como1;
                 matrix2(count,:,:)=como2;
                 matrix3(count,:,:)=como3;
                 count=count+1;
              end
          end
       end
   end
end
figure(a)
[n m l]=size(matrix1);
subplot(1,3,1)
COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix1),[m l]));
subplot(1,3,2)
COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix2),[m l]));
subplot(1,3,3)
COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix3),[m l]));
nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC_DELTA\Open\parts\',int2str(prot),'_',name(a,:),'_COMODO_mean_speeds');
saveas(a,nombre,'png')
save(nombre,'matrix1','matrix2','matrix3','-mat')

