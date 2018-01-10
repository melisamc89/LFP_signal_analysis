%% Before starting, load everything

clear all
close all
files=Data_Listing();

lowf1=1;
lowf2=6;
highf1=5;
highf2=20;
fstep=0.5;
filter_width=[0.5 0.25];
n_fs=400;
sf=50;
n_osc=3;

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
clear matrix1 matrix2
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
                 pos=Loading_Pos(register,session_number(se));
                 vel=pos.data.v;
                 time=floor(pos.data.t*n_fs)+1;
                 counter=zeros(2,1);
                 episodes1=1;episodes2=1;
                 for k=1:length(vel)
                     if vel(k)<10
                        if counter(1)==0
                            counter(1)=1;
                            if counter(2)>n_osc*sf;
                                 st_2(episodes2)=start;
                                 et_2(episodes2)=time(k);
                                 episodes2=episodes2+1;
                            end
                            counter(2)=0;
                            start=time(k);
                        else
                            counter(1)=counter(1)+1;
                        end
                     else
                        if counter(2)==0
                             counter(2)=1;
                             if counter(1)>n_osc*sf
                                 st_1(episodes1)=start;
                                 et_1(episodes1)=time(k);
                                 episodes1=episodes1+1;
                             end
                                 counter(1)=0;
                                 start=time(k);
                        else
                                 counter(2)=counter(2)+1;
                        end
                  end
             end       
             Y1=[];Y2=[];
             for trial=1:length(st_1)
                 Y1=[Y1,eeg(st_1(trial):et_1(trial))'];% slow
             end
             for trial=1:length(st_2)
                 Y2=[Y2,eeg(st_2(trial):et_2(trial))'];% fast
             end
             figure(i)
             subplot(1,2,1)
             como1=comodulogram(Y1,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
             COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como1);
             subplot(1,2,2)
             como2=comodulogram(Y2,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
             COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como2);
             como=[como1,como2];
             color=[max(max(como))];
             for j=1:2
                  subplot(1,2,j)
                  caxis([0 color])
             end
             saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC_DELTA\Open\',name(a,:),'_',lfp{a}(i,:)),'png')
             close(i)
             matrix1(count,:,:)=como1;
             matrix2(count,:,:)=como2;
             count=count+1;
          end
       end
   end
end
figure(a)
[n m l]=size(matrix1);
subplot(1,2,1)
COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix1),[m l]));
subplot(1,2,2)
COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix2),[m l]));
nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC_DELTA\Open\',int2str(prot),'_',name(a,:),'_COMODO_mean_speeds');
saveas(a,nombre,'png')
save(nombre,'matrix1','matrix2','-mat')

