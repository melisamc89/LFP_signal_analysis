%% Before starting, load everything

clear all
close all
files=Data_Listing();

n_fs=400;
sf=50;
n_osc=1;

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
rat=3;
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

protocol={'A','B','C'};
[n m]=size(area);
a=2;
prot=1;
count=1;
inf_limit=[0,5,10,15,20,25,30,35];
sup_limit=[5,10,15,20,25,30,35,40];
Y={0,0,0,0,0,0,0,0};
for i=1:length(lfp{a})
    st={0,0,0,0,0,0,0,0,0};
    et={0,0,0,0,0,0,0,0,0};
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
                 counter=zeros(8,1);
                 episodes=ones(8,1);
                 for k=1:length(vel)
                     for limit=1:8
                         if vel(k)>inf_limit(limit) && vel(k)<sup_limit(limit)
                                if counter(limit)==0
                                    last=find(counter);
                                    if counter(last)>n_osc*sf;
                                         st{last}(episodes(last))=start;
                                         et{last}(episodes(last))=time(k);
                                         episodes(last)=episodes(last)+1;
                                    end
                                    counter=zeros(5,1);
                                    counter(limit)=1;
                                    start=time(k);
                                else
                                    counter(limit)=counter(limit)+1;
                                end
                         end
                     end
                 end
                 for limit=1:8
                     if (st{limit})
                     for ep=1:length(st{limit})
                         x1=st{limit}(ep); x2=et{limit}(ep);
                         Y{limit}=[Y{limit},eeg(x1:x1+(x2-x1))'];
                     end
                     end
                 end
             end
         end
     end
end

figure(a)
subplot(2,2,1)
[meanp1 peakp1]=COMODOphase(Y{1},n_fs,9,3,18,strcat(int2str(0),'-',int2str(5)));
COMODOphase(Y{1},n_fs,8,3,18,strcat(int2str(0),'-',int2str(5)));
subplot(2,2,2)
[meanp2 peakp2]=COMODOphase(Y{2},n_fs,9,3,18,strcat(int2str(5),'-',int2str(15)));
COMODOphase(Y{2},n_fs,8,3,18,strcat(int2str(5),'-',int2str(15)));
subplot(2,2,3)
[meanp3 peakp3]=COMODOphase(Y{3},n_fs,9,3,18,strcat(int2str(15),'-',int2str(30)));
COMODOphase(Y{3},n_fs,8,3,18,strcat(int2str(15),'-',int2str(30)));
subplot(2,2,4)
[meanp4 peakp4]=COMODOphase(Y{4},n_fs,9,3,18,strcat(int2str(30),'-',int2str(50)));
COMODOphase(Y{4},n_fs,8,3,18,strcat(int2str(30),'-',int2str(50)));

nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\Open\',name(a,:),'_',lfp{a}(i,:));
save(nombre,'como1','como2','como3','-mat')
saveas(a,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\Open\',name(a,:),'_',lfp{a}(i,:)),'png')
close(a)


%% Before starting, load everything

clear all
close all
files=Data_Listing();

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

load('2CA1','-mat');
load('2EC1','-mat');
load('2EC2','-mat');
rat=3;
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

protocol={'A','B','C'};
[n m]=size(area);
a=2;
prot=1;
count=1;
inf_limit=[0,5,15,30];
sup_limit=[5,15,30,50];
Y={0,0,0,0};
for i=1:length(lfp{a})
    st={0,0,0,0};
    et={0,0,0,0};
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
                 counter=zeros(4,1);
                 episodes=ones(4,1);
                 for k=1:length(vel)
                     for limit=1:4
                         if vel(k)>inf_limit(limit) && vel(k)<sup_limit(limit)
                                if counter(limit)==0
                                    last=find(counter);
                                    if counter(last)>n_osc*sf;
                                         st{last}(episodes(last))=start;
                                         et{last}(episodes(last))=time(k);
                                         episodes(last)=episodes(last)+1;
                                    end
                                    counter=zeros(5,1);
                                    counter(limit)=1;
                                    start=time(k);
                                else
                                    counter(limit)=counter(limit)+1;
                                end
                         end
                     end
                 end
                 for limit=1:4
                     if (st{limit})
                     for ep=1:length(st{limit})
                         x1=st{limit}(ep); x2=et{limit}(ep);
                         Y{limit}=[Y{limit},eeg(x1:x1+(x2-x1))'];
                     end
                     end
                 end
             end
         end
     end
end

figure(a)
subplot(2,2,1)
[meanp1 peakp1]=COMODOphase(Y{1},n_fs,9,3,18,strcat(int2str(0),'-',int2str(5)));
COMODOphase(Y{1},n_fs,8,3,18,strcat(int2str(0),'-',int2str(5)));
subplot(2,2,2)
[meanp2 peakp2]=COMODOphase(Y{2},n_fs,9,3,18,strcat(int2str(5),'-',int2str(15)));
COMODOphase(Y{2},n_fs,8,3,18,strcat(int2str(5),'-',int2str(15)));
subplot(2,2,3)
[meanp3 peakp3]=COMODOphase(Y{3},n_fs,9,3,18,strcat(int2str(15),'-',int2str(30)));
COMODOphase(Y{3},n_fs,8,3,18,strcat(int2str(15),'-',int2str(30)));
subplot(2,2,4)
[meanp4 peakp4]=COMODOphase(Y{4},n_fs,9,3,18,strcat(int2str(30),'-',int2str(50)));
COMODOphase(Y{4},n_fs,8,3,18,strcat(int2str(30),'-',int2str(50)));

nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\Open\',name(a,:),'_',lfp{a}(i,:));
save(nombre,'como1','como2','como3','-mat')
saveas(a,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\Open\',name(a,:),'_',lfp{a}(i,:)),'png')
close(a)
