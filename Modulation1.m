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
protocol={'A','B','C'};
[n m]=size(area);
a=2;
clear matrix
prot=3;
count=1;
n_fs=400;

[n m]=size(area);
totaltime=[20,20,42];
non_mov_time=6;

for i=[1:19,21:length(lfp{a})]
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
                       [sr sl]=ReadStartTime(pos);
                        s1=sr-non_mov_time;
                        s2=sl-non_mov_time;
                        [srr err srl erl]=RestingTime(pos,totaltime(rat),non_mov_time);
                        s1=floor(s1*n_fs)+1;
                        s2=floor(s2*n_fs)+1;
                        er=floor((sr+totaltime(rat))*n_fs);
                        el=floor((sl+totaltime(rat))*n_fs);                
                        sr=floor(sr*n_fs);
                        sl=floor(sl*n_fs);
                        srr=floor(srr*n_fs)+1;
                        err=floor(err*n_fs);
                        srl=floor(srl*n_fs)+1;
                        erl=floor(erl*n_fs);

                        Y1=[];Y2=[];Y3=[];Y4=[];Y5=[];Y6=[];
                        x1=[];x2=[];x3=[];
                        for trial=1:length(s1)-1
                            Y1=[Y1,eeg(s1(trial):sr(trial))']; %expecting period
                            Y2=[Y2,eeg(s2(trial):sl(trial))'];
                            Y3=[Y3,eeg(sr(trial):er(trial))']; %running period
                            Y4=[Y4,eeg(sl(trial):el(trial))'];
                            Y5=[Y5,eeg(srr(trial):err(trial))']; % resting period
                            Y6=[Y6,eeg(srl(trial):erl(trial))'];
                        end
                        x2=[Y1,Y2];
                        x3=[Y3,Y4];
                        x1=[Y5,Y6];
                        figure(i)
                        subplot(2,2,1)
                        [meanp peakp]=COMODOphase(eeg',n_fs,9,3,18,'All');
                        COMODOphase(eeg',n_fs,10,3,18,'All');
                        subplot(2,2,2)
                        [meanp1 peakp1]=COMODOphase(x1,n_fs,9,3,18,'Resting');
                        COMODOphase(x1,n_fs,10,3,18,'Resting');
                        subplot(2,2,3)
                        [meanp2 peakp2]=COMODOphase(x2,n_fs,9,3,18,'Waiting');
                        COMODOphase(x2,n_fs,10,3,18,'Waiting');
                        subplot(2,2,4)
                        [meanp3 peakp3]=COMODOphase(x3,n_fs,9,3,18,'Running');
                        COMODOphase(x3,n_fs,10,3,18,'Running');
                        saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\States\',name(a,:),'_',lfp{a}(i,:)),'png')
                        close(i)
                        vector1=[meanp meanp1 meanp2 meanp3];
                        vector2=[peakp peakp1 peakp2 peakp3];
                        matrix1(count,:)=vector1;
                        matrix2(count,:)=vector2;
                        count=count+1;
                end
            end
        end
    end
end
clear i
figure(a)
vect1=mean((exp(i*matrix1)));
var_phase1=1-abs(vect1);
mean_phase1=atan2(imag(vect1),real(vect1));
subplot(2,1,1)
h=errorbar([1:4],mean_phase1,var_phase1,'k');
set(h,'Linewidth',2)
xlabel('States: all - Resting - Waiting - Running')
ylabel('Mean Mean Phase')

vect2=mean((exp(i*matrix2)));
var_phase2=1-abs(vect2);
mean_phase2=atan2(imag(vect2),real(vect2));
subplot(2,1,2)
h=errorbar([1:4],mean_phase2,var_phase2,'b');
set(h,'Linewidth',2)
xlabel('States: all - Resting - Waiting - Running')
ylabel('Mean Peak Phase')

nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\States\',int2str(prot),'_',name(a,:),'_COMODO_mean_states');
save(nombre,'matrix1','matrix2','-mat')
saveas(a,nombre,'png')
