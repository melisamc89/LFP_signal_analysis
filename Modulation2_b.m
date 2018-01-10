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

%% rat 14566

fast_time=5;
slow_time=15;

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
                       [sr sl]=ReadStartTime(pos);
                            s1=sr-non_mov_time;
                            s2=sl-non_mov_time;
                            [srr err srl erl]=RestingTime(pos,totaltime(rat),non_mov_time);

                            s1=floor(s1*n_fs)+1;
                            s2=floor(s2*n_fs)+1;
                            er=floor((sr+totaltime(rat))*n_fs);
                            el=floor((sl+totaltime(rat))*n_fs);                
                            srr=floor(srr*n_fs)+1;
                            err=floor(err*n_fs);
                            srl=floor(srl*n_fs)+1;
                            erl=floor(erl*n_fs);

                            er1=floor((sr+slow_time)*n_fs);
                            el1=floor((sl+fast_time)*n_fs);
                            sr=floor(sr*n_fs);
                            sl=floor(sl*n_fs);

                        Y1=[];Y2=[];Y3=[];Y4=[];Y5=[];Y6=[];Y7=[];Y8=[];
                        for trial=1:length(s1)-1
                            Y3=[Y3,eeg(s1(trial):sr(trial))']; %expecting period
                            Y4=[Y4,eeg(s2(trial):sl(trial))'];
                            Y1=[Y1,eeg(srr(trial):err(trial))']; % resting period
                            Y2=[Y2,eeg(srl(trial):erl(trial))'];
                            Y5=[Y5,eeg(sr(trial):er1(trial))'];%running period slow
                            Y6=[Y6,eeg(el1(trial):el(trial))'];
                            Y7=[Y7,eeg(er1(trial):er(trial))'];%running period fast
                            Y8=[Y8,eeg(sl(trial):el1(trial))'];
                        end
                        x1=[Y1,Y2];
                        x2=[Y3,Y4];
                        x3=[Y5,Y6];
                        x4=[Y7,Y8];
                        
                       figure(i)
                        subplot(2,2,1)
                        [meanp1 peakp1]=COMODOphase2(x1,n_fs,9,3,18,'REsting');
                        COMODOphase2(x1,n_fs,10,3,18,'Expecting');
                        subplot(2,2,2)
                        [meanp2 peakp2]=COMODOphase2(x2,n_fs,9,3,18,'Waiting');
                        COMODOphase2(x2,n_fs,10,3,18,'Resting');
                        subplot(2,2,3)
                        [meanp3 peakp3]=COMODOphase2(x3,n_fs,9,3,18,'Running Slow');
                        COMODOphase2(x3,n_fs,10,3,18,'Slow');
                        subplot(2,2,4)
                        [meanp4 peakp4]=COMODOphase2(x4,n_fs,9,3,18,'Running Fast');
                        COMODOphase2(x4,n_fs,10,3,18,'Fast');
                        saveas(i,strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/Modulation2/Speeds/',name(a,:),'_',lfp{a}(i,:)),'png')
                        close(i)
                        vector1=[meanp1 meanp2 meanp3 meanp4];
                        vector2=[peakp1 peakp2 peakp3 peakp4];
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
xlabel('States: Resting - Waiting - Running Slow - Running Fast')
ylabel('Mean Mean Phase')

vect2=mean((exp(i*matrix2)));
var_phase2=1-abs(vect2);
mean_phase2=atan2(imag(vect2),real(vect2));
subplot(2,1,2)
h=errorbar([1:4],mean_phase2,var_phase2,'b');
set(h,'Linewidth',2)
xlabel('States: Resting - Waiting - Running Slow - Running Fast')
ylabel('Mean Peak Phase')

nombre=strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/Modulation2/Speeds/',int2str(prot),'_',name(a,:),'_COMODO_mean_states');
save(nombre,'matrix1','matrix2','-mat')
saveas(a,nombre,'png')

%% rat 14570
non_mov_time=6;
fast_time=6;
slow_time=36;

for i=21:length(lfp{a})
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
                            srr=floor(srr*n_fs)+1;
                            err=floor(err*n_fs);
                            srl=floor(srl*n_fs)+1;
                            erl=floor(erl*n_fs);
                            
                            er1=floor((sr+fast_time)*n_fs);
                            el1=floor((sl+slow_time)*n_fs);
                            sr=floor(sr*n_fs);
                            sl=floor(sl*n_fs);

                        Y1=[];Y2=[];Y3=[];Y4=[];Y5=[];Y6=[];Y7=[];Y8=[];
                        for trial=1:length(s1)-1
                            Y3=[Y3,eeg(s1(trial):sr(trial))']; %expecting period
                            Y4=[Y4,eeg(s2(trial):sl(trial))'];
                            Y1=[Y1,eeg(srr(trial):err(trial))']; % resting period
                            Y2=[Y2,eeg(srl(trial):erl(trial))'];                            
                            Y5=[Y5,eeg(er1(trial):er(trial))'];%running period slow
                            Y6=[Y6,eeg(sl(trial):el1(trial))'];
                            Y7=[Y7,eeg(sr(trial):er1(trial))'];%running period fast
                            Y8=[Y8,eeg(el1(trial):el(trial))'];
                        end
                        x1=[Y1,Y2];
                        x2=[Y3,Y4];
                        x3=[Y5,Y6];
                        x4=[Y7,Y8];
                        
                        figure(i)
                        subplot(2,2,1)
                        [meanp1 peakp1]=COMODOphase2(x1,n_fs,9,3,18,'REsting');
                        COMODOphase2(x1,n_fs,10,3,18,'Expecting');
                        subplot(2,2,2)
                        [meanp2 peakp2]=COMODOphase2(x2,n_fs,9,3,18,'Waiting');
                        COMODOphase2(x2,n_fs,10,3,18,'Resting');
                        subplot(2,2,3)
                        [meanp3 peakp3]=COMODOphase2(x3,n_fs,9,3,18,'Running Slow');
                        COMODOphase2(x3,n_fs,10,3,18,'Slow');
                        subplot(2,2,4)
                        [meanp4 peakp4]=COMODOphase2(x4,n_fs,9,3,18,'Running Fast');
                        COMODOphase2(x4,n_fs,10,3,18,'Fast');
                        saveas(i,strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/Modulation2/Speeds/',name(a,:),'_',lfp{a}(i,:)),'png')
                        close(i)
                        vector1=[meanp1 meanp2 meanp3 meanp4];
                        vector2=[peakp1 peakp2 peakp3 peakp4];
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
xlabel('States: Resting - Waiting - Running Slow - Running Fast')
ylabel('Mean Mean Phase')

vect2=mean((exp(i*matrix2)));
var_phase2=1-abs(vect2);
mean_phase2=atan2(imag(vect2),real(vect2));
subplot(2,1,2)
h=errorbar([1:4],mean_phase2,var_phase2,'b');
set(h,'Linewidth',2)
xlabel('States: Resting - Waiting - Running Slow - Running Fast')
ylabel('Mean Peak Phase')

nombre=strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/Modulation2/Speeds/',int2str(prot),'_',name(a,:),'_COMODO_mean_states');
save(nombre,'matrix1','matrix2','-mat')
saveas(a,nombre,'png')