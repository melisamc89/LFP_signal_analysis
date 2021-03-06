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
totaltime=20;
jump_time=5*n_fs;

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
                            [sr sl]=ReadStartTime(pos);
                            s1=sr-non_mov_time;
                            s2=sl-non_mov_time;
                            [srr err srl erl]=RestingTime(pos,totaltime,non_mov_time);

                            s1=floor(s1*n_fs)+1;
                            s2=floor(s2*n_fs)+1;
                            er=floor((sr+totaltime)*n_fs);
                            el=floor((sl+totaltime)*n_fs);                
                            srr=floor(srr*n_fs)+1;
                            err=floor(err*n_fs);
                            srl=floor(srl*n_fs)+1;
                            erl=floor(erl*n_fs);
                    
                            sr=floor(sr*n_fs);
                            sl=floor(sl*n_fs);
                            
                        Y1=[];Y2=[];Y3=[];Y4=[];Y5=[];Y6=[];Y7=[];Y8=[];Y9=[];Y10=[];Y11=[];Y12=[];

                            for trial=1:length(s1)-1
                                Y1=[Y1,eeg(srr(trial):err(trial))'];% resting period
                                Y2=[Y2,eeg(srl(trial):erl(trial))'];
                                Y3=[Y3,eeg(s1(trial):sr(trial))'];%expecting period
                                Y4=[Y4,eeg(s2(trial):sl(trial))'];
                                
                                Y5=[Y5,eeg(sr(trial):sr(trial)+jump_time)'];%running1
                                Y6=[Y6,eeg(sr(trial)+jump_time:sr(trial)+2*jump_time)'];
                                Y7=[Y7,eeg(sr(trial)+2*jump_time:sr(trial)+3*jump_time)'];
                                Y8=[Y8,eeg(sr(trial)+3*jump_time:sr(trial)+4*jump_time)'];
                                
                                Y9=[Y9,eeg(sl(trial):sl(trial)+jump_time)'];%running1
                                Y10=[Y10,eeg(sl(trial)+jump_time:sl(trial)+2*jump_time)'];
                                Y11=[Y11,eeg(sl(trial)+2*jump_time:sl(trial)+3*jump_time)'];
                                Y12=[Y12,eeg(sl(trial)+3*jump_time:sl(trial)+4*jump_time)'];
                            end
                        x1=[Y1,Y2];
                        x2=[Y3,Y4];
                        x3=[Y5,Y9];
                        x4=[Y6,Y10];
                        x5=[Y7,Y11];
                        x6=[Y8,Y12];
                        
                        figure(i)
                        subplot(2,3,1)
                        [meanp1 peakp1]=COMODOphase(x1,n_fs,9,3,18,'REsting');
                        COMODOphase(x1,n_fs,10,3,18,'Resting');
                        subplot(2,3,2)
                        [meanp2 peakp2]=COMODOphase(x2,n_fs,9,3,18,'Waiting');
                        COMODOphase(x2,n_fs,10,3,18,'Waiting');
                        subplot(2,3,3)
                        [meanp3 peakp3]=COMODOphase(x3,n_fs,9,3,18,'v1');
                        COMODOphase(x3,n_fs,10,3,18,'v1');
                        subplot(2,3,4)
                        [meanp4 peakp4]=COMODOphase(x4,n_fs,9,3,18,'v2');
                        COMODOphase(x4,n_fs,10,3,18,'v2');
                        subplot(2,3,5)
                        [meanp5 peakp5]=COMODOphase(x5,n_fs,9,3,18,'v3');
                        COMODOphase(x5,n_fs,10,3,18,'v3');
                        subplot(2,3,6)
                        [meanp6 peakp6]=COMODOphase(x6,n_fs,9,3,18,'v4');
                        COMODOphase(x6,n_fs,10,3,18,'v4');
                        saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\Speeds2\',name(a,:),'_',lfp{a}(i,:)),'png')
                        close(i)
                        vector1=[meanp1 meanp2 meanp3 meanp4 meanp5 meanp6];
                        vector2=[peakp1 peakp2 peakp3 peakp4 peakp5 peakp6];
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
h=errorbar([1:6],mean_phase1,var_phase1,'k');
set(h,'Linewidth',2)
xlabel('States: Resting - Waiting - v1 -v2 - v3- v4')
ylabel('Mean Mean Phase')

vect2=mean((exp(i*matrix2)));
var_phase2=1-abs(vect2);
mean_phase2=atan2(imag(vect2),real(vect2));
subplot(2,1,2)
h=errorbar([1:6],mean_phase2,var_phase2,'b');
set(h,'Linewidth',2)
xlabel('States: Resting - Waiting - v1 - v2- v3- v4')
ylabel('Mean Peak Phase')

nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Modulation\Speeds2\',int2str(prot),'_',name(a,:),'_COMODO_mean_states');
save(nombre,'matrix1','matrix2','-mat')
saveas(a,nombre,'png')
