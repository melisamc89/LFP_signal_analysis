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
prot='C';

%% COMODO inspect

[n m]=size(area);
prot='A';
amplitude_frequency=10;
phase_frequency=4;

for a=1:m
    signal=[];
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs)); 
                pos=Loading_Pos(register,session_number(se));
                tf = isfield(pos.data, 'log');
                if ~(prot=='C') tf=1; end
                if tf
                    signal=[signal,eeg'];
                end
            end
        end
    end
    figure(a)
    COMODOamplitudephasedist(eeg',n_fs,amplitude_frequency,phase_frequency,18);
    amplitude_frequency=10;
    phase_frequency=3;
    COMODOinspect(eeg',n_fs,amplitude_frequency,phase_frequency,18);
    saveas(a,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO\','COMODOinspect',name(a,:)),'png')
    close(a)
end

%% all data coupling

lowf1=2;
lowf2=2;
highf1=6;
highf2=20;
fstep=1;
filter_width=[2 1];
n_fs=400;


[n m]=size(area);
prot='A';
for a=1:m
    count=1;
    clear matrix
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs)); 
                pos=Loading_Pos(register,session_number(se));
                tf = isfield(pos.data, 'log');
                if ~(prot=='C') tf=1; end
                if tf==1
                    figure(i)
                    %como=comodulogram(eeg',lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    como=comodulogram(eeg',lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterWidth',filter_width);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como);
                    %caxis([0 0.003])
                    saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\Rat14570\',name(a,:),'_',lfp{a}(i,:)),'png')
                    close(i)
                    matrix(count,:,:)=como;
                    count=count+1;
                end
            end
        end
    end
    figure(10*a)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\Rat14570\','COMODO_mean_',name(a,:));
    save(nombre,'matrix','-mat')
    saveas(10*a,nombre,'png')
    close(10*a)
end

%% trial concatenation

[n m]=size(area);
prot='C';
totaltime=20;
non_mov_time=6;

lowf1=2;
lowf2=2;
highf1=6;
highf2=20;
fstep=1;
filter_width=[2 1];
n_fs=400;

for a=1:m
    count=1;
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs));
                pos=Loading_Pos(register,session_number(se));
                tf = isfield(pos.data, 'log');
                if tf==1
                    [sr sl]=ReadStartTime(pos);
                    s1=sr-non_mov_time;
                    s2=sl-non_mov_time;
                    [srr err srl erl]=RestingTime(pos,totaltime,non_mov_time);

                    s1=floor(s1*n_fs)+1;
                    s2=floor(s2*n_fs)+1;
                    er=floor((sr+totaltime)*n_fs);
                    el=floor((sl+totaltime)*n_fs);                
                    sr=floor(sr*n_fs);
                    sl=floor(sl*n_fs);
                    srr=floor(srr*n_fs)+1;
                    err=floor(err*n_fs);
                    srl=floor(srl*n_fs)+1;
                    erl=floor(erl*n_fs);
                    
                    Y1=[];Y2=[];Y3=[];Y4=[];Y5=[];Y6=[];

                    for trial=1:length(s1)-1
                        Y1=[Y1,eeg(s1(trial):sr(trial))']; %expecting period
                        Y2=[Y2,eeg(s2(trial):sl(trial))'];
                        Y3=[Y3,eeg(sr(trial):er(trial))']; %running period
                        Y4=[Y4,eeg(sl(trial):el(trial))'];
                        Y5=[Y5,eeg(srr(trial):err(trial))']; % resting period
                        Y6=[Y6,eeg(srl(trial):erl(trial))'];
                    end
                    clear s1 s2 sr sl er el
                    x1=[Y1,Y2];
                    x2=[Y3,Y4];
                    x3=[Y5,Y6];
                    figure(i)
                    subplot(2,2,1)
                    como=comodulogram(eeg',lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como);
                    %caxis([0 0.02])
                    %saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO\ProtocolC\Rat1\','COMODO_',lfp{a}(i,:)),'png')
                    %close(i)
                    %figure(i)
                    subplot(2,2,2)
                    como1=comodulogram(x1,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como1);
                    subplot(2,2,3)
                    como2=comodulogram(x2,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como2);
                    subplot(2,2,4)
                    como3=comodulogram(x3,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                    color=[max(max(como1)),max(max(como2)),max(max(como3)),max(max(como))];
                    for j=1:4
                      subplot(2,2,j)
                       caxis([0 max(color)])
                    end
                    saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\Rat14566_cat\',lfp{a}(i,:)),'png')
                    close(i)
                    matrix(count,:,:)=como;
                    matrix1(count,:,:)=como1;
                    matrix2(count,:,:)=como2;
                    matrix3(count,:,:)=como3;
                    count=count+1;
                end
            end
        end
    end
    figure(10*a)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix),[10 25]));
    figure(a)
    subplot(2,2,2)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix1),[10 25]));
    subplot(2,2,3)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix2),[10 25]));
    subplot(2,2,4)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix3),[10 25]));
    color=[max(max(como1)),max(max(como2)),max(max(como3))];
    for j=2:4
        subplot(2,2,j)
        caxis([0 max(color)])
    end
    name=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\','COMODO_mean');
    save(name,'matrix','-mat')
    name=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO\ProtocolC\Rat1\','COMODO_mean_states');
    save(name,'matrix1','matrix2','matrix3','-mat')
    clear matrix matrix1 matrix2 matrix3
end

%% different state and speed conditions

[n m]=size(area);
prot='C';
totaltime=20;
non_mov_time=6;
fast_time=5;
slow_time=15;

for a=1:m
    count=1;
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs));
                pos=Loading_Pos(register,session_number(se));
                tf = isfield(pos.data, 'log');
                if tf==1
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
                    
                    er1=floor((sr+fast_time)*n_fs);
                    el1=floor((sl+slow_time)*n_fs);
                    sr=floor(sr*n_fs);
                    sl=floor(sl*n_fs);
                    
                    Y1=[];Y2=[];Y3=[];Y4=[];Y5=[];Y6=[];Y7=[];Y8=[];

                    for trial=1:length(s1)-1
                        Y1=[Y1,eeg(srr(trial):err(trial))'];% resting period
                        Y2=[Y2,eeg(srl(trial):erl(trial))'];
                        Y3=[Y3,eeg(s1(trial):sr(trial))'];%expecting period
                        Y4=[Y4,eeg(s2(trial):sl(trial))'];                        
                        Y5=[Y5,eeg(er1(trial):er(trial))'];%running period slow
                        Y6=[Y6,eeg(sl(trial):el1(trial))'];
                        Y7=[Y7,eeg(sr(trial):er1(trial))'];%running period fast
                        Y8=[Y8,eeg(el1(trial):el(trial))'];

                    end
                    x1=[Y1,Y2];
                    x2=[Y3,Y4];
                    x3=[Y5,Y6];
                    x4=[Y7,Y8];
                  
                    lowf2=40;
                    fstep=4;
                    highf2=100;
                    lowf1=5;
                    highf1=15;
                    
                    figure(i)
                    subplot(2,2,1)
                    como1=comodulogram(x1,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como1);
                    subplot(2,2,2)
                    como2=comodulogram(x2,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como2);
                    subplot(2,2,3)
                    como3=comodulogram(x3,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                    subplot(2,2,4)
                    como4=comodulogram(x4,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                    color=[max(max(como1)),max(max(como2)),max(max(como3)),max(max(como4))];
                    for j=1:4
                       subplot(2,2,j)
                       caxis([0 max(color)])
                    end
                    saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\Rat14566_cat_delta\',name(a,:),'_',lfp{a}(i,:)),'png')
                    close(i)
                    matrix1(count,:,:)=como1;
                    matrix2(count,:,:)=como2;
                    matrix3(count,:,:)=como3;
                    matrix4(count,:,:)=como4;
                    count=count+1;
                end
            end
        end
    end
    figure(a)
    subplot(2,2,1)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix1),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    subplot(2,2,2)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix2),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    subplot(2,2,3)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix3),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    subplot(2,2,4)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix4),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    for j=1:4
        subplot(2,2,j)
        caxis([0 0.0007])
    end
    nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\Rat14566_cat_delta\','COMODO_mean_speads_',name(a,:));
    save(nombre,'matrix1','matrix2','matrix3','matrix4','-mat')
    saveas(a,nombre,'png')
    clear matrix1 matrix2 matrix3 matrix4
end


%% for delta

%% different state and speed conditions

[n m]=size(area);
prot='C';
totaltime=40;
non_mov_time=6;
fast_time=5;
slow_time=35;

for a=1:m
    count=1;
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs));
                pos=Loading_Pos(register,session_number(se));
                tf = isfield(pos.data, 'log');
                if tf==1
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
                    
                    er1=floor((sr+fast_time)*n_fs);
                    el1=floor((sl+slow_time)*n_fs);
                    sr=floor(sr*n_fs);
                    sl=floor(sl*n_fs);
                    
                    Y1=[];Y2=[];Y3=[];Y4=[];Y5=[];Y6=[];Y7=[];Y8=[];

                    for trial=1:length(s1)-1
                        Y1=[Y1,eeg(srr(trial):err(trial))'];% resting period
                        Y2=[Y2,eeg(srl(trial):erl(trial))'];
                        Y3=[Y3,eeg(s1(trial):sr(trial))'];%expecting period
                        Y4=[Y4,eeg(s2(trial):sl(trial))'];
                        
                        Y5=[Y5,eeg(er1(trial):er(trial))'];%running period slow
                        Y6=[Y6,eeg(sl(trial):el1(trial))'];
                        Y7=[Y7,eeg(sr(trial):er1(trial))'];%running period fast
                        Y8=[Y8,eeg(el1(trial):el(trial))'];

                    end
                    x1=[Y1,Y2];
                    x2=[Y3,Y4];
                    x3=[Y5,Y6];
                    x4=[Y7,Y8];
                  
                    lowf2=2;
                    fstep=2;
                    highf2=20;
                    lowf1=2;
                    highf1=6;
                    
                    figure(i)
                    subplot(2,2,1)
                    como1=comodulogram(x1,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterWidth',[1 1]);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como1);
                    subplot(2,2,2)
                    como2=comodulogram(x2,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterWidth',[1 1]);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como2);
                    subplot(2,2,3)
                    como3=comodulogram(x3,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterWidth',[1 1]);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                    subplot(2,2,4)
                    como4=comodulogram(x4,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterWidth',[1 1]);
                    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                    color=[max(max(como1)),max(max(como2)),max(max(como3)),max(max(como4))];
                    for j=1:4
                       subplot(2,2,j)
                       caxis([0 max(color)])
                    end
                    saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\Rat14570_cat_delta\',name(a,:),'_',lfp{a}(i,:)),'png')
                    close(i)
                    matrix1(count,:,:)=como1;
                    matrix2(count,:,:)=como2;
                    matrix3(count,:,:)=como3;
                    matrix4(count,:,:)=como4;
                    count=count+1;
                end
            end
        end
    end
    figure(a)
    subplot(2,2,1)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix1),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    subplot(2,2,2)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix2),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    subplot(2,2,3)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix3),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    subplot(2,2,4)
    COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix4),[length(lowf1:fstep/2:highf1) length(lowf2:fstep:highf2)]));
    for j=1:4
        subplot(2,2,j)
        caxis([0 0.0007])
    end
    nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\COMODO2\Rat14570_cat_delta\','COMODO_mean_speads_',name(a,:));
    save(nombre,'matrix1','matrix2','matrix3','matrix4','-mat')
    saveas(a,nombre,'png')
    clear matrix1 matrix2 matrix3 matrix4
end
