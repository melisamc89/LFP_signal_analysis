%% Before starting, load everything

clear all
close all
files=Data_Listing();

load('2CA1','-mat');
load('2EC1','-mat');
load('2EC2','-mat');
rat=2;
area={CA1_2,EC1_2,EC2_2};
lfp={Days(CA1_2),Days(EC1_2),Days(EC2_2)};
name=['2CA1';'2EC1';'2EC2'];

n_fs=400;
prot='C';

%%
directory='/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/24-05/Speeds/';

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

%% Just for rat 14566
load('1CA1','-mat');
load('1EC1','-mat');
load('1EC2','-mat');
rat=1;
area={CA1,EC1,EC2};
lfp={Days(CA1),Days(EC1),Days(EC2)};
name=['1CA1';'1EC1';'1EC2'];

[n m]=size(area);
totaltime=[20,20,42];
non_mov_time=6;

[n m]=size(area);
prot='C';
non_mov_time=6;
fast_time=5;
slow_time=15;

% matrix1= resting
% matrix2= expecting
% matrix3= running

protocol={'A','B','C'};
[n m]=size(area);

for a=2:2
        clear matrix1 matrix2 matrix3 matrix4
        prot=3;
        %cargar=strcat(directory,'3_',name(a,:),'_COMODO_mean_states');
        %load(cargar,'-mat') 
        count=1;
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
                            como1=comodulogram(x1,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como1);
                            colormap(jet)
                            subplot(2,2,2)
                            como2=comodulogram(x2,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como2);
                            colormap(jet)
                            subplot(2,2,3)
                            como3=comodulogram(x3,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como3);
                            colormap(jet)
                            subplot(2,2,4)
                            como4=comodulogram(x4,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como4);
                            colormap(jet)
                            %color=[max(max(matrix1(count,:))),...
                             %   max(max(matrix2(count,:))),...
                              %  max(max(como3)),...
                               % max(max(como4))];
                           color=[max(max(como1)),max(max(como2)),max(max(como3)),...
                                max(max(como4))];
                            for j=1:4
                               subplot(2,2,j)
                               caxis([0 max(color)])
                            end
                            saveas(i,strcat(directory,name(a,:),'_',lfp{a}(i,:)),'png')
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
        end
        [n m l]=size(matrix1);
        figure(a)
        subplot(2,2,1)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix1),[m l]));
        colormap(jet)
        subplot(2,2,2)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix2),[m l]));
        colormap(jet)
        subplot(2,2,3)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix3),[m l]));
        colormap(jet)
        subplot(2,2,4)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix4),[m l]));
        colormap(jet)
        color=[max(max(mean(matrix1))),max(max(mean(matrix2))),max(max(mean(matrix3))),max(max(mean(matrix4)))];
        for j=1:4
            subplot(2,2,j)
            caxis([0 max(color)])
        end
        nombre=strcat(directory,int2str(prot),'_',name(a,:),'_COMODO_mean_speeds');
        saveas(a,nombre,'png')
        save(nombre,'matrix1','matrix2','matrix3','matrix4','freq_data','-mat')
end


%% Just for rat 14570

load('3CA1','-mat');
load('3EC1','-mat');
load('3EC2','-mat');
rat=3;
area={CA1_3,EC1_3,EC2_3};
lfp={Days(CA1_3),Days(EC1_3),Days(EC2_3)};
name=['3CA1';'3EC1';'3EC2'];

[n m]=size(area);
totaltime=[20,20,42];
non_mov_time=6;
protocol={'A','B','C'};
[n m]=size(area);

[n m]=size(area);
prot='C';
non_mov_time=6;
fast_time=6;
slow_time=30;
n_fs=400;
totaltime=[20,20,36];

for a=2:2
        clear matrix1 matrix2 matrix3 matrix4
        prot=3;
        %cargar=strcat(directory,'3_',name(a,:),'_COMODO_mean_states');
        %load(cargar,'-mat') 
        count=1;
        for i=[1:9,11:19,22:length(lfp{a})]
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
                            como1=comodulogram(x1,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como1);
                            colormap(jet)
                            subplot(2,2,2)
                            como2=comodulogram(x2,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como2);
                            colormap(jet)
                            subplot(2,2,3)
                            como3=comodulogram(x3,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como3);
                            colormap(jet)
                            subplot(2,2,4)
                            como4=comodulogram(x4,lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como4);
                            colormap(jet)
                            %color=[max(max(matrix1(count,:))),...
                             %   max(max(matrix2(count,:))),...
                              %  max(max(como3)),...
                               % max(max(como4))];
                           color=[max(max(como1)),max(max(como2)),max(max(como3)),...
                                max(max(como4))];
                            for j=1:4
                               subplot(2,2,j)
                               caxis([0 max(color)])
                            end
                            saveas(i,strcat(directory,name(a,:),'_',lfp{a}(i,:)),'png')
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
        end
        [n m l]=size(matrix1);
        figure(a)
        subplot(2,2,1)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix1),[m l]));
        colormap(jet)
        subplot(2,2,2)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix2),[m l]));
        colormap(jet)
        subplot(2,2,3)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix3),[m l]));
        colormap(jet)
        subplot(2,2,4)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix4),[m l]));
        colormap(jet)
        color=[max(max(mean(matrix1))),max(max(mean(matrix2))),max(max(mean(matrix3))),max(max(mean(matrix4)))];
        for j=1:4
            subplot(2,2,j)
            caxis([0 max(color)])
        end
        nombre=strcat(directory,int2str(prot),'_',name(a,:),'_COMODO_mean_speeds');
        saveas(a,nombre,'png')
        save(nombre,'matrix1','matrix2','matrix3','matrix4','freq_data','-mat')
end
