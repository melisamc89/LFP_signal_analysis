%% Before starting, load everything

clear all
close all
files=Data_Listing();

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
totaltime=[20,20,42];
non_mov_time=6;
totaltime=20;
jump_time=5*n_fs;

protocol={'A','B','C'};
[n m]=size(area);
for a=2:2
        clear matrix5 matrix6 matrix7 matrix8 matrix9 matrix10 matrix11 matrix12
        prot=3;
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
                            
                            
                            %for speed=1:4
                             %   for track=1:length(sr)
                              %      t=floor(sr(i)/400*50);
                               %     v(speed,:,track)=pos.data.vx(t+(speed-1)*5*50:t+speed*5*50);
                               % end
                            %end
    
                            
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

                            figure(i)
                            %subplot(2,6,1)
                            %como1=comodulogram(Y1,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como1);
                            %subplot(2,6,7)
                            %como2=comodulogram(Y2,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como2);
                            %subplot(2,5,1)
                            %como3=comodulogram(Y3,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                            %subplot(2,5,6)
                            %como4=comodulogram(Y4,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como4);
                            
                            subplot(2,2,1)
                            como5=comodulogram([Y5,Y12],lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como5);
                            colormap(jet)
                            subplot(2,2,2)
                            como6=comodulogram([Y6,Y11],lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como6);
                            colormap(jet)
                            subplot(2,2,3)
                            como7=comodulogram([Y7,Y10],lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como7);
                            colormap(jet)                            
                            subplot(2,2,4)
                            como8=comodulogram([Y8,Y9],lowf2:fstep:highf2,logspace(-0.01,0.8,number),n_fs,'FilterMethod','wavelet','FilterWidth',[5 3]);
                            COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),como8);
                            colormap(jet)
                            
                            %subplot(2,4,5)
                            %como9=comodulogram(Y12,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterMethod','wavelet');
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como9);
                            %colormap(jet)
                            %subplot(2,4,6)
                            %como10=comodulogram(Y11,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterMethod','wavelet');
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como10);
                            %colormap(jet)
                            %subplot(2,4,7)
                            %como11=comodulogram(Y10,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterMethod','wavelet');
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como11);
                            %colormap(jet)
                            %subplot(2,4,8)
                            %como12=comodulogram(Y9,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs,'FilterMethod','wavelet');
                            %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como12);
                            %colormap(jet)
                            
                            %como=[como5,como6,como7,como8,como9,como10,como11,como12];
                            como=[como5,como6,como7,como8];
                            color=[max(max(como))];
                            for j=1:4
                               subplot(2,2,j)
                               caxis([0 color])
                            end
                            saveas(i,strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/10-05/Speeds2/',name(a,:),'_',lfp{a}(i,:)),'png')
                            close(i)
          matrix5(count,:,:)=como5;matrix6(count,:,:)=como6;
          matrix7(count,:,:)=como7;matrix8(count,:,:)=como8; 
          %matrix9(count,:,:)=como9;matrix10(count,:,:)=como10;
          %matrix11(count,:,:)=como11;matrix12(count,:,:)=como12;

                            count=count+1;
                        end
                    end
                end
            end
        end
        [n m l]=size(matrix5);
        figure(a)
        subplot(2,2,1)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix5),[m l]));
        colormap(jet)
        subplot(2,2,2)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix6),[m l]));
        colormap(jet)
        subplot(2,2,3)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix7),[m l]));
        colormap(jet)        
        subplot(2,2,4)
        COMODOplot(lowf2:fstep:highf2,logspace(-0.01,0.8,number),reshape(mean(matrix8),[m l]));
        colormap(jet)        
        %subplot(2,4,5)
        %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix9),[m l]));
        %subplot(2,4,6)
        %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix10),[m l]));
        %subplot(2,4,7)
        %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix11),[m l]));
        %subplot(2,4,8)
        %COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix12),[m l]));
        como=[mean(matrix5),mean(matrix6),mean(matrix7),mean(matrix8)];
        color=[max(max(como))];
        for j=1:4
            subplot(2,2,j)
            caxis([0 color])
        end
        nombre=strcat('/home/melisa/Escritorio/Melisa/Doctorado/BandsCoupling/10-05/Speeds2/',int2str(prot),'_',name(a,:),'_COMODO_mean_speeds');
        saveas(a,nombre,'png')
        %save(nombre,'matrix5','matrix6','matrix7','matrix8','matrix9','matrix10','matrix11','matrix12','-mat')
        save(nombre,'matrix5','matrix6','matrix7','matrix8','freq_data','-mat')

end

%% 14570

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
slow_time=36;

lowf1=1;
lowf2=6;
highf1=6;
highf2=20;
fstep=0.5;
filter_width=[0.5 0.25];
n_fs=400;

a=2;
        clear matrix1 matrix2 matrix3 matrix4
        prot=3;
        count=1;
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
                                Y5=[Y5,eeg(er1(trial):er(trial))'];%running period slow
                                Y6=[Y6,eeg(sl(trial):el1(trial))'];
                                Y7=[Y7,eeg(sr(trial):er1(trial))'];%running period fast
                                Y8=[Y8,eeg(el1(trial):el(trial))'];
                            end

                            figure(i)
                            subplot(2,2,1)
                            como1=comodulogram(Y5,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como1);
                            subplot(2,2,3)
                            como2=comodulogram(Y6,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como2);
                            subplot(2,2,2)
                            como3=comodulogram(Y7,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como3);
                            subplot(2,2,4)
                            como4=comodulogram(Y8,lowf2:fstep:highf2,lowf1:fstep/2:highf1,n_fs);
                            COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,como4);
                            como=[como1,como2,como3,como4];
                            for j=1:4
                               subplot(2,2,j)
                               caxis([0 max(max(como))])
                            end
                            saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC_DELTA\Prueba\',name(a,:),'_',lfp{a}(i,:)),'png')
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
        COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix1),[m l]));
        subplot(2,2,2)
        COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix2),[m l]));
        subplot(2,2,3)
        COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix3),[m l]));
        subplot(2,2,4)
        COMODOplot(lowf2:fstep:highf2,lowf1:fstep/2:highf1,reshape(mean(matrix4),[m l]));
        color=[max(max(mean(matrix1))),max(max(mean(matrix2))),max(max(mean(matrix3))),max(max(mean(matrix4)))];
        for j=1:4
            subplot(2,2,j)
            caxis([0 max(color)])
        end
        nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\PAC_DELTA\Prueba\',int2str(prot),'_',name(a,:),'_COMODO_mean_speeds');
        saveas(a,nombre,'png')
        save(nombre,'matrix1','matrix2','matrix3','matrix4','freq_data','-mat')


