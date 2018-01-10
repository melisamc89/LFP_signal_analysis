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

%%


[n m]=size(area);
totaltime=[20,20,42];
non_mov_time=6;

protocol={'A','B','C'};
[n m]=size(area);

        clear matrix1 matrix2 matrix3
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
                        BW2=3;
                        filt_eeg= COMODOfilter(eeg',n_fs,9,BW2,'eegfilt');
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
                                Y1=[Y1,filt_eeg(s1(trial):sr(trial))]; %expecting period
                                Y2=[Y2,filt_eeg(s2(trial):sl(trial))];
                                Y3=[Y3,filt_eeg(sr(trial):er(trial))]; %running period
                                Y4=[Y4,filt_eeg(sl(trial):el(trial))];
                                Y5=[Y5,filt_eeg(srr(trial):err(trial))]; % resting period
                                Y6=[Y6,filt_eeg(srl(trial):erl(trial))];
                            end
                            x2=[Y1,Y2];
                            x3=[Y3,Y4];
                            x1=[Y5,Y6];
                            figure(i)
                            WINDOWSIZE=1;
                            SRATE=n_fs;
                            [~,peaks] = findpeaks(real(x1));
                            peaks=downsample(peaks,37); 
                            [aaa,bbb] = meshgrid(-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs),peaks);
                            ppp = aaa+bbb;
                            ppp(ppp<=0 | ppp>length(x1)) = length(x1+1);
                            vect1 = mean(x1(ppp));
                            plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,vect1,'b','linewidth',2);
                            hold on
                            xlim([-WINDOWSIZE +WINDOWSIZE]);
                            xlabel('Time (s)');
                            
                            [~,peaks] = findpeaks(real(x2));
                            peaks=downsample(peaks,37); 
                            [aaa,bbb] = meshgrid(-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs),peaks);
                            ppp = aaa+bbb;
                            ppp(ppp<=0 | ppp>length(x2)) = length(x2+1);
                            vect2 = mean(x2(ppp));
                            plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,vect2,'r','linewidth',2);
                            
                            [~,peaks] = findpeaks(real(x3));
                            peaks=downsample(peaks,37); 
                            [aaa,bbb] = meshgrid(-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs),peaks);
                            ppp = aaa+bbb;
                            ppp(ppp<=0 | ppp>length(x3)) = length(x3+1);
                            vect3 = mean(x3(ppp));
                            plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,vect3,'k','linewidth',2);
                            legend('Resting','Waiting','Running')
                            
                            saveas(i,strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Phase\Delta\',name(a,:),'_',lfp{a}(i,:)),'png')
                            close(i)
                            matrix1(count,:,:)=vect1;
                            matrix2(count,:,:)=vect2;
                            matrix3(count,:,:)=vect3;
                            count=count+1;
                        end
                    end
                end
            end
        end
        figure(a)
        [n m l]=size(matrix1);
        plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,reshape(mean(matrix1),[l 1]),'b','linewidth',2);
        hold on
        plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,reshape(mean(matrix2),[l 1]),'r','linewidth',2);
        plot((-round(WINDOWSIZE*n_fs):round(WINDOWSIZE*n_fs))/n_fs,reshape(mean(matrix3),[l 1]),'k','linewidth',2);
        xlabel('Time [s]')
        legend('Resting','Waiting','Running')
        nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\BandsCoupling\Phase\Delta\',int2str(prot),'_',name(a,:),'_COMODO_mean_states');
        saveas(a,nombre,'png')
        save(nombre,'matrix1','matrix2','matrix3','-mat')


