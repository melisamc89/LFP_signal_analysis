%% Before starting, load everything

clear all
close all
files=Data_Listing();

rat=1;
load('1CA1','-mat');
load('1EC1','-mat');
load('1EC2','-mat');
area={CA1,EC1,EC2};
lfp={Days(CA1),Days(EC1),Days(EC2)};
name=['1CA1';'1EC1';'1EC2'];


rat=2;
load('2CA1','-mat');
load('2EC1','-mat');
load('2EC2','-mat');
area={CA1_2,EC1_2,EC2_2};
lfp={Days(CA1_2),Days(EC1_2),Days(EC2_2)};
name=['2CA1';'2EC1';'2EC2'];

rat=3;
load('3CA1','-mat');
load('3EC1','-mat');
load('3EC2','-mat');
area={CA1_3,EC1_3,EC2_3};
lfp={Days(CA1_3),Days(EC1_3),Days(EC2_3)};
name=['3CA1';'3EC1';'3EC2'];

n_fs=400;
Fs = n_fs;            % Sampling frequency
T = 1/Fs;             % Sampling period

%% Mean fft spectrum
protocol={'A','B','C','IC'};
[n m]=size(area);
for a=1:3
    for prot=1:3
    counter=1;
    clear matrix
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),protocol{prot});
        if length(session_number)
            order=Wished_Register_Order(area{a},lfp{a}(i,:));
            register=Rat_Register(files,area{a},rat,order);
            [eeg,fs]=ReadEEG(register,session_number(1));
            [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
            if qual>0.5
                n_fs=floor(fs/floor(fs/n_fs));
                eeg=resample(eeg,1,floor(fs/n_fs));
                n = 2^10;
                Y=fft(eeg,n);
                L=length(Y);
                f = Fs*(0:(L-1))/L;
                matrix(counter,:)=abs(Y/sqrt(L));
                counter=counter+1;
            end
        end
    end

    [n1 n2]=size(matrix);
    if matrix
    P2 = mean(matrix);
    vP2=sqrt(var(matrix));
    mP2=vP2/sqrt(n1);
    P1 = P2(1:L/2+1);
    vP1=vP2(1:L/2+1);
    mP1=mP2(1:L/2+1);
    P1(2:end-1)= 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    P1=P1'; f=f'; vP1=vP1'; mP1=mP1';
    data=[f,P1,vP1,mP1];
    
    %figure!
    figure(a)
    errorbar(f,P1,vP1,'Color',[192 192 192]/225)
    hold on
    errorbar(f,P1,mP1,'Color',[128 128 128]/225)
    axis([0 150 0 max(P1(10:end))])
    plot(f,P1,'k','Linewidth',2);
    xlabel('Frequency [Hz]','FontSize',15);
    ylabel('Mean PSD','FontSize',15);
    box on

    axes('Position',[.5 .5 .35 .35])
    box on
    errorbar(f,P1,vP1,'Color',[192 192 192]/225)
    hold on
    errorbar(f,P1,mP1,'Color',[128 128 128]/225)
    axis([0.5 15 0 max(P1(10:end))])
    plot(f,P1,'k','Linewidth',2);
    xlabel('Frequency [Hz]','FontSize',10);
    ylabel('Mean PSD','FontSize',10);
    nombre=strcat('/home/melisa/Escritorio/Melisa/Doctorado/Spectrogram/Figuras2/',int2str(rat),'_',name(a,:),'_',protocol{prot});
    save_data=[f,P1,vP1,mP1];
    save(nombre,'save_data','-ASCII')
    saveas(a,nombre,'png')
    close(a)
    end
    end
end

%% Comparaci�n de espectros en diferentes din�micas. Carrito!

prot='C';
[n m]=size(area);
totaltime=[20,20,42];
non_mov_time=6;

for a=1:m
    counter=0;
    clear matrix1 matrix2 matrix3 matrix4 matrix5 matrix6 matrix
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
                if qual>0.5
                    eeg=resample(eeg,1,floor(fs/n_fs));
                    pos=Loading_Pos(register,session_number(se));
                    tf = isfield(pos.data, 'log');
                    if tf==1
                        [sr sl]=ReadStartTime(pos);
                        s1=sr-non_mov_time;
                        s2=sl-non_mov_time;
                        [srr err srl erl]=RestingTime(pos,totaltime(rat),non_mov_time);

                        s1=floor(s1*n_fs)+1;
                        s2=floor(s2*n_fs);
                        er=floor((sr+totaltime(rat))*n_fs);
                        el=floor((sl+totaltime(rat))*n_fs);                
                        sr=floor(sr*n_fs);
                        sl=floor(sl*n_fs);
                        srr=floor(srr*n_fs)+1;
                        err=floor(err*n_fs);
                        srl=floor(srl*n_fs)+1;
                        erl=floor(erl*n_fs);

                        n= 2^10;

                        for trial=1:length(s1)-1
                            counter=counter+1;
                            Y1=fft(eeg(s1(trial):sr(trial)),n);
                            Y2=fft(eeg(s2(trial):sl(trial)),n);
                            Y3=fft(eeg(sr(trial):er(trial)),n);
                            Y4=fft(eeg(sl(trial):el(trial)),n);
                            Y5=fft(eeg(srr(trial):err(trial)),n);
                            Y6=fft(eeg(srl(trial):erl(trial)),n);
                            L1=length(Y1);
                            L2=length(Y3);
                            L3=length(Y5);
                            f1 = Fs*(0:(L1-1))/L1;
                            f2 = Fs*(0:(L2-1))/L2;                    
                            f3 = Fs*(0:(L3-1))/L3;    
                            matrix1(counter,1:length(Y1))=abs(Y1/sqrt(L1));
                            matrix2(counter,1:length(Y2))=abs(Y2/sqrt(L1));
                            matrix3(counter,1:length(Y3))=abs(Y3/sqrt(L2));
                            matrix4(counter,1:length(Y4))=abs(Y4/sqrt(L2));
                            matrix5(counter,1:length(Y5))=abs(Y5/sqrt(L3));
                            matrix6(counter,1:length(Y6))=abs(Y6/sqrt(L3));
                        end
                    end
                end
            end
        end
    end
    
    matrix{2}=[matrix1;matrix2]; % expecting
    matrix{3}=[matrix3;matrix4]; %moving
    matrix{1}=[matrix5;matrix6]; %quiet
    L=[L3,L1,L2];
    colors=['k','b','r'];
    figure(a)
    for state=1:3
        [n1 m1]=size(matrix{state});
        P2 = mean(matrix{state});
        vP2=sqrt(var(matrix{state}));
        mP2=vP2/sqrt(n1);
        P1 = P2(1:L(state)/2+1);
        vP1=P2(1:L(state)/2+1);
        P1(2:end-1)= 2*P1(2:end-1);
        f = Fs*(0:(L(state)/2))/L(state);
        hold on
        plot(f,P1,colors(state),'Linewidth',2);
    end
    axis([1 150 0 2])
    xlabel('Frequency [Hz]','FontSize',20);
    ylabel('Mean PSD','FontSize',20);
    legend('Quiet','Waiting','Moving','FontSize',20)
    box on
    axes('Position',[.5 .35 .35 .35])
    box on
    for state=1:3
        [n1 m1]=size(matrix{state});
        P2 = mean(matrix{state});
        vP2=sqrt(var(matrix{state}));
        mP2=vP2/sqrt(n1);
        P1 = P2(1:L(state)/2+1);
        vP1=P2(1:L(state)/2+1);
        P1(2:end-1)= 2*P1(2:end-1);
        f = Fs*(0:(L(state)/2))/L(state);
        plot(f,P1,colors(state),'Linewidth',2);
        hold on
        xlabel('Frequency [Hz]','FontSize',10);
        ylabel('Mean PSD','FontSize',10);
        axis([0 25 0 2])
    end
    nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\Spectrogram\Figuras2\','states_',int2str(rat),'_',name(a,:),'_',prot);
    saveas(a,nombre,'png')
    close(a)
end
    
%% Comparaci�n de espectros en estados de movimiento (dos velocidades). Esto es solo para la rata 3!!!!!

rat=3;
load('3CA1','-mat');
load('3EC1','-mat');
load('3EC2','-mat');
area={CA1_3,EC1_3,EC2_3};
lfp={Days(CA1_3),Days(EC1_3),Days(EC2_3)};
name=['3CA1';'3EC1';'3EC2'];

[n m]=size(area);
totaltime=[20,20,42];
% rat 3
fast_time=6;
slow_time=36;
non_mov_time=6;
prot='C';

for a=1:m
    counter=0;
    clear matrix1 matrix2 matrix3 matrix4 matrix5 matrix6 matrix7 matrix8 matrix
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
                if qual>0.5
                    eeg=resample(eeg,1,floor(fs/n_fs));
                    pos=Loading_Pos(register,session_number(se));
                    tf = isfield(pos.data, 'log');
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
                        n= 2^10;

                        for trial=1:length(sl)-1
                            counter=counter+1;
                            Y1=fft(eeg(srr(trial):err(trial)),n);
                            Y2=fft(eeg(srl(trial):erl(trial)),n);
                            Y3=fft(eeg(s1(trial):sr(trial)),n);
                            Y4=fft(eeg(s2(trial):sl(trial)),n);

                            Y5=fft(eeg(sr(trial):er1(trial)),n);
                            Y6=fft(eeg(el1(trial):el(trial)),n);
                            Y7=fft(eeg(er1(trial):er(trial)),n);
                            Y8=fft(eeg(sl(trial):el1(trial)),n);

                            L1=length(Y1);
                            L2=length(Y3);
                            L3=length(Y5);
                            L4=length(Y7);
                            f1 = Fs*(0:(L1-1))/L1;
                            f2 = Fs*(0:(L2-1))/L2;
                            f3 = Fs*(0:(L3-1))/L3;
                            f4 = Fs*(0:(L4-1))/L4;                    
                            matrix1(counter,1:length(Y1))=abs(Y1/sqrt(L1));
                            matrix2(counter,1:length(Y2))=abs(Y2/sqrt(L1));
                            matrix3(counter,1:length(Y3))=abs(Y3/sqrt(L2));
                            matrix4(counter,1:length(Y4))=abs(Y4/sqrt(L2));
                            matrix5(counter,1:length(Y5))=abs(Y5/sqrt(L3));
                            matrix6(counter,1:length(Y6))=abs(Y6/sqrt(L3));
                            matrix7(counter,1:length(Y7))=abs(Y7/sqrt(L4));
                            matrix8(counter,1:length(Y8))=abs(Y8/sqrt(L4));
                        end
                    end
                end
            end
        end
    end
    matrix{1}=[matrix1;matrix2]; %quiet
    matrix{2}=[matrix3;matrix4]; %expecting
    matrix{3}=[matrix7;matrix8]; %slow
    matrix{4}=[matrix5;matrix6]; %fast
    
    L=[L1,L2,L4,L3];
    colors=['k','b','r','g'];
    figure(a)
    for state=1:4
        [n1 m1]=size(matrix{state});
        P2 = mean(matrix{state});
        vP2=sqrt(var(matrix{state}));
        mP2=vP2/sqrt(n1);
        P1 = P2(1:L(state)/2+1);
        vP1=P2(1:L(state)/2+1);
        P1(2:end-1)= 2*P1(2:end-1);
        f = Fs*(0:(L(state)/2))/L(state);
        hold on
        plot(f,P1,colors(state),'Linewidth',1);
    end
    axis([1 150 0 2])
    xlabel('Frequency [Hz]','FontSize',20);
    ylabel('Mean PSD','FontSize',20);    box on
    legend('Quiet','Waiting','Slow','Fast','FontSize',20)
    box on
    axes('Position',[.5 .35 .35 .35])
    box on
    for state=1:4
        [n1 m1]=size(matrix{state});
        P2 = mean(matrix{state});
        vP2=sqrt(var(matrix{state}));
        mP2=vP2/sqrt(n1);
        P1 = P2(1:L(state)/2+1);
        vP1=P2(1:L(state)/2+1);
        P1(2:end-1)= 2*P1(2:end-1);
        f = Fs*(0:(L(state)/2))/L(state);
        plot(f,P1,colors(state));
        hold on
        xlabel('Frequency [Hz]','FontSize',10);
        ylabel('Mean PSD','FontSize',10);
        axis([0 15 0 2])
    end
    nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\Spectrogram\Figuras2\','speeds_',int2str(rat),'_',name(a,:),'_',prot);
    saveas(a,nombre,'png')
    close(a)
    
end

%% Many speeds analysis. Esto es solo para la rata 1 o 2

rat=1;
load('1CA1','-mat');
load('1EC1','-mat');
load('1EC2','-mat');
area={CA1,EC1,EC2};
lfp={Days(CA1),Days(EC1),Days(EC2)};
name=['1CA1';'1EC1';'1EC2'];

[n m]=size(area);
totaltime=[20,20,42];
non_mov_time=6;
jump_time=4;

for a=1:m
    counter=1;
    clear matrix matrix1 matrix2 matrix3 matrix4 matrix5 matrix6 matrix7 matrix8 matrix9 matrix10 matrix11 matrix12 matrix13 matrix14
    for i=1:length(lfp{a})
        session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),prot);
        if length(session_number)
            for se=1:length(session_number)
                order=Wished_Register_Order(area{a},lfp{a}(i,:));
                register=Rat_Register(files,area{a},rat,order);
                [eeg,fs]=ReadEEG(register,session_number(se));
                n_fs=floor(fs/floor(fs/n_fs));
                [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
                if qual>0.5
                    eeg=resample(eeg,1,floor(fs/n_fs));
                    pos=Loading_Pos(register,session_number(se));
                    tf = isfield(pos.data, 'log');
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

                        er1=floor((sr+jump_time)*n_fs);
                        el1=floor((sl+jump_time)*n_fs);
                        er2=floor(er1+jump_time*n_fs);
                        el2=floor(el1+jump_time*n_fs);
                        er3=floor(er2+jump_time*n_fs);
                        el3=floor(el2+jump_time*n_fs);
                        er4=floor(er3+jump_time*n_fs);
                        el4=floor(el3+jump_time*n_fs);
                        sr=floor(sr*n_fs)+1;
                        sl=floor(sl*n_fs);
                        n= 2^10;

                        for trial=1:length(sl)-1
                            Y1=fft(eeg(srr(trial):err(trial)),n);
                            Y2=fft(eeg(srl(trial):erl(trial)),n);
                            Y3=fft(eeg(s1(trial):sr(trial)),n);
                            Y4=fft(eeg(s2(trial):sl(trial)),n);

                            Y5=fft(eeg(sr(trial):er1(trial)),n);
                            Y6=fft(eeg(er1(trial):er2(trial)),n);
                            Y7=fft(eeg(er2(trial):er3(trial)),n);
                            Y8=fft(eeg(er3(trial):er4(trial)),n);
                            Y9=fft(eeg(er4(trial):er(trial)),n);

                            Y14=fft(eeg(sl(trial):el1(trial)),n);
                            Y13=fft(eeg(el1(trial):el2(trial)),n);
                            Y11=fft(eeg(el2(trial):el3(trial)),n);
                            Y12=fft(eeg(el3(trial):el4(trial)),n);
                            Y10=fft(eeg(el4(trial):el(trial)),n);

                            L1=length(Y1);
                            L2=length(Y3);
                            L3=length(Y5);
                            L4=length(Y6);
                            L5=length(Y7);
                            L6=length(Y8);
                            L7=length(Y9);

                            matrix1(counter,1:length(Y1))=abs(Y1/sqrt(L1));
                            matrix2(counter,1:length(Y2))=abs(Y2/sqrt(L1));
                            matrix3(counter,1:length(Y3))=abs(Y3/sqrt(L2));
                            matrix4(counter,1:length(Y4))=abs(Y4/sqrt(L2));

                            matrix5(counter,1:length(Y5))=abs(Y5/sqrt(L3));
                            matrix6(counter,1:length(Y6))=abs(Y6/sqrt(L4));            
                            matrix7(counter,1:length(Y7))=abs(Y7/sqrt(L5));
                            matrix8(counter,1:length(Y8))=abs(Y8/sqrt(L6));
                            matrix9(counter,1:length(Y9))=abs(Y9/sqrt(L7));

                            matrix10(counter,1:length(Y10))=abs(Y10/sqrt(L3));
                            matrix11(counter,1:length(Y11))=abs(Y11/sqrt(L4));            
                            matrix12(counter,1:length(Y12))=abs(Y12/sqrt(L5));
                            matrix13(counter,1:length(Y13))=abs(Y13/sqrt(L6));
                            matrix14(counter,1:length(Y14))=abs(Y14/sqrt(L7));
                            counter=counter+1;
                        end
                    end
                end
            end
        end
    end

    matrix{1}=[matrix1;matrix2]; %quiet
    matrix{2}=[matrix3;matrix4]; %waiting
    matrix{3}=[matrix5;matrix10]; %v1
    matrix{4}=[matrix6;matrix11]; %v2
    matrix{5}=[matrix7;matrix12]; %v3
    matrix{6}=[matrix8;matrix13]; %v4
    matrix{7}=[matrix9;matrix14]; %v5

    L=[L1,L2,L3,L4,L5,L6,L7];
    c = colormap(jet(7))
    figure(a)
    for state=1:length(L)
        P2 = mean(matrix{state});
        P1 = P2(1:L(state)/2+1);
        P1(2:end-1)= 2*P1(2:end-1);
        f = Fs*(0:(L(state)/2))/L(state);
        hold on
        plot(f,P1,'Color',c(state,:),'Linewidth',1.5);
        xlabel('Frequency [Hz]');
        ylabel('Mean PSD');
        axis([1 150 0 2.7])
    end
    xlabel('Frequency [Hz]','FontSize',20);
    ylabel('Mean PSD','FontSize',20);    
    box on
    axes('Position',[.5 0.6 .3 .3])
    box on
    for state=1:7
        [n1 m1]=size(matrix{state});
        P2 = mean(matrix{state});
        vP2=sqrt(var(matrix{state}));
        mP2=vP2/sqrt(n1);
        P1 = P2(1:L(state)/2+1);
        vP1=P2(1:L(state)/2+1);
        P1(2:end-1)= 2*P1(2:end-1);
        f = Fs*(0:(L(state)/2))/L(state);
        plot(f,P1,'Color',c(state,:),'Linewidth',2);
        hold on
        xlabel('Frequency [Hz]','FontSize',10);
        ylabel('Mean PSD','FontSize',10);
        axis([1 15 0 2.7])
    end
    axes('Position',[.5 0.2 .3 .3])
    box on
    for state=1:7
        [n1 m1]=size(matrix{state});
        P2 = mean(matrix{state});
        vP2=sqrt(var(matrix{state}));
        mP2=vP2/sqrt(n1);
        P1 = P2(1:L(state)/2+1);
        vP1=P2(1:L(state)/2+1);
        P1(2:end-1)= 2*P1(2:end-1);
        f = Fs*(0:(L(state)/2))/L(state);
        plot(f,P1,'Color',c(state,:),'Linewidth',2);
        hold on
        xlabel('Frequency [Hz]','FontSize',10);
        ylabel('Mean PSD','FontSize',10);
        axis([40 120 0.05 0.15])
    end
    nombre=strcat('C:\Users\meli__000\Desktop\Melisa\Doctorado\Spectrogram\Figuras2\','speeds_',int2str(rat),'_',name(a,:),'_',prot);
    saveas(a,nombre,'png')
    close(a)
end

