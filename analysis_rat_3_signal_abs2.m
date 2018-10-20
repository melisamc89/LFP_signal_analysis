clear all
load('EC1_3','-mat')
%load('EC1_3_gamma','-mat')
load('EC1_3_vx','-mat')
load('EC1_3_signal','-mat')
load('EC1_3_amplitude','-mat')

counter=0;
for index=1:length(signal_low)
    for trial=1:length(signal_low{index})
        counter=counter+1;
        matrix_low(:,counter)=signal_low{index}{trial};
        matrix_high(:,counter)=signal_high{index}{trial};
        matrix_signal(:,counter)=signal{index}{trial};
        matrix_pos(:,counter)=signal_pos{index}{trial};
        matrix_amplitude_low(:,counter)=amp_low{index}{trial};
        matrix_amplitude_high(:,counter)=amp_high{index}{trial};        
    end
end

n_fs=400;
init=[1,2,7,12,42,47,48,54,84,89,94]*n_fs;

for i=1:length(init)-1
    signal1=[];
    signal2=[];
    signal3=[];
    signal4=[];
    for j=1:size(matrix_low,2)
        x=matrix_low(init(i):init(i+1)-1,j);
        y=matrix_low(init(i)+1:init(i+1),j);
        index=find(x<0 & y>0);
        index=init(i)+index;
        if isempty(index)~=1
        signal1=[signal1 matrix_low(index(1):index(end),j)'];
        signal2=[signal2 matrix_high(index(1):index(end),j)']; 
        signal3=[signal3 matrix_amplitude_low(index(1):index(end),j)']; 
        signal4=[signal4 matrix_amplitude_high(index(1):index(end),j)']; 
        end
    end
    time_signal_phase{i}=signal1;
    time_signal_amp{i}=signal2;
    time_signal_amp_low{i}=signal3;
    time_signal_amp_high{i}=signal4;
end

for i=1:length(time_signal_phase)
    subplot(2,5,i)
    histogram(time_signal_amp_low{i},'Binwidth',0.1,'Normalization','probability')
    hold on
    histogram(time_signal_amp_high{i},'Binwidth',0.1,'Normalization','probability')
    xlim([0 5])
    mean_low_amp(i)=mean(time_signal_amp_low{i});
end

clear time_signal_phase time_signal_amp
clear time_signal_amp_low time_signal_amp_high
n_fs=400;
init=[1,2,7,12,42,47,48,54,84,89,94]*n_fs;

for i=1:length(init)-1
    signal1=[];
    signal2=[];
    signal3=[];
    signal4=[];
    signal5=[];
    signal6=[];
    signal7=[];
    signal8=[];
    for j=1:size(matrix_low,2)
        x=matrix_low(init(i):init(i+1)-1,j);
        y=matrix_low(init(i)+1:init(i+1),j);
        index=find(x<0 & y>0);
        index=init(i)+index;
        if isempty(index)~=1
            if mean(matrix_amplitude_low(index(1):index(end),j))<mean_low_amp(i)
                signal1=[signal1 matrix_low(index(1):index(end),j)'];
                signal2=[signal2 matrix_high(index(1):index(end),j)']; 
                signal3=[signal3 matrix_amplitude_low(index(1):index(end),j)']; 
                signal4=[signal4 matrix_amplitude_high(index(1):index(end),j)']; 
            else
                signal5=[signal5 matrix_low(index(1):index(end),j)'];
                signal6=[signal6 matrix_high(index(1):index(end),j)']; 
                signal7=[signal7 matrix_amplitude_low(index(1):index(end),j)']; 
                signal8=[signal8 matrix_amplitude_high(index(1):index(end),j)']; 
            end
        end
    end
    time_signal_phase{i}{1}=signal1;
    time_signal_amp{i}{1}=signal2;
    time_signal_amp_low{i}{1}=signal3;
    time_signal_amp_high{i}{1}=signal4;
    time_signal_phase{i}{2}=signal5;
    time_signal_amp{i}{2}=signal6;
    time_signal_amp_low{i}{2}=signal7;
    time_signal_amp_high{i}{2}=signal8;
end

for i=1:length(time_signal_phase)
    for j=1:2
     [phase amplitude]=MakeMIHistogram(time_signal_phase{i}{j},time_signal_amp{i}{j},20);
     h_low(i,j)=mean(time_signal_amp_low{i}{j});
     h_low_std(i,j)=var(time_signal_amp_low{i}{j});
     h_high(i,j)=mean(time_signal_amp_high{i}{j});
     h_high_std(i,j)=var(time_signal_amp_low{i}{j});     
     mi(i,j)=ModularityIndex(phase,amplitude);
    end
end

speeds=[0,0,36,6,0,0,0,-6,-36,0];
plot(speeds,mi(:,1),'k')
hold on
errorbar(speeds,h_low(:,1)*0.01,h_low_std(:,1)*0.01,'r')
plot(speeds,mi(:,2),'b')
hold on
errorbar(speeds,h_low(:,2)*0.01,h_high_std(:,2)*0.01,'r')

plot(mi(:,1),'k')
hold on
errorbar(h_low(:,1)*0.01,h_low_std(:,1)*0.01,'r')
plot(mi(:,2),'b')
hold on
errorbar(h_low(:,2)*0.01,h_high_std(:,2)*0.01,'r')

dots1=[0.1,0.1,36,6,0.1];
dots2=[-0.1,-0.1,-6,-36,-0.1];
c=colormap(jet(10));
for i=1:2
    subplot(2,1,i)
plot([dots1(1),dots2(1)],[mi(1,i),mi(6,i)],'o','Color',c(1,:));
hold on
plot([dots1(1),dots2(1)],[h_low(1,i),h_low(6,i)]*0.01,'o','Color',c(2,:));
plot([dots1(2),dots2(2)],[mi(2,i),mi(7,i)],'o','Color',c(3,:));
plot([dots1(2),dots2(2)],[h_low(2,i),h_low(7,i)]*0.01,'o','Color',c(4,:));
plot([dots1(3),dots2(4)],[mi(3,i),mi(9,i)],'o','Color',c(5,:));
plot([dots1(3),dots2(4)],[h_low(3,i),h_low(9,i)]*0.01,'o','Color',c(6,:));
plot([dots1(4),dots2(3)],[mi(4,i),mi(8,i)],'o','Color',c(7,:));
plot([dots1(4),dots2(3)],[h_low(4,i),h_low(8,i)]*0.01,'o','Color',c(8,:));
plot([dots1(5),dots2(5)],[mi(5,i),mi(10,i)],'o','Color',c(9,:));
plot([dots1(5),dots2(5)],[h_low(5,i),h_low(10,i)]*0.01,'o','Color',c(10,:));
xlabel('Speed (cm/s)')
ylabel('MI')
box on
%legend('BeforeAlarm','AfterAlarm','FastRunning','SlowRunning','AfterRunning')
end

%% Pooling analysis

clear all
load('EC1_3','-mat')
%load('EC1_3_gamma','-mat')
load('EC1_3_vx','-mat')
load('EC1_3_signal','-mat')
load('EC1_3_amplitude','-mat')

counter=0;
for index=1:length(signal_low)
    for trial=1:length(signal_low{index})
        counter=counter+1;
        matrix_low(:,counter)=signal_low{index}{trial};
        matrix_high(:,counter)=signal_high{index}{trial};
        matrix_signal(:,counter)=signal{index}{trial};
        matrix_pos(:,counter)=signal_pos{index}{trial};
        matrix_amplitude_low(:,counter)=amp_low{index}{trial};
        matrix_amplitude_high(:,counter)=amp_high{index}{trial};        
    end
end

n_fs=400;
init=[1,2,7,12,42,47,48,54,84,89,94]*n_fs;
pools=175;

for i=1:length(init)-1
   for j=1:floor(size(matrix_low,2)/pools)
        signal1=[];
        signal2=[];
        signal3=[];
        signal4=[];
        for k=1:pools
            x=matrix_low(init(i):init(i+1)-1,(j-1)*pools+k);
            y=matrix_low(init(i)+1:init(i+1),(j-1)*pools+k);
            index=find(x<0 & y>0);
            index=init(i)+index;
            if isempty(index)~=1
            signal1=[signal1 matrix_low(index(1):index(end),(j-1)*pools+k)'];
            signal2=[signal2 matrix_high(index(1):index(end),(j-1)*pools+k)']; 
            signal3=[signal3 matrix_amplitude_low(index(1):index(end),(j-1)*pools+k)']; 
            signal4=[signal4 matrix_amplitude_high(index(1):index(end),(j-1)*pools+k)']; 
            end
        end
        time_signal_phase{i}{j}=signal1;
        time_signal_amp{i}{j}=signal2;
        time_low_amo{i}{j}=signal3;
        time_high_amo{i}{j}=signal4;
   end
end


for j=1:floor(size(matrix_low,2)/pools)
    new_time{1}{j}=[time_low_amo{1}{j},time_low_amo{6}{j}];
    new_time{2}{j}=[time_low_amo{2}{j},time_low_amo{7}{j}];
    new_time{3}{j}=[time_low_amo{3}{j},time_low_amo{9}{j}];
    new_time{4}{j}=[time_low_amo{4}{j},time_low_amo{8}{j}];
    new_time{5}{j}=[time_low_amo{5}{j},time_low_amo{10}{j}];    
end

for j=1:floor(size(matrix_low,2)/pools)
    new_time_h{1}{j}=[time_high_amo{1}{j},time_high_amo{6}{j}];
    new_time_h{2}{j}=[time_high_amo{2}{j},time_high_amo{7}{j}];
    new_time_h{3}{j}=[time_high_amo{3}{j},time_high_amo{9}{j}];
    new_time_h{4}{j}=[time_high_amo{4}{j},time_high_amo{8}{j}];
    new_time_h{5}{j}=[time_high_amo{5}{j},time_high_amo{10}{j}];    
end

for i=1:length(new_time)
    for j=1:3
        h_new(i,j)=mean(new_time{i}{j});
        %h_new_h(i,j)=mean(new_time_h{i}{j});
    end
end

for i=1:5
    subplot(1,5,i)
    histogram(new_time{i}{1},'Binwidth',0.1,'Normalization','probability')
    hold on
    histogram(new_time{i}{2},'Binwidth',0.1,'Normalization','probability')
    histogram(new_time{i}{3},'Binwidth',0.1,'Normalization','probability')
    xlim([0 5])
    ylim([0 0.14])
end


for i=1:length(time_signal_phase)
    for j=1:5
        [phase amplitude]=MakeMIHistogram(time_signal_phase{i}{j},time_signal_amp{i}{j},20);
        h_low(i,j)=mean(time_low_amo{i}{j});
        mi(i,j)=ModularityIndex(phase,amplitude);
    end
end

color=colormap(jet(3));
for i=1:3
    plot(h_low(:,i),'o','Color',color(i,:))
    hold on
end
