clear all
load('EC1_1','-mat')
load('EC1_1_vx','-mat')
load('EC1_1_amplitude','-mat')

counter=0;
for index=1:length(signal_low)
    for trial=1:length(signal_low{index})
        counter=counter+1;
        matrix_low(:,counter)=signal_low{index}{trial};
        matrix_high(:,counter)=signal_high{index}{trial};
        matrix_amplitude_low(:,counter)=amp_low{index}{trial};
        matrix_amplitude_high(:,counter)=amp_high{index}{trial};
        matrix_pos(:,counter)=signal_pos{index}{trial};
    end
end

n_fs=400;
init=[1,2,7,11,15,19,23,27,31,32,38,42,46,50,54,58,62]*n_fs;

pools=100;
pool_trial=1:pools:size(matrix_low,2);

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
        time_signal_amp_low{i}{j}=signal3;
        time_signal_amp_high{i}{j}=signal4;
   end
end

for j=1:floor(size(matrix_low,2)/pools)
    new_time{1}{j}=[time_signal_amp_low{1}{j},time_signal_amp_low{9}{j}];
    new_time{2}{j}=[time_signal_amp_low{2}{j},time_signal_amp_low{10}{j}];
    new_time{3}{j}=[time_signal_amp_low{3}{j},time_signal_amp_low{15}{j}];
    new_time{4}{j}=[time_signal_amp_low{4}{j},time_signal_amp_low{14}{j}];
    new_time{5}{j}=[time_signal_amp_low{5}{j},time_signal_amp_low{13}{j}];
    new_time{6}{j}=[time_signal_amp_low{6}{j},time_signal_amp_low{12}{j}];
    new_time{7}{j}=[time_signal_amp_low{7}{j},time_signal_amp_low{11}{j}];
    new_time{8}{j}=[time_signal_amp_low{8}{j},time_signal_amp_low{16}{j}];
end

for i=1:length(new_time)
    for j=1:floor(size(matrix_low,2)/pools)
        h_low_new(i,j)=mean(new_time{i}{j});
    end
end

for i=1:length(time_signal_phase)
    for j=1:floor(size(matrix_low,2)/pools)
        [phase amplitude]=MakeMIHistogram(time_signal_phase{i}{j},time_signal_amp{i}{j},20);
        mi(i,j)=ModularityIndex(phase,amplitude);
        h_low(i,j)=mean(time_signal_amp_low{i}{j});
        h_high(i,j)=mean(time_signal_amp_high{i}{j});
    end
end

color=colormap(jet(7));
for i=1:7
%plot(0:0.5:length(mi)*0.5-0.5,mi(:,i),'Linewidth',2,'Color',color(i,:))
hold on
%plot(0:0.5:length(mi)*0.5-0.5,mi(:,i),'Linewidth',2,'Color',color(i,:))
plot(0:0.5:length(mi)*0.5-0.5,h_low(1:length(mi),i),'Linewidth',2,'Color',color(i,:))
%plot(0:0.5:length(mi)*0.5-0.5,h_high(1:length(mi),i)*0.01,'b','Linewidth',2)
end
