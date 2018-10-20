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
        h_new_h(i,j)=mean(new_time_h{i}{j});
    end
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
