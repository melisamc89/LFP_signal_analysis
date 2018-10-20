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

pools=320;
pool_trial=1:pools:size(matrix_low,2);

n_fs=400;
bin_size=2*n_fs;
overlap=bin_size/4;
windows=floor(size(matrix_low,1)/bin_size);
init=1:overlap:size(matrix_low,1);

for i=1:length(init)-2
    for j=1:floor(size(matrix_low,2)/pools)
        signal1=[];
        signal2=[];
        signal3=[];
        signal4=[];
        for k=1:pools
            x=matrix_low(init(i):init(i)+bin_size-1,(j-1)*pools+k);
            y=matrix_low(init(i)+1:init(i)+bin_size,(j-1)*pools+k);
            index=find(x<0 & y>0);
            index=init(i)+index;
            if isempty(index)~=1
            signal1=[signal1 matrix_low(index(1):index(end),(j-1)*pools+k)'];
            signal2=[signal2 matrix_high(index(1):index(end),(j-1)*pools+k)']; 
            signal3=[signal3 matrix_amplitude_low(index(1):index(end),(j-1)*pools+k)'];
            signal4=[signal4 matrix_amplitude_high(index(1):index(end),(j-1)*pools+k)']; 
            end
            time_signal_phase{i}{j}=signal1;
            time_signal_amp{i}{j}=signal2;
            time_signal_amp_low{i}{j}=signal3;
            time_signal_amp_high{i}{j}=signal4;
        end
    end
end

for i=1:length(time_signal_amp_low)
    for j=1:floor(size(matrix_low,2)/pools)
        h_low(i,j)=mean(time_signal_amp_low{i}{j});
        [phase amplitude]=MakeMIHistogram(time_signal_phase{i}{j},time_signal_amp{i}{j},20);
        mi(i,j)=ModularityIndex(phase,amplitude);
    end
end

for j=1:floor(size(matrix_low,2)/pools)
    h_low_mean(j)=mean(h_low(:,j));
    h_low_std(j)=var(h_low(:,j));
end
errorbar(h_low_mean(1:end),h_low_std(1:end))

color=[0,0,1;0,1,0;1,0,0];
tiempo=0:0.5:length(h_low)*0.5-0.5;

for i=1:floor(size(matrix_low,2)/pools)
    subplot(3,1,i)
    plot(tiempo,h_low(:,i)/mean(h_low(:,i)),'--','Linewidth',2,'Color',color(i,:))
    hold on
    plot(tiempo,mi(:,i)/mean(mi(:,i)),'-','Linewidth',2,'Color',color(i,:))
end

for i=1:floor(size(matrix_low,2)/pools)
    subplot(3,1,i)
    histogram(mi(:,i),'Binwidth',0.0001,'Normalization','probability')
    hold on
    xlim([0 0.003])
    ylim([0 0.2])
end

cat1{1}=[];cat1{2}=[];cat1{3}=[];
for i=1:121
    cat1{1}=[cat1{1},time_signal_amp_low{i}{1}];
    cat1{2}=[cat1{2},time_signal_amp_low{i}{3}];
    cat1{3}=[cat1{3},time_signal_amp_low{i}{3}];
end

for i=1:floor(size(matrix_low,2)/pools)
    subplot(3,1,i)
    histogram(cat1{i},'Binwidth',0.01,'Normalization','probability')
    hold on
    %xlim([0.5 1.5])
    %ylim([0 0.1])
end

