clear all
load('EC1_1','-mat')
load('EC1_1_vx','-mat')
load('EC1_1_amplitude','-mat')

load('EC1_3','-mat')
load('EC1_3_vx','-mat')
load('EC1_3_amplitude','-mat')

counter=0;
for index=1:length(signal_low)
    for trial=1:length(signal_low{index})
        counter=counter+1;
        matrix_low(:,counter)=signal_low{index}{trial};
        matrix_high(:,counter)=signal_high{index}{trial};
        matrix_amplitude_low(:,counter)=amp_low{index}{trial};
        matrix_pos(:,counter)=signal_pos{index}{trial};
    end
end

%rat=3
n_fs=400;
init=[1,2,7,12,42,47,48,54,84,89,94]*n_fs;

%rat=1
n_fs=400;
init=[1,2,7,11,15,19,23,27,31,32,38,42,46,50,54,58,62]*n_fs;



for i=1:length(init)-1
    signal1{i}=matrix_low(init(i):init(i+1)-1,:);
    signal2{i}=matrix_high(init(i):init(i+1)-1,:);
    signal3{i}=matrix_amplitude_low(init(i):init(i+1)-1,:);
end


for i=1:length(signal1)
    for j=1:size(signal1{i},2)
        [phase amplitude]=MakeMIHistogram(signal1{i}(:,j),signal2{i}(:,j),20);
        amp(i,j,:)=amplitude;
        mi(i,j)=ModularityIndex(phase,amplitude);
    end
end

color=colormap(jet(floor(size(signal1{1},2)/40)));
for i=1:length(signal1)
    subplot(2,8,i)
    for j=1:length(color)
        value=reshape(mean(amp(i,(j-1)*40+1:j*40,:)),[20 1]);
        plot([phase-2*pi,phase,phase+2*pi],[value;value;value],'color',color(j,:))
        hold on
        mi_mean(i,j)=ModularityIndex(phase,value);
    end
    xlim([-3*pi/2 3*pi/2])
    ylim([0.035 0.065])
end

figure(2)
color=colormap(jet(floor(size(signal1{1},2)/40)));
for i=1:length(signal1)
    subplot(2,8,i)
    scatter(1:length(mi_mean(i,:)),mi_mean(i,:),'filled')
    xlim([0 30])
    ylim([0 0.0035])    
end

speeds=[0,0,36,6,0,0,0,-6,-36,0];
matrix=[speeds;h_low;h_high;mi*1000];
matrix_error=[zeros(1,length(speeds));h_low_std;h_high_std;zeros(1,length(speeds))];
barwitherr(matrix_error,matrix);
for i=1:4
    subplot(2,2,i)
    barwitherr(matrix_error(i,:)',matrix(i,:)')
end
