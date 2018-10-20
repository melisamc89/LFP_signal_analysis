clear all
load('EC1_1','-mat')
load('EC1_1_vx','-mat')
%load('EC1_1_signal','-mat')

counter=0;
for index=1:length(signal_low)
    for trial=1:length(signal_low{index})
        counter=counter+1;
        matrix_low(:,counter)=signal_low{index}{trial};
        matrix_high(:,counter)=signal_high{index}{trial};
        %matrix_signal(:,counter)=signal{index}{trial};
        matrix_pos(:,counter)=signal_pos{index}{trial};
    end
end

%% Analysis of MI in pools of days
clear mi time_signal_phase time_signal_amp
pools=100;
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
        for k=1:pools
            x=matrix_low(init(i):init(i)+bin_size-1,(j-1)*pools+k);
            y=matrix_low(init(i)+1:init(i)+bin_size,(j-1)*pools+k);
            index=find(x<0 & y>0);
            index=init(i)+index;
            if isempty(index)~=1
            signal1=[signal1 matrix_low(index(1):index(end),(j-1)*pools+k)'];
            signal2=[signal2 matrix_high(index(1):index(end),(j-1)*pools+k)']; 
            end
            time_signal_phase{i}{j}=signal1;
            time_signal_amp{i}{j}=signal2;
        end
    end
end

for i=1:length(time_signal_phase)
    for j=1:7 
        h_low(i,j)=mean(time_signal_amp{i}{j});
        [phase amplitude]=MakeMIHistogram(time_signal_phase{i}{j},time_signal_amp{i}{j},20);
        mi(i,j)=ModularityIndex(phase,amplitude);
    end
end
  
x=resample(mean(matrix_pos'),1,25);
y(1:floor(length(mi)/2))=x(1:floor(length(mi)/2));
y(floor(length(mi)/2)+1:length(mi))=-x(floor(length(mi)/2)+1:length(mi));

plot(0:0.5:length(mi)*0.5-0.5,y(1:length(mi))*0.0001,'k','Linewidth',2)
hold on
color=colormap(jet(4));
for i=1:4
plot(0:0.5:length(mi)*0.5-0.5,mi(:,i),'Linewidth',2,'Color',color(i,:))
hold on
end

%% Speed catting
clear mi time_signal_phase time_signal_amp
n_fs=400;
init=[1,2,7,11,15,19,23,27,31,32,38,42,46,50,54,58,62]*n_fs;
for i=1:length(init)-1
   for j=1:floor(size(matrix_low,2)/pools)
        signal1=[];
        signal2=[];
        for k=1:pools
            x=matrix_low(init(i):init(i+1)-1,(j-1)*pools+k);
            y=matrix_low(init(i)+1:init(i+1),(j-1)*pools+k);
            index=find(x<0 & y>0);
            index=init(i)+index;
            if isempty(index)~=1
            signal1=[signal1 matrix_low(index(1):index(end),(j-1)*pools+k)'];
            signal2=[signal2 matrix_high(index(1):index(end),(j-1)*pools+k)']; 
            end
        end
        time_signal_phase{i}{j}=signal1;
        time_signal_amp{i}{j}=signal2;
   end
end


for i=1:length(time_signal_phase)
    for j=1:4
        [phase amplitude]=MakeMIHistogram(time_signal_phase{i}{j},time_signal_amp{i}{j},20);
        mi(i,j)=ModularityIndex(phase,amplitude);
    end
end


