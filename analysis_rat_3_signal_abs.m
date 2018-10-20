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
bin_size=2*n_fs;
overlap=bin_size/4;
windows=floor(size(matrix_low,1)/bin_size);
init=1:overlap:size(matrix_low,1);

for i=1:length(init)-2
    signal1=[];
    signal2=[];
    signal3=[];
    signal4=[];
    for j=1:size(matrix_low,2)
            x=matrix_low(init(i):init(i)+bin_size-1,j);
            y=matrix_low(init(i)+1:init(i)+bin_size,j);
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
     [phase amplitude]=MakeMIHistogram(time_signal_phase{i},time_signal_amp{i},20);
     h_low(i)=mean(time_signal_amp_low{i});
     h_high(i)=mean(time_signal_amp_high{i});
     mi(i)=ModularityIndex(phase,amplitude);
end
  
x=resample(mean(matrix_pos'),1,25);
y(1:floor(length(mi)/2))=x(1:floor(length(mi)/2));
y(floor(length(mi)/2)+1:length(mi))=-x(floor(length(mi)/2)+1:length(mi));

plot(0:0.5:length(mi)*0.5-0.5,mi,'Linewidth',2,'Color','g')
hold on
plot(0:0.5:length(mi)*0.5-0.5,y(1:length(mi))*0.0001,'k','Linewidth',2)
plot(0:0.5:length(mi)*0.5-0.5,h_low(1:length(mi))*0.005,'r','Linewidth',2)
plot(0:0.5:length(mi)*0.5-0.5,h_high(1:length(mi))*0.005,'b','Linewidth',2)




