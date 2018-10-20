
%rat=3
n_fs=400;
init=[1,2,7,12,42,47,48,54,84,89,94]*n_fs;

number=[9,8,7,2,4,3];
for i=1:length(number)
    x=[];
    y=[];
    for j=1:size(matrix_low,2)
        index=find(matrix_low(init(number(i)):init(number(i)+1)-1,j)<0 & matrix_low(init(number(i))+1:init(number(i)+1),j)>0);
        index=init(number(i))+index;
        if isempty(index)~=1
            x=[x,matrix_low(index(1):index(end),j)'];
            y=[y,matrix_high(index(1):index(end),j)'];
        end
    end
    [phase amplitude]=MakeMIHistogram(x,y,20);
    mi_all(i)=ModularityIndex(phase,amplitude);
    mi_phase_all(i,:)=phase;
    mi_amplitude_all(i,:)=amplitude;
end

for i=1:6
    x=[];
    y=[];
    for j=1:size(matrix_low,2)
        index=find(matrix_low(init(number(i)):init(number(i)+1)-1,j)<0 & matrix_low(init(number(i))+1:init(number(i)+1),j)>0);
        index=init(number(i))+index;
        if isempty(index)~=1
            x=[x,matrix_low(index(1):index(end),j)'];
            y=[y,matrix_high(index(1):index(end),j)'];
        end
    end
    for k=1:100
        position=round(rand(1,1)*length(x));
        new_position=length(x)-position;
        auxiliar=x;
        new_signal1(1:new_position)=auxiliar(position+1:length(x));
        new_signal1(new_position+1:length(x))=auxiliar(1:position);
        
        position=round(rand(1,1)*length(x));
        new_position=length(x)-position;
        auxiliar=y;
        new_signal2(1:new_position)=auxiliar(position+1:length(x));
        new_signal2(new_position+1:length(x))=auxiliar(1:position);
        
        [phase amplitude]=MakeMIHistogram(new_signal1,new_signal2,20);
        mi_test(i,k)=ModularityIndex(phase,amplitude);
     end
    clear new_signal1 new_signal2
end


for i=1:6
    auxiliar=sort(mi_test(i,:));
    surrogate(i)=auxiliar(100);
end

subplot(2,2,1)
speeds=[-36,-6,-0.5,0.5,6,36];
bar(speeds,mi_all)
hold on
bar(speeds,surrogate)

subplot(2,2,2)
histogram(mi_test)

clear all
%rat exponential RAT14566
load('EC1_1','-mat')
%load('EC1_1_gamma','-mat')
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

number=[15,14,13,12,11,10,2,3,4,5,6,7];
for i=1:length(number)
    x=[];
    y=[];
    for j=1:size(matrix_low,2)
        index=find(matrix_low(init(number(i)):init(number(i)+1)-1,j)<0 & matrix_low(init(number(i))+1:init(number(i)+1),j)>0);
        index=init(number(i))+index;
        if isempty(index)~=1
            x=[x,matrix_low(index(1):index(end),j)'];
            y=[y,matrix_high(index(1):index(end),j)'];
        end
    end
    [phase amplitude]=MakeMIHistogram(x,y,20);
    mi_all(i)=ModularityIndex(phase,amplitude);
    mi_phase_all(i,:)=phase;
    mi_amplitude_all(i,:)=amplitude;
end


for i=1:12
    x=[];
    y=[];
    for j=1:size(matrix_low,2)
        index=find(matrix_low(init(number(i)):init(number(i)+1)-1,j)<0 & matrix_low(init(number(i))+1:init(number(i)+1),j)>0);
        index=init(number(i))+index;
        if isempty(index)~=1
            x=[x,matrix_low(index(1):index(end),j)'];
            y=[y,matrix_high(index(1):index(end),j)'];
        end
    end
    for k=1:100
        position=round(rand(1,1)*length(x));
        new_position=length(x)-position;
        auxiliar=x;
        new_signal1(1:new_position)=auxiliar(position+1:length(x));
        new_signal1(new_position+1:length(x))=auxiliar(1:position);
        
        position=round(rand(1,1)*length(x));
        new_position=length(x)-position;
        auxiliar=y;
        new_signal2(1:new_position)=auxiliar(position+1:length(x));
        new_signal2(new_position+1:length(x))=auxiliar(1:position);
        
        [phase amplitude]=MakeMIHistogram(new_signal1,new_signal2,20);
        mi_test(i,k)=ModularityIndex(phase,amplitude);
     end
    clear new_signal1 new_signal2
end


for i=1:12
    auxiliar=sort(mi_test(i,:));
    surrogate(i)=auxiliar(100);
end

subplot(2,2,3)
speeds=[-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5];
bar(speeds,mi_all)
hold on
bar(speeds,surrogate)

subplot(2,2,4)
histogram(mi_test)