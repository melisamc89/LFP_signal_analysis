clear all
n_fs=400;
fs=4800;
[16262,16347,16456,16653];
count=ones(5,1);
count2=ones(10,1);
for rat=[16262,16347,16456,16653];
    files=Data_Listing2(rat);
    abs_speed_vector=[0,20,40,60,80];
    speed_vector=[0,20,40,60,80,0,-20,-40,-60,-80];
    [n1 n2]=size(files);
    for index=1:n1
        load(['RandomProtocol/' files(index,:)],'-mat');
        eeg_vector=resample(db.V.ecL.eeg,1,floor(fs/n_fs));
        n_fs=fs/floor(fs/n_fs);
        [qual4,qual,qual2]=egf_theta_quality_2016(db.V.ecL.eeg);
        if qual>0.5
             [phase1 amplitude1]=Hilbert(eeg_vector,1.5,4.5,n_fs);
             [phase2 amplitude2]=Hilbert(eeg_vector,6,12,n_fs);
             
             amp1=amplitude1/mean(amplitude1);
             amp2=amplitude2/mean(amplitude2);
             %[phase1 amplitude1]=Hilbert(eeg_vector,6,12,n_fs);
             %[phase2 amplitude2]=Hilbert(eeg_vector,40,80,n_fs);
        
            for i=1:length(db.V.pos.log.alarm_r)
                time((i-1)*2+1,1)=floor(db.V.pos.log.alarm_r(i)*n_fs);
                time(i*2,:)=floor(db.V.pos.log.alarm_l(i)*n_fs);
            end
            
            random_t=reshape(db.V.pos.log.speed.t(2:end),...
                [4 floor((length(db.V.pos.log.speed.t)-1)/4)])'/1000*n_fs;
            
            random_time=floor([time,random_t]);
            random_speed=reshape(db.V.pos.log.speed.log_speed(1:4*floor(length(db.V.pos.log.speed.t)/4)),...
                [4 floor(length(db.V.pos.log.speed.t)/4)])';
        
            %random_speed=[random_speed1(:,2),random_speed1(:,3)];
            %random_time=[random_time1(:,2),random_time1(:,3),random_time1(:,4)];
            
            for i=1:size(random_speed,1)
                for j=2:3               
                    for k=1:length(abs_speed_vector)
                        if abs(random_speed(i,j))==abs_speed_vector(k)
                            ecL_phase{k}{count(k)}=phase1(random_time(i,j):random_time(i,j+1));
                            ecL_amplitude{k}{count(k)}=amplitude2(random_time(i,j):random_time(i,j+1));
                            ecR_amplitude_low{k}{count(k)}=amp1(random_time(i,j):random_time(i,j+1));
                            ecR_amplitude_high{k}{count(k)}=amp2(random_time(i,j):random_time(i,j+1));                            
                            count(k)=count(k)+1;
                        end
                    end
                    for k=1:length(speed_vector)
                        if random_speed(i,j)==speed_vector(k)
                            ecL_phase2{k}{count2(k)}=phase1(random_time(i,j):random_time(i,j+1));
                            ecL_amplitude2{k}{count2(k)}=amplitude2(random_time(i,j):random_time(i,j+1));
                            ecR_amplitude_low2{k}{count2(k)}=amp1(random_time(i,j):random_time(i,j+1));
                            ecR_amplitude_high2{k}{count2(k)}=amp2(random_time(i,j):random_time(i,j+1));
                            count2(k)=count2(k)+1;
                        end
                    end                    
                end
            end
        end
        
    end
end
