%% KHPD 20160520

clearvars

%% set up params for eventual CCM (relevant to how much data to grab)

E = 10;
tau = E;
maxL = 400;
Lskip = 2;

%% get data from EC100 and EC108 from a couple days each

% EC100 9/14
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC100/EC100_f9fb03e7/mat_allChannels/EC100_f9fb03e7_all_002.mat','data');
EC100_f9 = tmp.data;

% EC100 9/16
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC100/EC100_6cff1ed2/mat_allChannels/EC100_6cff1ed2_all_002.mat','data');
EC100_6c = tmp.data;

% EC100 9/17
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC100/EC100_d663e5c0/mat_allChannels/EC100_d663e5c0_all_002.mat','data');
EC100_d6 = tmp.data;

% EC100 9/18
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC100/EC100_3f450082/mat_allChannels/EC100_3f450082_all_002.mat','data');
EC100_3f = tmp.data;

% EC100 9/19
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC100/EC100_1f1a51b1/mat_allChannels/EC100_1f1a51b1_all_002.mat','data');
EC100_1f = tmp.data;

% EC100 9/20
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC100/EC100_bd59e758/mat_allChannels/EC100_bd59e758_all_002.mat','data');
EC100_bd1 = tmp.data;

% EC100 9/21
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC100/EC100_bd59e758/mat_allChannels/EC100_bd59e758_all_090.mat','data');
EC100_bd2 = tmp.data;


%%%%%%%%%%%%%%%%%%%%

% EC108 12/8
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC108/EC108_5408ebc9/mat_allChannels/EC108_5408ebc9_all_002.mat','data');
EC108_54 = tmp.data;

% EC108 12/9
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC108/EC108_2fd7f1ac/mat_allChannels/EC108_2fd7f1ac_all_002.mat','data');
EC108_2f = tmp.data;

% EC108 12/10
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC108/EC108_9ae22432/mat_allChannels/EC108_9ae22432_all_002.mat','smallData');
EC108_9a = tmp.smallData;

% EC108 12/11
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC108/EC108_ff348c69/mat_allChannels/EC108_ff348c69_all_002.mat','data');
EC108_ff = tmp.data;

% EC108 12/12
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC108/EC108_8337e405/mat_allChannels/EC108_8337e405_all_002.mat','data');
EC108_83 = tmp.data;

% EC108 12/13
tmp=load('/sabes_data3/kderosier/Chang_lab_data/EC108/EC108_cd686c40/mat_allChannels/EC108_cd686c40_all_002.mat','data');
EC108_cd = tmp.data;

clear tmp

%% set indices for each area for each subject

% from EC108_montage spreadsheet, accounting for dropping the language grid
% during preprocessing
EC108_idx.ofc = [65:96]-64;
EC108_idx.am = [149:158]-64;
EC108_idx.hp = [159:168]-64;
EC108_idx.acc = [179:188]-64;

% from EC100 montage .mat file, accounting for dropping the language grid
% during preprocessing
EC100_idx.ofc = [73:80]-64;
EC100_idx.am = [89:98]-64;
EC100_idx.hp = [99:108]-64;
EC100_idx.acc = [129:138]-64;


%% for each recording, for each area, look to find bad electrodes
for i=0 % hide this code since i only needed it once and it's long
% figure()
% %EC108
% hold on
% for chan = EC108_idx.ofc
%     plot(EC108_2f(chan, 1:512)+(chan-EC108_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC108_2f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.ofc
%     plot(EC108_54(chan, 1:512)+(chan-EC108_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC108_54',chan));
%     pause()
%  end
% hold off
% clf
% hold on
% for chan = EC108_idx.ofc   
%     plot(EC108_83(chan, 1:512)+(chan-EC108_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC108_83',chan));
%     pause()
%  end
% hold off
% clf
% hold on
% for chan = EC108_idx.ofc   
%     plot(EC108_9a(chan, 1:512)+(chan-EC108_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC108_9a',chan));
%     pause()
%   end
% hold off
% clf
% hold on
% for chan = EC108_idx.ofc  
%     plot(EC108_cd(chan, 1:512)+(chan-EC108_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC108_cd',chan));
%     pause()
%  end
% hold off
% clf
% hold on
% for chan = EC108_idx.ofc   
%     plot(EC108_ff(chan, 1:512)+(chan-EC108_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC108_ff',chan));
%     pause()
% end
% hold off
% clf
% 
% 
% hold on
% for chan = EC108_idx.am
%     plot(EC108_2f(chan, 1:512)+(chan-EC108_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC108_2f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.am
%     plot(EC108_54(chan, 1:512)+(chan-EC108_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC108_54',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.am
%     plot(EC108_83(chan, 1:512)+(chan-EC108_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC108_83',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.am    
%     plot(EC108_9a(chan, 1:512)+(chan-EC108_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC108_9a',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.am    
%     plot(EC108_cd(chan, 1:512)+(chan-EC108_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC108_cd',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.am    
%     plot(EC108_ff(chan, 1:512)+(chan-EC108_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC108_ff',chan));
%     pause()
% end
% hold off
% clf
% 
% hold on
% for chan = EC108_idx.hp
%     plot(EC108_2f(chan, 1:512)+(chan-EC108_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC108_2f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.hp
%     plot(EC108_54(chan, 1:512)+(chan-EC108_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC108_54',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.hp    
%     plot(EC108_83(chan, 1:512)+(chan-EC108_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC108_83',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.hp    
%     plot(EC108_9a(chan, 1:512)+(chan-EC108_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC108_9a',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.hp    
%     plot(EC108_cd(chan, 1:512)+(chan-EC108_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC108_cd',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.hp    
%     plot(EC108_ff(chan, 1:512)+(chan-EC108_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC108_ff',chan));
%     pause()
%  
% end
% hold off
% clf
% 
% hold on
% for chan = EC108_idx.acc
%     plot(EC108_2f(chan, 1:512)+(chan-EC108_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC108_2f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.acc
%     plot(EC108_54(chan, 1:512)+(chan-EC108_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC108_54',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.acc    
%     plot(EC108_83(chan, 1:512)+(chan-EC108_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC108_83',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.acc    
%     plot(EC108_9a(chan, 1:512)+(chan-EC108_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC108_9a',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.acc    
%     plot(EC108_cd(chan, 1:512)+(chan-EC108_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC108_cd',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC108_idx.acc    
%     plot(EC108_ff(chan, 1:512)+(chan-EC108_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC108_ff',chan));
%     pause()
%  
% end
% hold off
% clf
% 
% 
% % EC100
% hold on
% for chan = EC100_idx.ofc
%     plot(EC100_1f(chan, 1:512)+(chan-EC100_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC100_1f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.ofc
%     plot(EC100_3f(chan, 1:512)+(chan-EC100_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC100_3f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.ofc    
%     plot(EC100_6c(chan, 1:512)+(chan-EC100_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC100_6c',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.ofc    
%     plot(EC100_bd1(chan, 1:512)+(chan-EC100_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC100_bd1',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.ofc    
%     plot(EC100_bd2(chan, 1:512)+(chan-EC100_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC100_bd2',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.ofc    
%     plot(EC100_d6(chan, 1:512)+(chan-EC100_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC100_d6',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.ofc    
%     plot(EC100_f9(chan, 1:512)+(chan-EC100_idx.ofc(1))*100)
%     title(sprintf('%s ofc channel %d','EC100_f9',chan));
%     pause()
%  
% end
% hold off
% clf
% 
% hold on
% for chan = EC100_idx.am
%     plot(EC100_1f(chan, 1:512)+(chan-EC100_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC100_1f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.am
%     plot(EC100_3f(chan, 1:512)+(chan-EC100_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC100_3f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.am    
%     plot(EC100_6c(chan, 1:512)+(chan-EC100_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC100_6c',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.am    
%     plot(EC100_bd1(chan, 1:512)+(chan-EC100_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC100_bd1',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.am    
%     plot(EC100_bd2(chan, 1:512)+(chan-EC100_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC100_bd2',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.am    
%     plot(EC100_d6(chan, 1:512)+(chan-EC100_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC100_d6',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.am    
%     plot(EC100_f9(chan, 1:512)+(chan-EC100_idx.am(1))*100)
%     title(sprintf('%s am channel %d','EC100_f9',chan));
%     pause()
%  
% end
% hold off
% clf
% 
% hold on
% for chan = EC100_idx.hp
%     plot(EC100_1f(chan, 1:512)+(chan-EC100_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC100_1f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.hp
%     plot(EC100_3f(chan, 1:512)+(chan-EC100_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC100_3f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.hp    
%     plot(EC100_6c(chan, 1:512)+(chan-EC100_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC100_6c',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.hp    
%     plot(EC100_bd1(chan, 1:512)+(chan-EC100_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC100_bd1',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.hp    
%     plot(EC100_bd2(chan, 1:512)+(chan-EC100_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC100_bd2',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.hp
%     plot(EC100_d6(chan, 1:512)+(chan-EC100_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC100_d6',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.hp    
%     plot(EC100_f9(chan, 1:512)+(chan-EC100_idx.hp(1))*100)
%     title(sprintf('%s hp channel %d','EC100_f9',chan));
%     pause()
%  
% end
% hold off
% clf
% 
% hold on
% for chan = EC100_idx.acc
%     plot(EC100_1f(chan, 1:512)+(chan-EC100_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC100_1f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.acc
%     plot(EC100_3f(chan, 1:512)+(chan-EC100_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC100_3f',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.acc    
%     plot(EC100_6c(chan, 1:512)+(chan-EC100_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC100_6c',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.acc    
%     plot(EC100_bd1(chan, 1:512)+(chan-EC100_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC100_bd1',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.acc    
%     plot(EC100_bd2(chan, 1:512)+(chan-EC100_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC100_bd2',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.acc    
%     plot(EC100_d6(chan, 1:512)+(chan-EC100_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC100_d6',chan));
%     pause()
% end
% hold off
% clf
% hold on
% for chan = EC100_idx.acc    
%     plot(EC100_f9(chan, 1:512)+(chan-EC100_idx.acc(1))*100)
%     title(sprintf('%s acc channel %d','EC100_f9',chan));
%     pause()
%  
% end
% hold off
% clf
% close()
% 
% clear chan
end
clear i

%% record bad electrodes

% EC 108 ofc bad none
% EC 108 am bad 93 94
% EC 108 hp bad 98 104
% EC 108 acc bad 124
% EC108_ff needs to be thrown out (ACC is noise all channels)

% EC 100 ofc bad none
% EC 100 am bad 30
% EC 100 hp bad 43, 44
% EC 100 acc bad 70 74
% EC100_1f needs to be thrown out (stim in first 2nd)
% EC100_bd2 needs to be thrown out (Am, Hp, ACC noisey all channels)

%% remove bad files and electrodes

clear EC100_1f EC100_bd2 EC108_ff

EC108_idx.am = 85:92;
EC108_idx.hp = [95:97,99:103];
EC108_idx.acc = 115:123;

EC100_idx.am = [25:29,31:34];
EC100_idx.hp = 35:42;
EC100_idx.acc = [65:69,71:73];

%% create mean signal for each area for each recording

% note that now there are conveniently 5 datasets for each patient
% so for each patient for each brain area we can organize data in 5 rows

% EC108
EC108.ofc = [mean(EC108_2f(EC108_idx.ofc,:)); ...
    mean(EC108_54(EC108_idx.ofc,:)); ...
    mean(EC108_83(EC108_idx.ofc,:)); ...
    mean(EC108_9a(EC108_idx.ofc,:)); ...
    mean(EC108_cd(EC108_idx.ofc,:)) ];

EC108.am = [mean(EC108_2f(EC108_idx.am,:)); ...
    mean(EC108_54(EC108_idx.am,:)); ...
    mean(EC108_83(EC108_idx.am,:)); ...
    mean(EC108_9a(EC108_idx.am,:)); ...
    mean(EC108_cd(EC108_idx.am,:)) ];

EC108.hp = [mean(EC108_2f(EC108_idx.hp,:)); ...
    mean(EC108_54(EC108_idx.hp,:)); ...
    mean(EC108_83(EC108_idx.hp,:)); ...
    mean(EC108_9a(EC108_idx.hp,:)); ...
    mean(EC108_cd(EC108_idx.hp,:)) ];

EC108.acc = [mean(EC108_2f(EC108_idx.acc,:)); ...
    mean(EC108_54(EC108_idx.acc,:)); ...
    mean(EC108_83(EC108_idx.acc,:)); ...
    mean(EC108_9a(EC108_idx.acc,:)); ...
    mean(EC108_cd(EC108_idx.acc,:)) ];

% EC100
EC100.ofc = [mean(EC100_3f(EC100_idx.ofc,:)); ...
    mean(EC100_6c(EC100_idx.ofc,:)); ...
    mean(EC100_bd1(EC100_idx.ofc,:)); ...
    mean(EC100_d6(EC100_idx.ofc,:)); ...
    mean(EC100_f9(EC100_idx.ofc,:)) ];

EC100.am = [mean(EC100_3f(EC100_idx.am,:)); ...
    mean(EC100_6c(EC100_idx.am,:)); ...
    mean(EC100_bd1(EC100_idx.am,:)); ...
    mean(EC100_d6(EC100_idx.am,:)); ...
    mean(EC100_f9(EC100_idx.am,:)) ];

EC100.hp = [mean(EC100_3f(EC100_idx.hp,:)); ...
    mean(EC100_6c(EC100_idx.hp,:)); ...
    mean(EC100_bd1(EC100_idx.hp,:)); ...
    mean(EC100_d6(EC100_idx.hp,:)); ...
    mean(EC100_f9(EC100_idx.hp,:)) ];

EC100.acc = [mean(EC100_3f(EC100_idx.acc,:)); ...
    mean(EC100_6c(EC100_idx.acc,:)); ...
    mean(EC100_bd1(EC100_idx.acc,:)); ...
    mean(EC100_d6(EC100_idx.acc,:)); ...
    mean(EC100_f9(EC100_idx.acc,:)) ];

%% clear the old variables

clear EC100_3f EC100_6c EC100_bd1 EC100_d6 EC100_f9 EC108_2f EC108_54 EC108_83 EC108_9a EC108_cd EC100_idx EC108_idx

%% filter Am & Hp beta signals

butterBeta = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',13,'HalfPowerFrequency2',30, ...
    'SampleRate',512);

for k=1:5
    EC100.amBeta(k,:) = filtfilt(butterBeta, EC100.am(k,:));
    EC100.hpBeta(k,:) = filtfilt(butterBeta, EC100.hp(k,:));
    
    EC108.amBeta(k,:) = filtfilt(butterBeta, EC108.am(k,:));
    EC108.hpBeta(k,:) = filtfilt(butterBeta, EC108.hp(k,:));
end

%% compute A-H beta coh for each 10s bin
cohLength = 5200; % use 5200 samples to compute coherence

cohIdx = 1:cohLength:size(EC108.ofc,2)-cohLength; % gives start idx in EC10X.area(dataset,:) for each coh bin

for d = 1:5
    
    EC100.Bcoh(d,:) = getCohNSamp( EC100.amBeta(d,:), EC100.hpBeta(d,:), 512, cohLength );
    EC108.Bcoh(d,:) = getCohNSamp( EC108.amBeta(d,:), EC108.hpBeta(d,:), 512, cohLength );
end

%% assign 400 sample chunks to be high or low

for d = 1:5
    % compute median A-H beta coh for each patient
    EC100.medCoh(d,1) = median(EC100.Bcoh(d,:));
    EC108.medCoh(d,1) = median(EC108.Bcoh(d,:));
    
    % decide high/low coh thresholds for each patient
    EC100.hiThresh(d,1) = prctile(EC100.Bcoh(d,:),95);
    EC100.loThresh(d,1) = prctile(EC100.Bcoh(d,:),5);
    
    EC108.hiThresh(d,1) = prctile(EC108.Bcoh(d,:),95);
    EC108.loThresh(d,1) = prctile(EC108.Bcoh(d,:),5);
    
end

EC100.hiO = [];
EC100.hiA = [];
EC100.hiH = [];
EC100.hiC = [];

EC100.loO = [];
EC100.loA = [];
EC100.loH = [];
EC100.loC = [];

EC108.hiO = [];
EC108.hiA = [];
EC108.hiH = [];
EC108.hiC = [];

EC108.loO = [];
EC108.loA = [];
EC108.loH = [];
EC108.loC = [];

skip = 800; % pull out 9 chunks of 400 samples each from the middle of the coh bin
l = 9*400-1;
for d=1:5
    for b = 1:length(cohIdx)
        if EC100.Bcoh(d,b) >= EC100.hiThresh(d,1) % EC100 hi
            s=cohIdx(b)+skip;
            EC100.hiO = [EC100.hiO, EC100.ofc(d, s:s+l)];
            EC100.hiA = [EC100.hiA, EC100.am(d, s:s+l)];
            EC100.hiH = [EC100.hiH, EC100.hp(d, s:s+l)];
            EC100.hiC = [EC100.hiC, EC100.acc(d, s:s+l)];
            
        elseif EC100.Bcoh(d,b) <= EC100.loThresh(d,1) % EC100 lo
            s=cohIdx(b)+skip;
            EC100.loO = [EC100.loO, EC100.ofc(d, s:s+l)];
            EC100.loA = [EC100.loA, EC100.am(d, s:s+l)];
            EC100.loH = [EC100.loH, EC100.hp(d, s:s+l)];
            EC100.loC = [EC100.loC, EC100.acc(d, s:s+l)];
            
        end
        
         if EC108.Bcoh(d,b) >= EC108.hiThresh(d,1) % EC108 hi
            s=cohIdx(b)+skip;
            EC108.hiO = [EC108.hiO, EC108.ofc(d, s:s+l)];
            EC108.hiA = [EC108.hiA, EC108.am(d, s:s+l)];
            EC108.hiH = [EC108.hiH, EC108.hp(d, s:s+l)];
            EC108.hiC = [EC108.hiC, EC108.acc(d, s:s+l)];
            
        elseif EC108.Bcoh(d,b) <= EC108.loThresh(d,1) % EC108 lo
            s=cohIdx(b)+skip;
            EC108.loO = [EC108.loO, EC108.ofc(d, s:s+l)];
            EC108.loA = [EC108.loA, EC108.am(d, s:s+l)];
            EC108.loH = [EC108.loH, EC108.hp(d, s:s+l)];
            EC108.loC = [EC108.loC, EC108.acc(d, s:s+l)];
            
        end
    end
end

%% do CCM 

% note Lvals, startIdx depend on length & params, which are same for all
n=[];

% EC100 hi
[ hi.EC100.OxmapA, hi.EC100.AxmapO, ~, ~ ] = CCM( EC100.hiO, EC100.hiA, E, tau, n, maxL, Lskip );

[ hi.EC100.OxmapH, hi.EC100.HxmapO, ~, ~ ] = CCM( EC100.hiO, EC100.hiH, E, tau, n, maxL, Lskip );

[ hi.EC100.OxmapC, hi.EC100.CxmapO, ~, ~ ] = CCM( EC100.hiO, EC100.hiC, E, tau, n, maxL, Lskip );

[ hi.EC100.CxmapA, hi.EC100.AxmapC, ~, ~ ] = CCM( EC100.hiC, EC100.hiA, E, tau, n, maxL, Lskip );

[ hi.EC100.CxmapH, hi.EC100.HxmapC, ~, ~ ] = CCM( EC100.hiC, EC100.hiH, E, tau, n, maxL, Lskip );

[ hi.EC100.AxmapH, hi.EC100.HxmapA, hi.EC100.Lvals, hi.EC100.startIdx ] = CCM( EC100.hiA, EC100.hiH, E, tau, n, maxL, Lskip );

% EC100 lo
[ lo.EC100.OxmapA, lo.EC100.AxmapO, ~, ~ ] = CCM( EC100.loO, EC100.loA, E, tau, n, maxL, Lskip );

[ lo.EC100.OxmapH, lo.EC100.HxmapO, ~, ~ ] = CCM( EC100.loO, EC100.loH, E, tau, n, maxL, Lskip );

[ lo.EC100.OxmapC, lo.EC100.CxmapO, ~, ~ ] = CCM( EC100.loO, EC100.loC, E, tau, n, maxL, Lskip );

[ lo.EC100.CxmapA, lo.EC100.AxmapC, ~, ~ ] = CCM( EC100.loC, EC100.loA, E, tau, n, maxL, Lskip );

[ lo.EC100.CxmapH, lo.EC100.HxmapC, ~, ~ ] = CCM( EC100.loC, EC100.loH, E, tau, n, maxL, Lskip );

[ lo.EC100.AxmapH, lo.EC100.HxmapA, lo.EC100.Lvals, lo.EC100.startIdx ] = CCM( EC100.loA, EC100.loH, E, tau, n, maxL, Lskip );

% EC108 hi
[ hi.EC108.OxmapA, hi.EC108.AxmapO, ~, ~ ] = CCM( EC108.hiO, EC108.hiA, E, tau, n, maxL, Lskip );

[ hi.EC108.OxmapH, hi.EC108.HxmapO, ~, ~ ] = CCM( EC108.hiO, EC108.hiH, E, tau, n, maxL, Lskip );

[ hi.EC108.OxmapC, hi.EC108.CxmapO, ~, ~ ] = CCM( EC108.hiO, EC108.hiC, E, tau, n, maxL, Lskip );

[ hi.EC108.CxmapA, hi.EC108.AxmapC, ~, ~ ] = CCM( EC108.hiC, EC108.hiA, E, tau, n, maxL, Lskip );

[ hi.EC108.CxmapH, hi.EC108.HxmapC, ~, ~ ] = CCM( EC108.hiC, EC108.hiH, E, tau, n, maxL, Lskip );

[ hi.EC108.AxmapH, hi.EC108.HxmapA, hi.EC108.Lvals, hi.EC108.startIdx ] = CCM( EC108.hiA, EC108.hiH, E, tau, n, maxL, Lskip );

% EC108 lo
[ lo.EC108.OxmapA, lo.EC108.AxmapO, ~, ~ ] = CCM( EC108.loO, EC108.loA, E, tau, n, maxL, Lskip );

[ lo.EC108.OxmapH, lo.EC108.HxmapO, ~, ~ ] = CCM( EC108.loO, EC108.loH, E, tau, n, maxL, Lskip );

[ lo.EC108.OxmapC, lo.EC108.CxmapO, ~, ~ ] = CCM( EC108.loO, EC108.loC, E, tau, n, maxL, Lskip );

[ lo.EC108.CxmapA, lo.EC108.AxmapC, ~, ~ ] = CCM( EC108.loC, EC108.loA, E, tau, n, maxL, Lskip );

[ lo.EC108.CxmapH, lo.EC108.HxmapC, ~, ~ ] = CCM( EC108.loC, EC108.loH, E, tau, n, maxL, Lskip );

[ lo.EC108.AxmapH, lo.EC108.HxmapA, lo.EC108.Lvals, lo.EC108.startIdx ] = CCM( EC108.loA, EC108.loH, E, tau, n, maxL, Lskip );

% save

save('data/EC100_EC108_HiVsLoAHBCoh_E10_tau10_maxL400_strictThresh.mat')

%% make plots (copy of code from script_20160519.m)

%% EC100 hi
data = hi.EC100;
dataName='EC100 high Am-Hp beta coherence';

% f1 = plotCCM(data.Lvals,data.OxmapA,data.AxmapO,E,tau,dataName,'OFC','Am','EC100_hiCoh_E10_tau10_OFCAm_strictThresh',[0.75 1]);
% 
% f2 = plotCCM(data.Lvals,data.OxmapH,data.HxmapO,E,tau,dataName,'OFC','Hp','EC100_hiCoh_E10_tau10_OFCHp_strictThresh',[0.75 1]);
% 
% f3 = plotCCM(data.Lvals,data.OxmapC,data.CxmapO,E,tau,dataName,'OFC','ACC','EC100_hiCoh_E10_tau10_OFCACC_strictThresh',[0.75 1]);
% 
% f4 = plotCCM(data.Lvals,data.CxmapA,data.AxmapC,E,tau,dataName,'ACC','Am','EC100_hiCoh_E10_tau10_ACCAm_strictThresh',[0.75 1]);
% 
% f5 = plotCCM(data.Lvals,data.CxmapH,data.HxmapC,E,tau,dataName,'ACC','Hp','EC100_hiCoh_E10_tau10_ACCHp_strictThresh',[0.75 1]);

f6 = plotCCM(data.Lvals,data.AxmapH,data.HxmapA,E,tau,dataName,'Am','Hp','EC100_hiCoh_E10_tau10_AmHp_strictThresh',[0.75 1]);

%% EC100 lo
data = lo.EC100;
dataName='EC100 low Am-Hp beta coherence';

% f1 = plotCCM(data.Lvals,data.OxmapA,data.AxmapO,E,tau,dataName,'OFC','Am','EC100_loCoh_E10_tau10_OFCAm_strictThresh',[0.75 1]);
% 
% f2 = plotCCM(data.Lvals,data.OxmapH,data.HxmapO,E,tau,dataName,'OFC','Hp','EC100_loCoh_E10_tau10_OFCHp_strictThresh',[0.75 1]);
% 
% f3 = plotCCM(data.Lvals,data.OxmapC,data.CxmapO,E,tau,dataName,'OFC','ACC','EC100_loCoh_E10_tau10_OFCACC_strictThresh',[0.75 1]);
% 
% f4 = plotCCM(data.Lvals,data.CxmapA,data.AxmapC,E,tau,dataName,'ACC','Am','EC100_loCoh_E10_tau10_ACCAm_strictThresh',[0.75 1]);
% 
% f5 = plotCCM(data.Lvals,data.CxmapH,data.HxmapC,E,tau,dataName,'ACC','Hp','EC100_loCoh_E10_tau10_ACCHp_strictThresh',[0.75 1]);

f6 = plotCCM(data.Lvals,data.AxmapH,data.HxmapA,E,tau,dataName,'Am','Hp','EC100_loCoh_E10_tau10_AmHp_strictThresh',[0.75 1]);

%% EC108 hi
data = hi.EC108;
dataName='EC108 high Am-Hp beta coherence';

f1 = plotCCM(data.Lvals,data.OxmapA,data.AxmapO,E,tau,dataName,'OFC','Am','EC108_hiCoh_E10_tau10_OFCAm_strictThresh',[0.75 1]);

f2 = plotCCM(data.Lvals,data.OxmapH,data.HxmapO,E,tau,dataName,'OFC','Hp','EC108_hiCoh_E10_tau10_OFCHp_strictThresh',[0.75 1]);

f3 = plotCCM(data.Lvals,data.OxmapC,data.CxmapO,E,tau,dataName,'OFC','ACC','EC108_hiCoh_E10_tau10_OFCACC_strictThresh',[0.75 1]);

f4 = plotCCM(data.Lvals,data.CxmapA,data.AxmapC,E,tau,dataName,'ACC','Am','EC108_hiCoh_E10_tau10_ACCAm_strictThresh',[0.75 1]);

f5 = plotCCM(data.Lvals,data.CxmapH,data.HxmapC,E,tau,dataName,'ACC','Hp','EC108_hiCoh_E10_tau10_ACCHp_strictThresh',[0.75 1]);

f6 = plotCCM(data.Lvals,data.AxmapH,data.HxmapA,E,tau,dataName,'Am','Hp','EC108_hiCoh_E10_tau10_AmHp_strictThresh',[0.75 1]);



%% EC108 l0
data = lo.EC108;
dataName='EC108 low Am-Hp beta coherence';

f1 = plotCCM(data.Lvals,data.OxmapA,data.AxmapO,E,tau,dataName,'OFC','Am','EC108_loCoh_E10_tau10_OFCAm_strictThresh',[0.75 1]);

f2 = plotCCM(data.Lvals,data.OxmapH,data.HxmapO,E,tau,dataName,'OFC','Hp','EC108_loCoh_E10_tau10_OFCHp_strictThresh',[0.75 1]);

f3 = plotCCM(data.Lvals,data.OxmapC,data.CxmapO,E,tau,dataName,'OFC','ACC','EC108_loCoh_E10_tau10_OFCACC_strictThresh',[0.75 1]);

f4 = plotCCM(data.Lvals,data.CxmapA,data.AxmapC,E,tau,dataName,'ACC','Am','EC108_loCoh_E10_tau10_ACCAm_strictThresh',[0.75 1]);

f5 = plotCCM(data.Lvals,data.CxmapH,data.HxmapC,E,tau,dataName,'ACC','Hp','EC108_loCoh_E10_tau10_ACCHp_strictThresh',[0.75 1]);

f6 = plotCCM(data.Lvals,data.AxmapH,data.HxmapA,E,tau,dataName,'Am','Hp','EC108_loCoh_E10_tau10_AmHp_strictThresh',[0.75 1]);









