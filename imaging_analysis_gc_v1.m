%% load the analog signals from the MDF files and extract the event
analog = SettingTwoPhoton('RVKC314_020218_imaging');
for i = 1:length(analog)
    info = analog2p(analog(i).data,10000);
    analog(i).frame = info.frame;
    analog(i).tone  = info.tone;
    taste =  {'S', 'N','C','Q','W'};
    idx   = [isempty(info.S), isempty(info.N),isempty(info.C), isempty(info.Q), isempty(info.W)];
    analog(i).taste = taste(~idx);
    analog(i).taste_ts = [info.S info.N info.C info.Q info.W];
    analog(i).lick     = info.lick;
    analog(i).time     = info.time;
    analog(i).tone_fr  = min(find(analog(i).frame > analog(i).tone));
    analog(i).taste_fr = min(find(analog(i).frame > analog(i).taste_ts(1)));

end
clear info
%% check the licking signal
figure
for i = 1:3
    subplot(1,3,i)
    plot(analog(i).time, analog(i).data(:,end));
    hold on 
    scatter(analog(i).lick, 0.2*ones(size(analog(i).lick)));
end
%% load the imaging data and visualize the trace
load('H:\GC Project\RVKC314\020518\Analysis\G\RVKC314\20180205\1\F_RVKC314_20180205_plane1_proc.mat')
ind = [dat.stat.iscell];
info.totalNeuron = sum(ind);
F = dat.Fcell{1};
info.F = F(find(ind==1),:);
info.idex = find(ind==1)
coeff = [dat.stat(find(ind==1)).neuropilCoefficient];
neuroF = dat.FcellNeu{1};
info.neuroF = neuroF(find(ind==1),:);
for i = 1:length(coeff)
    Fcell(i,:) = info.F(i,:)-coeff(i)*info.neuroF(i,:);

end
% Fcell = info.F-coeff.*info.neuroF;
Fcell(28,:) =[];
figure;
F0 = median(Fcell,2);
dF = Fcell - repmat(F0,1,size(Fcell,2));
dF_F = dF./repmat(F0,1,size(Fcell,2));
imagesc(dF_F)

figure;
for i =6:size(Fcell,1)/2
    plot(dF_F(2*(i-6)+6,:)+i,'k')
    hold on
    
end
%%
% take 2s before tone and 4 s after tone  6.2 Hz
% take 4s before taste and 4 s after taste
% it correspond as 12 frames before and 24 frames after
% after 5 times average, each trials has 60 frames (I acqure 300 frames for each trial)
trial_fr = 300/5;
for i =1:length(analog)
    analog(i).F = info.F(:,(i-1)*trial_fr+1 :i*trial_fr);
    frame_start = floor(analog(i).tone_fr/5)-11;
    frame_end   = floor(analog(i).tone_fr/5)+24;
    if frame_end >trial_fr
        frame_end = trial_fr;
    end
    analog(i).F_tone = analog(i).F(:,frame_start:frame_end);
    % extract the taste event 4s before taste
    frame_start = floor(analog(i).taste_fr(1)/5)-23;
    frame_end   = floor(analog(i).taste_fr(1)/5)+24;
    if frame_end >trial_fr
        frame_end = trial_fr;
    end
    analog(i).F_taste = analog(i).F(:,frame_start:frame_end);
    timestep =  0.0323*5 ; % 5 time average;  30.96 Hz scanning
    analog(i).timepointTone  = -11 * timestep: timestep: 24*timestep;
    analog(i).timepointTaste = -23 * timestep: timestep: 24*timestep;

    % take 1s before the tone as the baseline
    F0 = mean(analog(i).F_tone(:, find(analog(i).timepointTone<0 & analog(i).timepointTone>-1)),2);
    analog(i).dF_tone = (analog(i).F_tone - repmat(F0,1,size(analog(i).F_tone,2)))./repmat(F0,1,size(analog(i).F_tone,2));
    analog(i).dF_taste = (analog(i).F_taste - repmat(F0,1,size(analog(i).F_taste,2)))./repmat(F0,1,size(analog(i).F_taste,2));
end

%%
save('summary.mat','analog')
%% trying to reorganize the data to extract response base on tastant
% let's extract the sucrose responsive neurons
j = 1;
for i = 1:length(analog)
    if strcmp(analog(i).taste,'S')
        % As there should be 36 frames, however some trials have few frame,
        % so I pad zero to make it as 36 frames
        sucrose = zeros(size(analog(i).dF_taste,1),48);
        sucrose(:,1: size(analog(i).dF_taste,2))= analog(i).dF_taste;
        trial.sucrose(:,:,j) = sucrose;
        trial.tone_S(:,:,j)    = analog(i).dF_tone;
        j = j+1;
    end
end
figure; imagesc(analog(1).timepointTaste,[],squeeze(mean(trial.sucrose,3)));caxis([-0.5, 1])
title('S')
% NaCl
j = 1;
for i = 1:length(analog)
    if strcmp(analog(i).taste,'N')
        % As there should be 36 frames, however some trials have few frame,
        % so I pad zero to make it as 36 frames
        nacl = zeros(size(analog(i).dF_taste,1),48);
        nacl(:,1: size(analog(i).dF_taste,2))= analog(i).dF_taste;
        trial.nacl(:,:,j) = nacl;
        trial.tone_N(:,:,j)    = analog(i).dF_tone;

        j = j+1;
    end
end
figure; imagesc(analog(1).timepointTaste,[],squeeze(mean(trial.nacl,3)));caxis([-0.5, 1])
title('N')
% CA
j = 1;
for i = 1:length(analog)
    if strcmp(analog(i).taste,'C')
        % As there should be 36 frames, however some trials have few frame,
        % so I pad zero to make it as 36 frames
        ca = zeros(size(analog(i).dF_taste,1),48);
        ca(:,1: size(analog(i).dF_taste,2))= analog(i).dF_taste;
        trial.ca(:,:,j) = ca;
        trial.tone_C(:,:,j)    = analog(i).dF_tone;

        j = j+1;
    end
end
figure; imagesc(analog(1).timepointTaste,[],squeeze(mean(trial.ca,3)));caxis([-0.5, 1])
title('C')
% Q
j = 1;
for i = 1:length(analog)
    if strcmp(analog(i).taste,'Q')
        % As there should be 36 frames, however some trials have few frame,
        % so I pad zero to make it as 36 frames
        q = zeros(size(analog(i).dF_taste,1),48);
        q(:,1: size(analog(i).dF_taste,2))= analog(i).dF_taste;
        trial.q(:,:,j) = q;
        trial.tone_Q(:,:,j)    = analog(i).dF_tone;

        j = j+1;
    end
end
figure; imagesc(analog(1).timepointTaste,[],squeeze(mean(trial.q,3)));caxis([-0.5, 1])
title('Q')
% W
j = 1;
for i = 1:length(analog)
    if strcmp(analog(i).taste,'W')
        % As there should be 36 frames, however some trials have few frame,
        % so I pad zero to make it as 36 frames
        w = zeros(size(analog(i).dF_taste,1),48);
        w(:,1: size(analog(i).dF_taste,2))= analog(i).dF_taste;
        trial.w(:,:,j) = w;
        trial.tone_W(:,:,j)    = analog(i).dF_tone;

        j = j+1;
    end
end
figure; imagesc(analog(1).timepointTaste,[],squeeze(mean(trial.w,3)));caxis([-0.5, 1])
title('W')
%%
save('trial_info.mat','trial')

%%
% load trial_new.mat which was calculted by the server to test the
% responsibility (bootstrap): code: tasteResp.m which need trial_info.mat
% as input
%%
for i = 1: size(trial.S_prob,1)
    ind{1,i} = find(trial.S_prob(i,:)<0.01 & trial.S_prob(i,:)>0 )';
end
for i = 1: size(trial.N_prob,1)
    ind{2,i} = find(trial.N_prob(i,:)<0.01 & trial.N_prob(i,:)>0 )';
end
for i = 1: size(trial.C_prob,1)
    ind{3,i} = find(trial.C_prob(i,:)<0.01 & trial.C_prob(i,:)>0 )';
end
for i = 1: size(trial.Q_prob,1)
    ind{4,i} = find(trial.Q_prob(i,:)<0.01 & trial.Q_prob(i,:)>0 )';
end
for i = 1: size(trial.W_prob,1)
    ind{5,i} = find(trial.W_prob(i,:)<0.01 & trial.W_prob(i,:)>0 )';
end
for i = 1: size(trial.Tone_prob,1)
    ind{6,i} = find(trial.Tone_prob(i,:)<0.01 & trial.Tone_prob(i,:)>0 )';
end
%% Taste response
for i = 1:4
    for j = 1:size(ind,2)
        if isempty(ind{i,j})
            Tresp_boot(i,j) =0;
        elseif length(find(ind{i,j}<=39))>1 % 1 s is 6 frame; 2 s is 12 frame, 3 s is 18 frames; start is 24 frames; 2.5 s is 39 frames
            Tresp_boot(i,j)=1;
        else
            Tresp_boot(i,j)=0;
        end
    end
end
a = sum(Tresp_boot,1);
%% Cue response
for i = 1:size(ind,2)
    if isempty(ind{end,i})
        CueRes_boot(i) =0;
    elseif length(find(ind{end,i}<=24))>1 % 1 s is 6 frame; 2 s is 12 frame, start is 12 frames; 2 s is 24 frames
        CueRes_boot(i) =1;
    else
        CueRes_boot(i)=0;
    end    
end

%%
%% statistics based on m+sd; only excitatory
sig = 4;
for i = 1 : size(trial.sucrose,1)
    resp = squeeze(trial.sucrose(i,:,:));
    tone = squeeze(trial.tone_S(i,:,:));
    m = mean(tone(7:12,:),1); % 1 s before the tone
    s = std(tone(7:12,:),1);
    for j = 1:size(trial.sucrose,3)
        if ~isempty(find(resp(25:39,j)>m(j)+sig*s(j))) % 2 s within the taste
        test(i,j) = 1;
        elseif ~isempty(find(resp(25:39,j)< m(j)-sig*s(j))) % 2 s within the taste
            test(i,j) = 0;
        else
            test(i,j) = 0;
        end
    end
end
taste_response(1,:) = sum(test,2)';

% for i = 1 : size(trial.sucrose,1)
%     for j = 1:size(trial.sucrose,3)
%         resp = smooth(squeeze(trial.sucrose(i,:,j)),3);
%         tone = smooth(squeeze(trial.tone_S(i,:,j)),3);   
%         m = mean(tone(7:12)); % 1 s before the tone
%         s = std(tone(7:12));
%         
%         if ~isempty(find(resp(25:39)>m+sig*s)) % 2 s within the taste
%             test(i,j) = 1;
%         elseif ~isempty(find(resp(25:39)< m-sig*s)) % 2 s within the taste
%             test(i,j) = 0;
%         else
%             test(i,j) = 0;
%         end
%     end
% end
% taste_response(2,:) = sum(test,2)';


%
for i = 1 : size(trial.nacl,1)
    resp = squeeze(trial.nacl(i,:,:));
    tone = squeeze(trial.tone_N(i,:,:));
    m = mean(tone(7:12,:),1); % 1 s before the tone
    s = std(tone(7:12,:),1);
    for j = 1:size(trial.nacl,3)
        if ~isempty(find(resp(25:39,j)>m(j)+sig*s(j))) % 2 s within the taste
        test(i,j) = 1;
        elseif ~isempty(find(resp(25:39,j)< m(j)-sig*s(j))) % 2 s within the taste
            test(i,j) = 0;
        else
            test(i,j) = 0;
        end
    end
end
taste_response(2,:) = sum(test,2)';
%
for i = 1 : size(trial.ca,1)
    resp = squeeze(trial.ca(i,:,:));
    tone = squeeze(trial.tone_C(i,:,:));
    m = mean(tone(7:12,:),1); % 1 s before the tone
    s = std(tone(7:12,:),1);
    for j = 1:size(trial.ca,3)
        if ~isempty(find(resp(25:39,j)>m(j)+sig*s(j))) % 2 s within the taste
        test(i,j) = 1;
        elseif ~isempty(find(resp(25:39,j)< m(j)-sig*s(j))) % 2 s within the taste
            test(i,j) = 0;
        else
            test(i,j) = 0;
        end
    end
end
taste_response(3,:) = sum(test,2)';
%
for i = 1 : size(trial.q,1)
    resp = squeeze(trial.q(i,:,:));
    tone = squeeze(trial.tone_Q(i,:,:));
    m = mean(tone(7:12,:),1); % 1 s before the tone
    s = std(tone(7:12,:),1);
    for j = 1:size(trial.ca,3)
        if ~isempty(find(resp(25:39,j)>m(j)+sig*s(j))) % 2 s within the taste
        test(i,j) = 1;
        elseif ~isempty(find(resp(25:39,j)< m(j)-sig*s(j))) % 2 s within the taste
            test(i,j) = 0;
        else
            test(i,j) = 0;
        end
    end
end
taste_response(4,:) = sum(test,2)';
%
for i = 1 : size(trial.w,1)
    resp = squeeze(trial.w(i,:,:));
    tone = squeeze(trial.tone_W(i,:,:));
    m = mean(tone(7:12,:),1); % 1 s before the tone
    s = std(tone(7:12,:),1);
    for j = 1:size(trial.w,3)
        if ~isempty(find(resp(25:39,j)>m(j)+sig*s(j))) % 2 s within the taste
        test(i,j) = 1;
        elseif ~isempty(find(resp(25:39,j)< m(j)-sig*s(j))) % 2 s within the taste
            test(i,j) = 0;
        else
            test(i,j) = 0;
        end
    end
end
taste_response(5,:) = sum(test,2)';
% tone response
clear test
tone = cat(3,trial.tone_S, trial.tone_N, trial.tone_C, trial.tone_Q, trial.tone_W);
for i = 1 : size(tone,1)
    resp = squeeze(tone(i,:,:));
    baseline = squeeze(tone(i,:,:));
    m = mean(baseline(7:12,:),1); % 1 s before the tone
    s = std(tone(7:12,:),1);
    for j = 1:size(baseline,2) % loop through each trial
        if ~isempty(find(resp(13:24,j)>m(j)+sig*s(j))) % 2 s within the tone
        test(i,j) = 1;
        elseif ~isempty(find(resp(13:24,j)< m(j)-sig*s(j))) % 2 s within the taste
            test(i,j) = 0;
        else
            test(i,j) = 0;
        end
    end
end
%%
% figure; imagesc(analog(1).timepointTone,[],squeeze(mean(tone,3)));caxis([-0.5, 1])
% title('Tone')
% taste_response(6,:) = sum(test,2)';
%
% I take that at 2 trials show significant higher response
% taste_response(find(taste_response<4))=0; % each row is a tastant
%
% for i = 1:size(taste_response,2)
%      re = find(taste_response(1:5,i)>0);
%      if isempty(re)
%          a(i) = 0;
%      elseif length(re) == 1;
%          a(i) = 1;
%      elseif length(re) == 2;
%          a(i) =2;
%      elseif length(re) == 3;
%          a(i) = 3;
%      elseif length(re) == 4;
%          a(i) = 4;
%      elseif length(re) ==5;
%          a(i) =5;
%      end
% end
%%
responeRatio = 1-length(find(a==0))./length(a);
responeRatio1 = length(find(a==1))./length(a);
responeRatio2 = length(find(a==2))./length(a);
responeRatio3 = length(find(a==3))./length(a);
responeRatio4 = length(find(a==4))./length(a);
% responeRatio5 = length(find(a==5))./length(a);
figure;
subplot(1,3,1)
plot([1,2,3,4],[responeRatio1,responeRatio2,responeRatio3,responeRatio4],'-o')
xlim([0,5])
set(gca,'XTickLabels',({'','1','2','3','4',''}))
box off
xlabel('Number of tastant');
ylabel('% Taste repsonive')
set(gca,'TickDir','out');
% plot the response ratio for each tastant
% for i =1:5
%     ratio(i) = length(find(taste_response(i,:)>0))/size(taste_response,2);
% end
for i =1:4
    ratio(i) = length(find(Tresp_boot(i,:)>0))/size(Tresp_boot,2);
end
subplot(1,3,2)
bar(ratio,'FaceColor',[1,1,1])
ylim([0,0.4])
xlim([0,6])
% yticks([0,0.25,0.5])
% xticks([0,1,2,3,4,5,6])
set(gca,'XTickLabels',({'S','N','C','Q'}))
set(gca,'YTick',([0,0.2]))
box off
ylabel('% Taste responsive')

subplot(1,3,3)
bar(sum(CueRes_boot)/length(CueRes_boot),'FaceColor',[0.5,0.5,0.5])
ylim([0,0.2])
xlim([0,2])
set(gca,'YTick',([0,0.2]))
ylabel('% Cue response')

%% summarize taste responsive neurons
%% To list
% visulize each active neuron very easily
% Visulize each trial very easilty
% overlapping; spatial analysis
% population coding; noise correlation
%
a = find(Tresp_boot(1,:) ==1);
b = find(Tresp_boot(2,:) ==1);
c = find(Tresp_boot(3,:) ==1);
d = find(Tresp_boot(4,:) ==1);
all = unique([a, b, c,d]);
figure; imagesc(squeeze(mean(trial.nacl(b,:,:),3)));
%%
dF_FS = mean(trial.sucrose,3);
dF_FN = mean(trial.nacl,3);
dF_FC = mean(trial.ca,3);
dF_FQ = mean(trial.q,3);
dF_FW = mean(trial.w,3);
%%
ex = 9;
figure; 
plot(analog(1).timepointTaste,dF_FS(ex,:))
hold on
plot(analog(1).timepointTaste,dF_FN(ex,:)+0.5)
plot(analog(1).timepointTaste,dF_FC(ex,:)+1)
plot(analog(1).timepointTaste,dF_FQ(ex,:)+1.5)
plot(analog(1).timepointTaste,dF_FW(ex,:)+2)
plot([0,0],[0,2.5],'r')
legend({'S','N','C','Q','W'})
title('response to all')
%%

ind = find(test(1,:)>0);
figure; imagesc(analog(1).timepointTaste,[],squeeze(mean(trial.sucrose(ind,:,:),3)));caxis([-0.5, 1])
title('S')
%%
clear test
sig = 5;
s = squeeze(mean(trial.sucrose,3));
s_tone = squeeze(mean(trial.tone_S,3));
for i = 1:size(s,1)
    m = mean(s_tone(i,7:12));
    st = std(s_tone(i,7:12));
    if ~isempty(find(s(i,:)>m+sig*st))
        test(i) = 1;
    else
        test(i) = 0;
    end
end

%%
figure; 
for i =1:size(tone,3)
    plot(analog(1).timepoint,squeeze(tone(10,:,i)))
    hold on
end

%% plot representative trace
% sucrose response trace

for i =1:88
    figure
cell = squeeze(trial.sucrose(i,:,:));
cell_avg = mean(cell,2);
plot(analog(1).timepointTaste,cell_avg)
hold on
cell = squeeze(trial.nacl(i,:,:));
cell_avg = mean(cell,2);
plot(analog(1).timepointTaste,cell_avg+1)
hold on
cell = squeeze(trial.ca(i,:,:));
cell_avg = mean(cell,2);
plot(analog(1).timepointTaste,cell_avg+2)
hold on
cell = squeeze(trial.q(i,:,:));
cell_avg = mean(cell,2);
plot(analog(1).timepointTaste,cell_avg+3)
end
%% plot the spatial map
% for sucrose
s_ind = find(Tresp_boot(1,:)==1);
ind = [dat.stat.iscell];
loc = dat.stat(find(ind ==1));
figure;imshow(image)
% hold on
% for i = 1:length(s_ind)
%     scatter(loc(s_ind(i)).xpix,512-loc(s_ind(i)).ypix,'.')
%     hold on
%     
% end
% image2 = zeros(512,512);
i = s_ind;
%  for j =1:length(i)
%       for k = 1:length(loc(i(j)).xpix)
%        image2(loc(i(j)).ypix(k),loc(i(j)).xpix(k))=1;
%       end
%  end
 hold on
 for j =1:length(i)    
       scatter(loc(i(j)).xpix,loc(i(j)).ypix,'.m')
       hold on
 end
xlim([1,512]);ylim([1,512])
% imagesc(image2)
[B,L] = bwboundaries(image,'noholes');
figure;imshow(image)
hold on
for k = 1:length(B)
   boundary = B{k};
   plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end