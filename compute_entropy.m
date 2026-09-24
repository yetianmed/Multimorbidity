clear all
close all
addpath ./cbrewer
gender=1;
if gender==1
    load('burden_paths_male.mat','T','cluster_name','dgroup');
else
    load('burden_paths_female.mat','T','cluster_name','dgroup');
end

% T is in a table format organised as below
% each subject has one or more rows depending on the number of states
% {'ID'}    {'Age'}    {'State'}    {'Cluster'}    {'Cardiovascular'}    {'Eye'}
%  {'Gastrointestinal'}    {'Kidney'}    {'Liver'}    {'Mental'}    {'Metabolic'}
%    {'Musculoskeletal'}    {'Neurological'}    {'Reproductive'}    {'Respiratory'}
%  {'Death'}

X=table2array(T);

system=unique(dgroup(:,2)); % dgroup comprises two columns. 1st individual diagnosis; 2nd disease category
for i=1:length(system)
    d(i)=length(find(contains(dgroup(:,2),system(i))));
end

Clusters=length(cluster_name);%4; %number of clusters
States=4;   %number of states

Boot=100; %number of boostrap samples

h_norm=zeros(States,Clusters,Boot); %entropy normalized to cluster-independent entropy
b_norm=zeros(States,Clusters,Boot); %burden normalized to cluster-independent burden
hb=zeros(States,Clusters,Boot); %ratio of burden to entropy (non-normalized)
hb_norm=zeros(States,Clusters,Boot); %b_norm/h_norm
hb_ref=zeros(States,1,Boot); %burden to entropy ratio independent of cluster (reference)
frst=0; 
for boot=1:Boot
    x=datasample(X,size(X,1),1,'Replace',true); %bootrap sample

    data=x(:,5:end-1); %isolate LTC counts 
    

    b=zeros(States,Clusters); %total burden
    h=zeros(States,Clusters); %entropy
    for n=1:Clusters
        for i=1:States
            ind=find(x(:,4)==n & x(:,3)==i);
            h(i,n)=jointEntropyDiscrete(data(ind,:));
            b(i,n)=sum(mean(data(ind,:)));
        end
    end

    %Global (independent of cluster type)
    h_all=zeros(States,1);
    b_all=zeros(States,1);
    for i=1:States
        ind_all=find(x(:,3)==i);
        h_all(i)=jointEntropyDiscrete(data(ind_all,:));
        b_all(i)=sum(mean(data(ind_all,:)));
    end

    h_norm(:,:,boot)=h./repmat(h_all,1,Clusters);
    b_norm(:,:,boot)=b./repmat(b_all,1,Clusters);

    hb(:,:,boot)=b./h;
    hb_norm(:,:,boot)=b_norm(:,:,boot)./h_norm(:,:,boot);
    hb_ref(:,1,boot)=b_all./h_all;

    show_progress(boot,Boot,frst); frst=1; 
end

%Draw figures
addpath ./cbrewer
cmap=cbrewer('qual','Set2',5); 
if gender==0
   cmap=cmap([2,5],:);
end
marker_type={'^','*','v','x','o','+'};
if gender==0
    marker_type=marker_type([2,5]);
end
alpha_rng=[1:States]; 
facealpha=0.9; 

%%%
%NORMALIZED burden-entropy ratio
%%%
hf=figure; hf.Position=[100,100,300,280]; hf.Color='w'; 
for i=1:Clusters
    mm=mean(hb_norm(:,i,:),3);
    sem=std(hb_norm(:,i,:),0,3);  
    upper_bound = mm + 1.96*sem;
    lower_bound = mm - 1.96*sem;
    hShade = fill([alpha_rng, fliplr(alpha_rng)], [upper_bound', fliplr(lower_bound')], cmap(i,:), ...
    'EdgeColor', 'none', ...      % Removes the border line around the shading
    'FaceAlpha', facealpha);        % Sets transparency (25% opacity)
    hold on;
    hs(i)=scatter(alpha_rng,mm,150, 'Marker', marker_type{i}, 'MarkerEdgeColor', cmap(i,:),'LineWidth', 1.5);
end
ha=gca; ha.Box='off'; ha.LineWidth=3;
ha.XTick=1:4; ha.XTickLabel={'I','II','III','IV'}; 
ha.XLabel.String='State'; ha.YLabel.String='Burden/Entropy'; ha.FontSize=14; 
axis tight; 

%Draw the legend
ShowLegend=1; 
if ShowLegend
    [hl,item_handles]=legend(hs,cluster_name,'Location','NorthEast'); 
    text_objects = findobj(item_handles, 'Type', 'text');
    legend('boxoff');

    legend_markers = findobj(item_handles, 'type', 'line', '-or', 'type', 'patch');

    % 3. Increase the radius size of the hollow markers
    set(legend_markers, 'MarkerSize', 12);

    % Make sure the legend uses normalized units for drawing bounds
    set(hl, 'Units', 'normalized');
    lgd_pos = hl.Position; % [left, bottom, width, height] of the legend box

    % 4. Loop backward through text objects to create the background blocks
    num_entries = length(text_objects);
    for idx = 1:num_entries
        txt = text_objects(idx);
        pos = txt.Position; % Internal text coordinate [X, Y, Z] relative to legend

        color_idx =idx;
        % FIX: Calculate positions relative to the entire Figure canvas
        box_width = lgd_pos(3) * 0.94;
        if gender==1
        box_height = lgd_pos(4) * 0.15;
        else
        box_height = lgd_pos(4) * 0.25;
        end
        box_left = lgd_pos(1) + (lgd_pos(3) * 0.03);
        box_bottom = lgd_pos(2) + (lgd_pos(4) * pos(2)) - (box_height / 2);

        % FIX: First argument changed to 'gcf' (get current figure)
        annotation(gcf, 'rectangle', [box_left, box_bottom, box_width, box_height], ...
            'FaceColor', cmap(color_idx, :), ...
            'EdgeColor', 'none','FaceAlpha',0.3);
    end

    % Crucial: Push the legend object back to the front layer over the annotations
    uistack(hl, 'top');
    %%%
end

