clear all
clc
close all

%% version 1 flowchart:
% 1. open maps
% 2. calculate Zscore
% 3. filter out extreme Zscore singletons
% 4. go back to the maps and filter out rows and columns >1.5 IQR
% 5. re-calculate Zscores and compare them


%% options:
distance = 5; % distance for correlation between neighbors and for variance analysis
binsize = 6000; % genomic size of the bin
startcoord = 99005149; % start coordinate of the first bin in the square matrix (after cutting with convert.sh, see below)
ZEROS = 'true'; % if 'true', Zscores are calculated keeping 0s in the data; otherwise if 'false' they are discarded
tablename = '20160822_5C-Samples.xlsx'; % Xls file with information on samples, format like the example provided
color=[0.8 0.8 0.8]; %color for NaNs [0.8 0.8 0.8]=grey
viewpoint=99021000; %viewpoint for virtual 4C profile
Inversion = 'false'; % Attention! IT IS VALID FOR ALL MAPS! true if you want the genome browser to display the correct (inverted) genomic region when you click on the maps
Zscore_range = [-2,2]; % range of plot for Zscore maps
log_ratio_range = [-3,3]; % range of plot for log ratio maps
Zscore_range_filtered = [-3,3]; % range of plot for Zscore maps after filters
top = 2; % percentage of top changing interaction from log2fc
saturation = 200; %max value in the heatmap plots

% samples:
wt_sample = 'E14-WT_pooled';
mut_sample = 'LinxCBS-inv_pooled';

% file names of WT and mutant sample (these are Nicolas' pirwise matrices, wither raw or ICEd):
wt_filename = '../FINAL-5C-MAPS_April2017/04_Pooled_maps/01_WT/E14-WT_pooled_binned.mat';
mut_filename = '../FINAL-5C-MAPS_April2017/04_Pooled_maps/03_Inversions/LinxCBS-inv_pooled_inversion_binned.mat';


%% 1. data input
disp(['** Hello Rafael, I am processing samples ', wt_sample, ' and ', mut_sample, ' **'])

% look up if the mutant is a deletion or inversion and find deletion coordinates
table=readtable(tablename); % imports the database table with all 5C samples
disp(['Importing sample information from table ', tablename])
%table(:,'Var11') = []; % get rid of 11th column (useless)

%needed for luca?s version of matlab:
%table(:,'Var10') = []; % get rid of 10th column (useless)
%table(:,'Var9') = []; % get rid of 9th column (useless)
%end of what is needed for luca?s version of matlab


names = table.Sample; % this is the sample names
mut_info = table{strcmp(mut_sample, names), :}; % finds sample in the table and extracts the corresp. row with info
if strcmp(mut_info(3),'Deletion')
     disp(['The sample ', mut_sample, ' was detected as a Deletion'])
     % determine position and size of the deletion (useful for calculating distance matrix below and obtain correct Zscore!)
     deletion_start = str2double(mut_info(7));
     deletion_end = str2double(mut_info(8));
     deletion_size = deletion_start-deletion_end;
     disp(['Deletion start Chrx:', num2str(deletion_start)])
     disp(['Deletion end Chrx:', num2str(deletion_end)])

     % determine start and end bin of the deletion, excluding both start and end bins:
     deletion_start_bin = fix((deletion_start -startcoord) / binsize +1);
     deletion_end_bin = ceil((deletion_end -startcoord) / binsize);
     deletion_size_bin = deletion_end_bin - deletion_start_bin+1;
     disp(['Corresponding to start bin ', num2str(deletion_start_bin)])
     disp(['and end bin ', num2str(deletion_end_bin)])

end

if strcmp('raw', wt_filename)
    disp('Working on raw data')
elseif strcmp('ice', wt_filename)
    disp('Working on raw data')
else
    disp('Attention! Not specified if raw or ice data')
end


%% Pairwise matrix conversion to square matrix format, cutting to uniform region
%     convert to matrix format and cut the matrices to ensure that they have the same size
%     using Zhan's 'conver.sh' Bash script
disp('Converting pairwise format into square matrix...')
unix(['./convert.sh ', wt_filename]);
unix(['./convert.sh ', mut_filename]);

% import WT 5C square matrix
wt = importdata([wt_filename, '.matrix'], ' ');
% import mutant 5C square matrix
mut = importdata([mut_filename, '.matrix'], ' ');

% extract matrix entries:
wt_data = wt.data; % contains the 'data' field of the structure wt
mut_data = mut.data; % idem
% simmetrize the matrix (only half of it is given as an input):
for i=1:size(wt_data,1)
    for j=1:size(wt_data,1)
        wt_data(j,i) = wt_data(i,j);
        mut_data(j,i) = mut_data(i,j);
    end
end
% genomic coordinates of beginning of bins (assuming 5kb binning):
gencoord = ([1:size(wt_data,1)]-1) * binsize + startcoord;

% normalize by the total number of reads and by a common read number:
wt_data = wt_data ./ sum(sum(wt_data)) *10e6;
mut_data = mut_data ./ sum(sum(mut_data)) *10e6;
disp('Matrices imported and normalized by the number of reads x 10million')

%% PRE-FILTERING:
% sum or rows:
row_sum_wt = sum(wt_data,2);
row_sum_mut = sum(mut_data,2);
% assign rows where sum of counts is striclty 0 to NaN (non-working primers)
wt_data(row_sum_wt==0,:)=NaN;
mut_data(row_sum_mut==0,:)=NaN;
wt_data(:,row_sum_wt==0)=NaN;
mut_data(:,row_sum_mut==0)=NaN;
disp('Assigning masked bins to Nan')

%% plot WT and mutant matrices 
disp(['% of Max for ',wt_sample,': ',num2str(1-sum(wt_data>saturation)/sum(not(isnan(wt_data))))])
figure('Name', 'WT')
    appo=wt_data;
    appo(isnan(appo))=-1;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [-3, saturation])
    
figure('Name', 'WT')
    appo=wt_data(209:385,209:385);
    appo(isnan(appo))=-1;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [-3, saturation])
disp(['% of Max for ',mut_sample,': ',num2str(1-sum(mut_data>saturation)/sum(not(isnan(mut_data))))])
figure('Name', 'mut')
    appo=mut_data;
    appo(isnan(appo))=-1;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [-3, saturation])
    axis square

figure('Name', 'mut')
    appo=mut_data(209:385,209:385);
    appo(isnan(appo))=-1;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [-3, saturation])
    axis square
     
% load WhiteRedBlack colormap:
map3nn=load('WhiteOrangeRedBrown');
% set colormap in figure 1 and 2
set(1:4,'Colormap',map3nn.map)
%saveas([3],mut_sample,'bmp')
%saveas([4],[mut_sample,'_zoom'],'bmp')

%% Correlation of neighborhoods
% build a matrix of distances for WT sample
disp('Doing correlation analysis of counts at different distances')
for i=1:size(wt_data,1)
    for j=1:size(wt_data,1)
        dist_wt(i,j) = i-j;
    end
end

%correlation matrix
correlation = [zeros(size(wt_data,1)-distance-10,1) zeros(size(wt_data,1)-distance-10,1)];
for i=1:(size(wt_data,1)-distance-10)
    x=wt_data(dist_wt==i)*i;
    x=x(1:length(x)-distance);
    y=wt_data(dist_wt==i+distance)*(i+distance);
    nan=unique([find(isnan(x)) ;find(isnan(y))]);
    x(nan)=[];
    y(nan)=[];
    correlation(i,1)=i;
    correlation(i,2)=corr(log2(x),log2(y),'Type','Spearman');
    
end
figure('Name',['Correlation between bead |i-j|=fixed and bead |i-j|+' num2str(distance)])
plot(correlation(:,1),correlation(:,2))
xlabel('|i-j|')
ylabel('Correlation value (spearman)')


%% neighborhood Coefficient of variation (CV)
disp('Doing coefficient of variation analysis to discard noisy 5C interactions')
normalised_wt=wt_data-wt_data;
n2=normalised_wt;
normalised_mut=n2;
n2_mut=n2;
for i=1:(size(wt_data,1))
    for j=1:(size(wt_data,1))
        conta_wt=0;
        conta_mut=0;
       
        if (abs(i-j)>=15)
            for h=-distance:distance
                for k=-distance:distance
                    if(i+h>0 & j+k>0 & i+h<=(size(wt_data,1)) & j+k<=(size(wt_data,1)))
                        if(~(isnan(wt_data(i+h,j+k))))
                            conta_wt=conta_wt+1;
                            normalised_wt(i,j)=normalised_wt(i,j)+wt_data(i+h,j+k);
                            n2(i,j)=n2(i,j)+wt_data(i+h,j+k)*wt_data(i+h,j+k);                        
                        end
                        if (~(isnan(mut_data(i+h,j+k))))
                            conta_mut=conta_mut+1;
                            normalised_mut(i,j)=normalised_mut(i,j)+mut_data(i+h,j+k);
                            n2_mut(i,j)=n2_mut(i,j)+mut_data(i+h,j+k)*mut_data(i+h,j+k);
                        end
                    end
                end
            end
        normalised_wt(i,j)=normalised_wt(i,j)/conta_wt;
        if(n2(i,j)/conta_wt-normalised_wt(i,j)*normalised_wt(i,j)>=0)
            n2(i,j)=sqrt(n2(i,j)/conta_wt-normalised_wt(i,j)*normalised_wt(i,j))/normalised_wt(i,j);
        else
            n2(i,j)=NaN;
        end
        normalised_mut(i,j)=normalised_mut(i,j)/conta_mut;
        if(n2_mut(i,j)/conta_mut-normalised_mut(i,j)*normalised_mut(i,j)>=0)
            n2_mut(i,j)=sqrt(n2_mut(i,j)/conta_mut-normalised_mut(i,j)*normalised_mut(i,j))/normalised_mut(i,j);
        else
            n2_mut(i,j)=NaN;
        end
        else
            n2(i,j)=NaN;
            n2_mut(i,j)=NaN;
        end
    end
end
figure('Name', 'variance/mean WT')
    imagesc(n2)
    set(gca, 'CLim', [0, 2])
    colorbar
figure('Name','Histogram of CV WT')    
    hist(log2(n2(:)),100)
    xlabel('log2(CV)')
    
disp('Click on the histogram plot of WT and set the cut-off for Hi-C interactions')
[x,y] = ginput(1);
cutoff_wt=2^(x);
disp(['CV cutoff of wt = ' num2str(cutoff_wt)])

figure('Name', 'variance/mean MUT')
    imagesc(n2_mut)
    set(gca, 'CLim', [0, 2])
    colorbar
figure('Name','Histogram of CV MUT')    
    hist(log2(n2_mut(:)),100)
    xlabel('log2(CV)')
    
disp('Click on the histogram plot of MUT and set the cut-off for Hi-C interactions')
[x,y] = ginput(1);
cutoff_mut=2^(x);
disp(['CV cutoff of mut = ' num2str(cutoff_wt)])

disp('Filtering maps with CV > cutoff ')
disp('% of reads discarded ')
disp(['     1-Wt: ' num2str(100*sum(sum(n2>cutoff_wt | isnan(n2)))/length(n2(:)))])
disp(['     2-Mut: ' num2str(100*sum(sum(n2_mut>cutoff_mut | isnan(n2_mut)))/length(n2_mut(:)))])

wt_data((n2>cutoff_wt | n2_mut>cutoff_mut))=NaN;
mut_data((n2>cutoff_wt | n2_mut>cutoff_mut))=NaN;


%% 2. Zscore calculation and filtering of singletons

% build a matrix of distances for WT sample
for i=1:size(wt_data,1)
    for j=1:size(wt_data,1)
        dist_wt(i,j) = i-j;
    end
end
% modify parts affected by deletion in the mutant, if the sample is a deletion:
dist_mut = dist_wt;
if strcmp(mut_info(3),'Deletion')
    disp('Correcting distance matrix for the presence of a deletion...')
    % correct distance matrix for the presence of a deletion
    for i = 1:deletion_start_bin-1
        for j = deletion_end_bin:size(wt_data,1)
            dist_mut(i,j)=dist_wt(i,j)+deletion_size_bin;
            dist_mut(j,i)=dist_wt(j,i)-deletion_size_bin;
        end
    end
end


% extract matrix elements at fixed distances
for i=1:size(wt_data,1)
    wt_diag{i} = wt_data(dist_wt==i);
    mut_diag{i} = mut_data(dist_mut==i);
end

switch ZEROS
    
    case 'true' % keeping zeros
        disp('NB: keeping zeros in Zscore calculation!')
        % calculate means over diagonals:
        for i=1:size(wt_data,1)
            mean_wt_diag(i) = nanmean(wt_diag{i});
            std_wt_diag(i) = nanstd(wt_diag{i});
            mean_mut_diag(i) = nanmean(mut_diag{i});
            std_mut_diag(i) = nanstd(mut_diag{i});
        end
        % calculate Zscores
        for i=1:size(wt_data,1)
            for j=1:size(wt_data,1)
                Zscore_wt(i,j) = (wt_data(i,j) - mean_wt_diag(abs(dist_wt(i,j))+1)) ./ std_wt_diag(abs(dist_wt(i,j))+1);
                Zscore_mut(i,j) = (mut_data(i,j) - mean_mut_diag(abs(dist_mut(i,j))+1)) ./ std_mut_diag(abs(dist_mut(i,j))+1) ;
            end
        end
        Zscore_wt(isinf(Zscore_wt))=NaN;
        Zscore_mut(isinf(Zscore_mut))=NaN;
        figure('Name','Zscore WT with zeros')
            appo=Zscore_wt;
            appo(appo<Zscore_range(1))=Zscore_range(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range(1)-0.1,Zscore_range(2)])
        figure('Name','Zscore Mut with zeros')
            appo=Zscore_mut;
            appo(appo<Zscore_range(1))=Zscore_range(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range(1)-0.1,Zscore_range(2)])
    
        
    case 'false' % excluding zeros
        disp('NB: discarding zeros in Zscore calculation')
        % calculate means over diagonals:
        for i=1:size(wt_data,1)
            mean_wt_diag(i) = nanmean(nonzeros(wt_diag{i}));
            std_wt_diag(i) = nanstd(nonzeros(wt_diag{i}));
            mean_mut_diag(i) = nanmean(nonzeros(mut_diag{i}));
            std_mut_diag(i) = nanstd(nonzeros(mut_diag{i}));
        end
        % calculate Zscores
        for i=1:size(wt_data,1)
            for j=1:size(wt_data,1)
                Zscore_wt(i,j) = (wt_data(i,j) - mean_wt_diag(abs(dist_wt(i,j))+1)) ./ std_wt_diag(abs(dist_wt(i,j))+1);
                Zscore_mut(i,j) = (mut_data(i,j) - mean_mut_diag(abs(dist_mut(i,j))+1)) ./ std_mut_diag(abs(dist_mut(i,j))+1) ;
            end
        end
        Zscore_wt(isinf(Zscore_wt))=NaN;
        Zscore_mut(isinf(Zscore_mut))=NaN;        
        figure('Name','Zscore WT without zeros')
            appo=Zscore_wt;
            appo(appo<Zscore_range(1))=Zscore_range(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range(1)-0.1,Zscore_range(2)])
        figure('Name','Zscore Mut without zeros')
            appo=Zscore_mut;
            appo(appo<Zscore_range(1))=Zscore_range(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range(1)-0.1,Zscore_range(2)])
end



% plot Zscore distributions:
figure('Name', 'distribution of Zscores')
subplot(2,1,1)    
    hist(nonzeros(Zscore_wt(:)), 20)
    xlabel 'Zscore'
    ylabel 'occurrence'
    axis([-5 10 0 8e4])
    legend 'WT'
subplot(2,1,2)    
    hist(nonzeros(Zscore_mut(:)), 20)
    xlabel 'Zscore'
    ylabel 'occurrence'
    axis([-5 10 0 8e4])
    legend 'mut'

%% filter out singletons:
disp('Filtering out Zscore sigletons outside .999 quantile')
threshold_wt = quantile(nonzeros(Zscore_wt(:)), 0.999);
threshold_mut = quantile(nonzeros(Zscore_mut(:)), 0.999);

wt_data_filtered = wt_data;
mut_data_filtered = mut_data;
wt_data_filtered(Zscore_wt>threshold_wt & Zscore_mut>threshold_mut) = NaN;
mut_data_filtered(Zscore_wt>threshold_wt & Zscore_mut>threshold_mut) = NaN;
    

%% 3. filter on row-column reads
disp('Filtering out rows and columns with sums > 3rd quantile + 1.5 IQR and below 1st quantile -1.5IQR')
switch ZEROS 
    case 'true' % keeping zeros
        disp('Keeping zeros in the row sum')
        % determine sum of rows on filtered matrix including zeros:
        row_sum_filtered_wt = nansum(wt_data_filtered,2);
        row_sum_filtered_mut = nansum(mut_data_filtered,2);
    case 'false' % removing zeros
        disp('Removing zeros in the row sum')
        % determine sum of rows on filtered matrix removing zeros:
        for i=1:size(wt_data,1)
            row_sum_filtered_wt(i,1) = nansum(nonzeros(wt_data_filtered(i,:)));
            row_sum_filtered_mut(i,1) = nansum(nonzeros(mut_data_filtered(i,:)));
        end
end

% definition of outliers
%quantile 3
q = quantile(row_sum_filtered_mut(row_sum_filtered_mut>0),[0.25 0.75]);
l_threshold=q(1)-(q(2)-q(1))*1.5;
u_threshold=q(2)+(q(2)-q(1))*1.5;
row_filter=(row_sum_filtered_mut>u_threshold | row_sum_filtered_mut<l_threshold);
mut_data_filtered(row_filter==1,:)=NaN;
mut_data_filtered(:,row_filter==1)=NaN;

q = quantile(row_sum_filtered_wt(row_sum_filtered_wt>0),[0.25 0.75]);
l_threshold=q(1)-(q(2)-q(1))*1.5;
u_threshold=q(2)+(q(2)-q(1))*1.5;
row_filter=(row_sum_filtered_wt>u_threshold | row_sum_filtered_wt<l_threshold);
wt_data_filtered(row_filter==1,:)=NaN;
wt_data_filtered(:,row_filter==1)=NaN;


%% 4. Zscore calculation on filtered maps
% extract matrix elements at fixed distances
for i=1:size(wt_data,1)
    wt_diag{i} = wt_data_filtered(dist_wt==i);
    mut_diag{i} = mut_data_filtered(dist_mut==i);
end

disp('Calculating Zscores on filtered contact maps')
switch ZEROS
    
    case 'true' % keeping zeros
        % calculate means over diagonals:
        for i=1:size(wt_data,1)
            mean_wt_diag(i) = nanmean(wt_diag{i});
            std_wt_diag(i) = nanstd(wt_diag{i});
            mean_mut_diag(i) = nanmean(mut_diag{i});
            std_mut_diag(i) = nanstd(mut_diag{i});
        end
        % calculate Zscores
        for i=1:size(wt_data,1)
            for j=1:size(wt_data,1)
                Zscore_wt(i,j) = (wt_data_filtered(i,j) - mean_wt_diag(abs(dist_wt(i,j))+1)) ./ std_wt_diag(abs(dist_wt(i,j))+1);
                Zscore_mut(i,j) = (mut_data_filtered(i,j) - mean_mut_diag(abs(dist_mut(i,j))+1)) ./ std_mut_diag(abs(dist_mut(i,j))+1) ;
            end
        end
        Zscore_wt(isinf(Zscore_wt))=NaN;
        Zscore_mut(isinf(Zscore_mut))=NaN;  
        figure('Name','Zscore WT with zeros after filters')
            appo=Zscore_wt;
            appo(appo<Zscore_range_filtered(1))=Zscore_range_filtered(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range_filtered(1)-0.1,Zscore_range_filtered(2)])
        figure('Name','Zscore Mut with zeros after filters')
            appo=Zscore_mut;
            appo(appo<Zscore_range_filtered(1))=Zscore_range_filtered(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range_filtered(1)-0.1,Zscore_range_filtered(2)])
    
        
    case 'false' % excluding zeros
        % calculate means over diagonals:
        for i=1:size(wt_data,1)
            mean_wt_diag(i) = nanmean(nonzeros(wt_diag{i}));
            std_wt_diag(i) = nanstd(nonzeros(wt_diag{i}));
            mean_mut_diag(i) = nanmean(nonzeros(mut_diag{i}));
            std_mut_diag(i) = nanstd(nonzeros(mut_diag{i}));
        end
        % calculate Zscores
        for i=1:size(wt_data,1)
            for j=1:size(wt_data,1)
                Zscore_wt(i,j) = (wt_data_filtered(i,j) - mean_wt_diag(abs(dist_wt(i,j))+1)) ./ std_wt_diag(abs(dist_wt(i,j))+1);
                Zscore_mut(i,j) = (mut_data_filtered(i,j) - mean_mut_diag(abs(dist_mut(i,j))+1)) ./ std_mut_diag(abs(dist_mut(i,j))+1) ;
            end
        end
        Zscore_wt(isinf(Zscore_wt))=NaN;
        Zscore_mut(isinf(Zscore_mut))=NaN;  
        figure('Name','Zscore WT without zeros after filters')
            appo=Zscore_wt;
            appo(appo<Zscore_range_filtered(1))=Zscore_range_filtered(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range_filtered(1)-0.1,Zscore_range_filtered(2)])
        figure('Name','Zscore Mut without zeros after filters')
            appo=Zscore_mut;
            appo(appo<Zscore_range_filtered(1))=Zscore_range_filtered(1);
            appo(isnan(appo))=-55;
            imagesc(appo)
            colorbar
            set(gca, 'CLim', [Zscore_range_filtered(1)-0.1,Zscore_range_filtered(2)])
end



% plot Zscore distributions:
% figure('Name', 'distribution of Zscores after filters')
% subplot(2,1,1)    
%     hist(nonzeros(Zscore_wt(:)), 300)
%     xlabel 'Zscore'
%     ylabel 'occurrence'
%     axis([-100 350 0 6e5])
%     legend 'WT'
% subplot(2,1,2)    
%     hist(nonzeros(Zscore_mut(:)), 300)
%     xlabel 'Zscore'
%     ylabel 'occurrence'
%     axis([-100 350 0 6e5])
%     legend 'mut'



%% calculation of Zscore differences
Zscore_diff = Zscore_mut - Zscore_wt;
figure('Name','Zscore difference (Mut - WT)')
    appo=Zscore_diff;
    appo(appo<Zscore_range(1))=Zscore_range(1);
    appo(isnan(appo))=-55;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [Zscore_range(1)-0.1,Zscore_range(2)])
    axis square
    
figure('Name','Zscore difference (Mut - WT)')
    appo=Zscore_diff(209:385,209:385);
    appo(appo<Zscore_range(1))=Zscore_range(1);
    appo(isnan(appo))=-55;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [Zscore_range(1)-0.1,Zscore_range(2)])
    axis square
    
    
    


%% Plot ratio between mut and wt filtered contact maps
figure('Name', 'log2 Ratio of filtered mut/WT maps')
    ratio=log2((mut_data_filtered./dist_mut)./(wt_data_filtered./dist_wt));
    ratio(isinf(ratio))=NaN;
    appo=ratio;
    appo(appo<log_ratio_range(1))=log_ratio_range(1);
    appo(isnan(appo))=-55;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [log_ratio_range(1)-0.1,log_ratio_range(2)])
    axis square
    
figure('Name', 'log2 Ratio of filtered mut/WT maps')
    ratio=log2((mut_data_filtered./dist_mut)./(wt_data_filtered./dist_wt));
    ratio(isinf(ratio))=NaN;
    appo=ratio(209:385,209:385);
    appo(appo<log_ratio_range(1))=log_ratio_range(1);
    appo(isnan(appo))=-55;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [log_ratio_range(1)-0.1,log_ratio_range(2)])    
    axis square
       

    % load WhiteBlueRed colormap:
map3=load('BlueWhiteRed.mat');
% set colormap in figure 4 -> 10
set([ 10 11 13 14 15 16 17 18],'Colormap',map3.BlueWhiteRed)
%set([ 3 4],'Colormap',map3.BlueWhiteRed)
% saveas([15],[mut_sample,'_Zscorediff'],'bmp')
% saveas([16],[mut_sample,'_Zscorediff_zoom'],'bmp')
% saveas([17],[mut_sample,'_log2ratio'],'bmp')
% saveas([18],[mut_sample,'_log2ratio_zoom'],'bmp')




%% top percentage of maximum gain or loss in interaction
map3=load('BlueWhiteRed.mat');
figure('Name', 'log2 Ratio top INCREASE in interaction (GREEN)')
    ratio=log2((mut_data_filtered./dist_mut)./(wt_data_filtered./dist_wt));
    ratio(isinf(ratio))=NaN;
    appo=ratio;
    appo(appo<log_ratio_range(1))=log_ratio_range(1);
    appo(isnan(appo))=-55;
    ts=quantile(appo(appo>-55),1-top/100);
    appo(appo>ts)=max(max(appo))+10;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [log_ratio_range(1)-0.1,log_ratio_range(2)])    
figure('Name', 'log2 Ratio top DECREASE in interaction (GREEN)')
    ratio=log2((mut_data_filtered./dist_mut)./(wt_data_filtered./dist_wt));
    ratio(isinf(ratio))=NaN;
    appo=ratio;
    appo(appo<log_ratio_range(1))=log_ratio_range(1);
    appo(isnan(appo))=-55;
    ts=quantile(appo(appo>-55),top/100);
    appo(appo<ts & appo>-55)=max(max(appo))+10;
    imagesc(appo)
    colorbar
    set(gca, 'CLim', [log_ratio_range(1)-0.1,log_ratio_range(2)])         
set([17,18],'Colormap',map3.BlueWhiteRed)



%% virtual 4C profiling
disp('Extracting virtual 4C profile')
disp(['viewpoint coordinate: ' num2str(viewpoint)])
viewpoint = int64((viewpoint-startcoord)/binsize)+1;
virtual_4C_wt(:,1)=([1:size(wt_data,1)]-1)*binsize+startcoord;
virtual_4C_mut(:,1)=([1:size(mut_data,1)]-1)*binsize+startcoord;
virtual_4C_wt(:,2)=wt_data(:,viewpoint);
virtual_4C_mut(:,2)=mut_data(:,viewpoint);
if (sum(isnan(virtual_4C_mut(:,2)))==size(mut_data,1) | sum(isnan(virtual_4C_wt(:,2)))==size(wt_data,1))
display('CANNOT EXTRACT 4C PROFILE BECAUSE EITHER MUTANT OR WT IS NAN')
else
    figure('Name','virtual_4C')
        plot(virtual_4C_wt(:,1),virtual_4C_wt(:,2))
        hold on
        plot(virtual_4C_mut(:,1),virtual_4C_mut(:,2))
        legend('Wt','Mut')
        xlabel('Genomic Coordinate')
        ylabel('virtual 4C counts')
end
%% interactive linking to Genome Browser
% enter a 'pause' state that can be quit with control+c
while 1
    disp('### Press any key and select interaction to visualize in UCSC')
    disp('### Press Ctrl+c to quit')
    disp('-------------------------------------------------------------')
    pause('on');
    pause;
    [x y]=(ginput(1));
    x=round(x);
    y=round(y);
    coordx = gencoord(x);
    coordy = gencoord(y);
 
    %ucsc_str = ['chr10:',num2str(us_bin_coord),'-',num2str(ds_bin_coord)]
    if (strcmp(mut_info(3),'Inversion') & strcmp(Inversion,'true'))
        disp('ATTENTION! It modifies the genomic region according to the inversion!')
        disp('If you want to display the wt maps with wt coordinate, set Inversion=false')
        inv_start = str2double(mut_info(7));
        inv_end = str2double(mut_info(8));
        if(coordx>inv_start & coordx<inv_end)
           coordx=inv_end-(coordx-inv_start); 
        end
        if(coordy>inv_start & coordy<inv_end)
           coordy=inv_end-(coordy-inv_start); 
        end
        if(coordy>coordx)
            ucsc_str = ['https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A',num2str(coordx),'-',num2str(coordy),'&hgsid=390798063_VxbwBzrDXaAtOM7WArKjA4whs0Rs', '-browser'];
        else
            ucsc_str = ['https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A',num2str(coordy),'-',num2str(coordx),'&hgsid=390798063_VxbwBzrDXaAtOM7WArKjA4whs0Rs', '-browser'];
        end
    else
    ucsc_str = ['https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A',num2str(coordx),'-',num2str(coordy),'&hgsid=390798063_VxbwBzrDXaAtOM7WArKjA4whs0Rs', '-browser'];    
    end
    web(ucsc_str, '-browser')
end
