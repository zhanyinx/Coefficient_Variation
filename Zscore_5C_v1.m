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
binsize = 5000; % genomic size of the bin
startcoord = 99011149; % start coordinate of the first bin in the square matrix (after cutting with convert.sh, see below)
ZEROS = 'true'; % if 'true', Zscores are calculated keeping 0s in the data; otherwise if 'false' they are discarded
tablename = '20160822_5C-Samples.xlsx'; % Xls file with information on samples
%% 1. data input
% samples:
wt_sample = 'B20_E14-2';
mut_sample = 'B29_51.13';
disp(['** Hello Rafael, I am processing samples ', wt_sample, ' and ', mut_sample, ' **'])

% look up if the mutant is a deletion or inversion and find deletion coordinates
table=readtable(tablename); % imports the database table with all 5C samples
disp(['Importing sample information from table ', tablename])
table(:,'Var11') = []; % get rid of 11th column (useless)
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
     deletion_end_bin = ceil((deletion_end -startcoord) / binsize +1);
     deletion_size_bin = deletion_end_bin - deletion_start_bin;
     disp(['Corresponding to start bin ', num2str(deletion_start_bin)])
     disp(['and end bin ', num2str(deletion_end_bin)])

end

% file names of WT and mutant sample (these are Nicolas' pirwise matrices, wither raw or ICEd):
wt_filename = [wt_sample, '_wo_20kb_4_raw.mat'];
mut_filename = [mut_sample, '_wo_20kb_4_raw.mat'];

if strcmp('raw', wt_filename)
    disp('Working on raw data')
elseif strcmp('ice', wt_filename)
    disp('Working on raw data')
end


%% Pairwise matrix conversion to square matrix format, cutting to uniform region
%     convert to matrix format and cut the matrices to ensure that they have the same size
%     using Zhan's 'conver.sh' Bash script
disp('Converting pairwise format into square matrix...')
unix(['./convert.sh ', wt_filename])
unix(['./convert.sh ', mut_filename])

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
row_sum = sum(wt_data,2);
% assign rows where sum of counts is striclty 0 to NaN (non-working primers)
wt_data(row_sum==0,:)=NaN;
mut_data(row_sum==0,:)=NaN;
wt_data(:,row_sum==0)=NaN;
mut_data(:,row_sum==0)=NaN;
disp('Assigning masked bins to Nan')

%% plot WT and mutant matrices 
figure('Name', 'WT')
    imagesc(wt_data)
    colorbar
    set(gca, 'CLim', [0, 200])

figure('Name', 'mut')
    imagesc(mut_data)
    colorbar
    set(gca, 'CLim', [0, 200])

% load WhiteRedBlack colormap:
map3=load('WhiteRedBlack');
% set colormap in figure 1 and 2
set(1:2,'Colormap',map3.WhiteRedBlack)

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
        figure('Name','Zscore WT with zeros')
            imagesc(Zscore_wt)
            colorbar
            set(gca, 'CLim', [-3, 3])
        figure('Name','Zscore Mut with zeros')
            imagesc(Zscore_mut)
            colorbar
            set(gca, 'CLim', [-3, 3])
    
        
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
        figure('Name','Zscore WT without zeros')
            imagesc(Zscore_wt)
            colorbar
            set(gca, 'CLim', [-3, 3])
        figure('Name','Zscore Mut without zeros')
            imagesc(Zscore_mut)
            colorbar
            set(gca, 'CLim', [-3, 3])
end



% plot Zscore distributions:
figure('Name', 'distribution of Zscores')
subplot(2,1,1)    
    hist(nonzeros(Zscore_wt(:)), 300)
    xlabel 'Zscore'
    ylabel 'occurrence'
    axis([-100 350 0 6e5])
    legend 'WT'
subplot(2,1,2)    
    hist(nonzeros(Zscore_mut(:)), 300)
    xlabel 'Zscore'
    ylabel 'occurrence'
    axis([-100 350 0 6e5])
    legend 'mut'

% filter out singletons:
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
        figure('Name','Zscore WT with zeros after filters')
            imagesc(Zscore_wt)
            colorbar
            set(gca, 'CLim', [-3, 3])
        figure('Name','Zscore Mut with zeros after filters')
            imagesc(Zscore_mut)
            colorbar
            set(gca, 'CLim', [-3, 3])
    
        
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
        figure('Name','Zscore WT without zeros after filters')
            imagesc(Zscore_wt)
            colorbar
            set(gca, 'CLim', [-3, 3])
        figure('Name','Zscore Mut without zeros after filters')
            imagesc(Zscore_mut)
            colorbar
            set(gca, 'CLim', [-3, 3])
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
    imagesc(Zscore_diff)
    colorbar
    set(gca, 'CLim', [-2, 2])
% load WhiteBlueRed colormap:
map3=load('BlueWhiteRed.mat');
% set colormap in figure 1 and 2
set(5,'Colormap',map3.BlueWhiteRed)

% % filter singletons in Zscore map as in Hnisz et al but based on histogram of nonzero Zscore diffs:
% Zscore_diff(Zscore_diff>10)
% 
% 
% % build vectors containing diagonal elements of the Z-score difference matrix:
% for i=1:size(wt_data,1)
%     Zscore_diag{i} = Zscore_diff(dist_wt==i);
%     mean_Zscore_diag(i) = mean(nonzeros(Zscore_diag{i}));
%     kurt_Zscore_diag(i) = kurtosis(nonzeros(Zscore_diag{i}));
% end


%% Plot ratio between mut and wt filtered contact maps
figure('Name', 'log2 Ratio of filtered mut/WT maps')
    imagesc(log2(mut_data_filtered./wt_data_filtered))
    colorbar

set([3 4 6 7 8 9],'Colormap',map3.BlueWhiteRed)
set(gca, 'CLim', [-2, 2])

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
    ucsc_str = ['https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm9&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chrX%3A',num2str(coordx),'-',num2str(coordy),'&hgsid=390798063_VxbwBzrDXaAtOM7WArKjA4whs0Rs', '-browser'];
    web(ucsc_str, '-browser')
end
