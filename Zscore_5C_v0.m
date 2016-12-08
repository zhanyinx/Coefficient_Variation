clear all
clc
close all

% options:
binsize = 5000; % genomic size of the bin
startcoord = 99001149; % start coordinate of the first bin in the square matrix (after cutting with convert.sh, see below)
ZEROS = 'false'; % if 'true', Zscores are calculated keeping 0s in the data; otherwise if 'false' they are discarded
 % size of the deletion (useful for calculating distance matrix!)
deletion_start = 100621959;
deletion_end = 100640125;
deletion_size = deletion_start-deletion_end; 

% determine start and end bin of the deletion, excluding both start and end bins:
deletion_start_bin = fix((deletion_start -startcoord) / binsize +1)-2;
deletion_end_bin = ceil((deletion_end -startcoord) / binsize +1)-2;
deletion_size_bin = deletion_end_bin - deletion_start_bin;

% file names of WT and mutant sample (these are Nicolas' pirwise matrices, wither raw or ICEd):
wt_filename = 'B20_E14-2_wo_20kb_4_iced.mat';
mut_filename = 'B29_C3_wo_20kb_4_iced.mat';

% INSERT CHECHPOINT to see if there is already a matrix file and avoid
% conversion:

% convert to matrix format and cut the matrices to ensure that they have the same size
% using Zhan's 'conver.sh' Bash script
% unix(['./convert.sh ', wt_filename])
% unix(['./convert.sh ', mut_filename])


% import WT 5C square matrix
wt = importdata([wt_filename, '.matrix'], ' ');
% import mutant 5C square matrix
mut = importdata([mut_filename, '.matrix'], ' ');

% extract matrix entries:
wt_data = wt.data; % contains the 'data' field of the structure wt
mut_data = mut.data; % idem

% genomic coordinates of beginning of bins (assuming 5kb binning):
gencoord = ([1:size(wt_data,1)]-1) * binsize + startcoord;

% normalize by the total number of reads and by a common read number:
wt_data = wt_data ./ sum(sum(wt_data)) *10e6;
mut_data = mut_data ./ sum(sum(mut_data)) *10e6;


% plot WT and mutant matrices 
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




%% Zscore calculation

% build a matrix of distances for WT sample
for i=1:size(wt_data,1)
    for j=1:size(wt_data,1)
        dist_wt(i,j) = i-j;
    end
end
% modify parts affected by deletion in the mutant:
dist_mut = dist_wt;
for i = 1:deletion_start_bin-1
    for j = deletion_end_bin:size(wt_data,1)
        dist_mut(i,j)=dist_wt(i,j)+deletion_size_bin;
        dist_mut(j,i)=dist_wt(j,i)-deletion_size_bin;
    end
end
% to make sure that distances within the deletion
% dist_mut(deletion_start_bin:deletion_end_bin-1,:)=55555;
% dist_mut(:,deletion_start_bin:deletion_end_bin-1)=55555;

% dist_mut=dist_wt;


% extract matrix elements at fixed distances
for i=1:size(wt_data,1)
    wt_diag{i} = wt_data(dist_wt==i);
    mut_diag{i} = mut_data(dist_mut==i);
end


switch ZEROS
    
    case 'true' % keeping zeros
        % calculate means over diagonals:
        for i=1:size(wt_data,1)
            mean_wt_diag(i) = mean(wt_diag{i});
            mean_mut_diag(i) = mean(mut_diag{i});
        end
        % calculate Zscores
        for i=1:size(wt_data,1)
            for j=1:size(wt_data,1)
                Zscore_wt(i,j) = wt_data(i,j) ./ mean_wt_diag(abs(dist_wt(i,j))+1);
                Zscore_mut(i,j) = mut_data(i,j) ./ mean_mut_diag(abs(dist_mut(i,j))+1);

            end
        end
        figure('Name','Zscore WT with zeros')
            imagesc(Zscore_wt)
            colorbar
            set(gca, 'CLim', [0, 3])
        figure('Name','Zscore Mut with zeros')
            imagesc(Zscore_mut)
            colorbar
            set(gca, 'CLim', [0, 3])
    
        
    case 'false' % excluding zeros
        % calculate means over diagonals:
        for i=1:size(wt_data,1)
            mean_wt_diag(i) = mean(nonzeros(wt_diag{i}));
            mean_mut_diag(i) = mean(nonzeros(mut_diag{i}));
        end
        % calculate Zscores
        for i=1:size(wt_data,1)
            for j=1:size(wt_data,1)
                Zscore_wt(i,j) = wt_data(i,j) ./ mean_wt_diag(abs(dist_wt(i,j))+1);
                Zscore_mut(i,j) = mut_data(i,j) ./ mean_mut_diag(abs(dist_mut(i,j))+1);
            end
        end
        figure('Name','Zscore WT without zeros')
            imagesc(Zscore_wt)
            colorbar
            set(gca, 'CLim', [0, 3])
        figure('Name','Zscore Mut without zeros')
            imagesc(Zscore_mut)
            colorbar
            set(gca, 'CLim', [0, 3])
end

%% calculation of Zscore differences
Zscore_diff = Zscore_mut - Zscore_wt;
figure('Name','Zscore difference (Mut - WT)')
    imagesc(Zscore_diff)
    colorbar
    set(gca, 'CLim', [-1.5, 1.5])
% load WhiteBlueRed colormap:
map3=load('BlueWhiteRed.mat');
% set colormap in figure 1 and 2
set(5,'Colormap',map3.BlueWhiteRed)

% filter singletons in Zscore map as in Hnisz et al but based on histogram of nonzero Zscore diffs:
Zscore_diff(Zscore_diff>10)


% build vectors containing diagonal elements of the Z-score difference matrix:
for i=1:size(wt_data,1)
    Zscore_diag{i} = Zscore_diff(dist_wt==i);
    mean_Zscore_diag(i) = mean(nonzeros(Zscore_diag{i}));
    kurt_Zscore_diag(i) = kurtosis(nonzeros(Zscore_diag{i}));
end

