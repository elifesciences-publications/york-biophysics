function [nPositionsInCluster, distancesWithinClusters, distancesBetweenClusters meanClusterPos] = nearestNeighbourClustering(positions, clusterThreshold)
% Nearest neighbour based clustering of positions.
% by Stephan Uphoff 12.09.2012
%
% positions = N x 2 matrix of [x y] coordinates for N localisations
% clusterThreshold = positions separated by less than this threshold are
% assigned to the same cluster.

nPos = length(positions(:,1));

remainingPosIndices = (1:nPos)'; % give each position an index

posClusterNumber = zeros(size(positions(:,1))); % will be used to assign each position to a cluster

% initialise
clusterNumber = 0;

figure; hold all;

while numel(remainingPosIndices) > 0
    
    clusterNumber = clusterNumber + 1;

    clusterIndices = remainingPosIndices(1);
    
    while numel(clusterIndices) > 0
    
    % choose currentClusterIndex as the first in the list of clusterIndices
    % it will be removed after this round so the next lowest index will be
    % at clusterIndices(1)
    currentClusterIndex = clusterIndices(1);
    
    distances = nan(nPos,1);
    
    distances(remainingPosIndices) = sqrt((positions(remainingPosIndices,1) - positions(currentClusterIndex,1)).^2 +...
        (positions(remainingPosIndices,2) - positions(currentClusterIndex,2)).^2);
    
    clusterIndices = [clusterIndices; find(distances < clusterThreshold)];
    
    %keep only unique indices
    clusterIndices = unique(clusterIndices);
    
    %remove the currentClusterIndex from clusterIndices so it is not
    %included again
    clusterIndices(find(clusterIndices == currentClusterIndex)) = [];
    
    %remove currentClusterIndex from clusterIndices
    remainingPosIndices(find(remainingPosIndices == currentClusterIndex)) = [];
    
    %label currentClusterIndex with clusterNumber
    posClusterNumber(currentClusterIndex) = clusterNumber;
    
    end
    
end


%%% analyse clusters

nPositionsInCluster = zeros(clusterNumber-1, 1);
distancesWithinClusters = [];
meanClusterPos = zeros(clusterNumber-1, 2);

% loop over all clusters
for ii = 1:clusterNumber
    
    % plot localisations colour coded according to their cluster assignment
    plotcolor = rand(3,1);
    
    index = find(posClusterNumber == ii);
    plot(positions(index,1), positions(index,2), 'o', 'Color', plotcolor, 'MarkerFaceColor', plotcolor)
    
    % number of positions in this cluster
    nPositionsInCluster(ii) = length(index);
    
    % calculate distances between positions within this cluster
    for jj = 1:length(index) 
        
       dist = sqrt((positions(index,1) - positions(index(jj),1)).^2 +...
        (positions(index,2) - positions(index(jj),2)).^2);
    
       % concatenate all distances (excluding dist = 0 between the position
       % and itself)
       distancesWithinClusters = [distancesWithinClusters; dist(dist>0)];
       
    end
    
    % center position of this cluster
    meanClusterPos(ii,:) = mean(positions(index,:),1);
    
end

% remove the pair-wise double counting of distances by keeping only unique
% distances
distancesWithinClusters = unique(distancesWithinClusters);

% calculate distances between cluster centers
distancesBetweenClusters = [];

for kk = 1:clusterNumber % calculate inter cluster distances
        
       dist = sqrt((meanClusterPos(kk,1) - meanClusterPos(:,1)).^2 +...
        (meanClusterPos(kk,2) - meanClusterPos(:,2)).^2);
    
       % concatenate all distances (excluding dist = 0 between the position
       % and itself)
       distancesBetweenClusters = [distancesBetweenClusters; dist(dist>0)];
       
end

% remove the pair-wise double counting of distances by keeping only unique
% distances
distancesBetweenClusters = unique(distancesBetweenClusters);

% plot all positions for reference
plot(positions(:,1), positions(:,2), 'kx');
hold off


figure; hist(distancesWithinClusters,100);
xlabel('distance')
title('Distances within clusters');

figure; hist(distancesBetweenClusters,100);
xlabel('distance')
title('Distances between clusters');

end