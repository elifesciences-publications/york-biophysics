%%DistanceCluster
% Links n dimensional data of m points into clusters within linkDist and returns a
% vector of m points with assigned cluster numbers.
%Input Mig1coords must be n x m but can have any dimensionality


function [ClusterMatrix]=DistanceCluster(Mig1coords,linkDist)
% linkDist=30;
%  Mig1coords=cat(1,Xmig1,Ymig1,Zmig1);
PWdist=pdist(Mig1coords');
PWdistM=squareform(PWdist);
ClusterMatrix=zeros(size(PWdistM,1),1);
ClusterNumber=1;
for i=1:size(PWdistM,1)
    for j=1:size(PWdistM,2)
        if i==j
        else
            if PWdistM(i,j)<linkDist
                if ClusterMatrix(i,1)+ClusterMatrix(j,1)==0
                    ClusterMatrix(i,1)=ClusterNumber;
                    ClusterMatrix(j,1)=ClusterNumber;
                    ClusterNumber=ClusterNumber+1;
                else
                    PossCN=cat(1,ClusterMatrix(i,1),ClusterMatrix(j,1));
                    ClusterMatrix(i,1)=min(PossCN(PossCN>0));
                    ClusterMatrix(j,1)=min(PossCN(PossCN>0));
                end
            end
        end
    end
end
