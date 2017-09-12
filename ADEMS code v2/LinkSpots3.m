function [spot_traj]=LinkSpots3(spots, framenum1, framenum0, d_01_max, Iratio_01_min,...
    Iratio_01_max, SigmaRatio_01_min, SigmaRatio_01_max)

if isempty(spots(spots(:,9)==framenum0))==0 && ... % at least 1 accepted spot in previous frame and
        isempty(spots(spots(:,9)==framenum1))==0  % at least 1 accepted spot in current frame.
    
    N1 = size(spots(spots(:,9)==framenum1),1);  % no. of accepted spots in current frame
    N0 = size(spots(spots(:,9)==framenum0),1);  % no. of accepted spots in previous frame
    
    %initialise a matrix of spot pairs
    Assignments=ones(N0,N1,5)*inf;
    %counters=bad programming
    m0=1;
    m1=1;
    % loop though accepted spots in previous frame.
    for q0 = min(find(spots(:,9)==framenum0)):max(find(spots(:,9)==framenum0))
        m1=1;
        for q1 = min(find(spots(:,9)==framenum1)):max(find(spots(:,9)==framenum1))
            % d_01: distance between spot Centres in previous and current frames:
            d_01 = sqrt((spots(q0,1)-spots(q1,1))^2+(spots(q0,2)-spots(q1,2))^2);
            % Iratio_01: ratio of intensities of spot Centre in previous and current frames:
            Iratio_01= spots(q0,5)/spots(q1,5);
            % SigmaRatio_01: ratio of widths of spots (Gaussian fits) in previous and current frames AVERAGE X AND Y:
            SigmaRatio_01= (spots(q0,6)+spots(q0,7))/(spots(q1,6)+spots(q1,7));
            
            %If the criteria are met, save the values and spot numbers to
            %the assignment matrix
            if d_01 < d_01_max && ...  % see PARAMETERS at start of this function.
                    Iratio_01_min <= Iratio_01 && Iratio_01 <= Iratio_01_max && ...
                    SigmaRatio_01_min <= SigmaRatio_01 && SigmaRatio_01 <= SigmaRatio_01_max
                Assignments(m0,m1,1)=d_01;
                Assignments(m0,m1,2)=Iratio_01;
                Assignments(m0,m1,3)=SigmaRatio_01;
                Assignments(m0,m1,4)=q0;
                Assignments(m0,m1,5)=q1;
            end
            %Arbitrarily remove duplicate assigments here MUST BE A BETTER WAY
            Assignments(Assignments(:,m1,1)>min(Assignments(:,m1,1)),m1,:)=inf;
            m1=m1+1;
        end
        %Arbitrarily remove duplicate assigments here MUST BE A BETTER WAY
        Assignments(m0,Assignments(m0,:,1)>min(Assignments(m0,:,1)),:)=inf;
        m0=m0+1;
    end
    
    [a, b]=find(Assignments(:,:,1)<inf);
    %If there are some successful assignments
    if isempty(a)==0
        %Find current largest trajectory number
        tr=max(spots(:,10));
        % Loop over all the assignments
        for i=1:size(a)
            
            %If the spot in the previous frame had an assignment, use that
            %trajectory number
            if spots(Assignments(a(i),b(i),4),10)>0
                spots(Assignments(a(i),b(i),5),10)=spots(Assignments(a(i),b(i),4),10);
            else
                tr=tr+1;
                %If new track, assign new trajectory number
                spots(Assignments(a(i),b(i),4),10)=tr;
                spots(Assignments(a(i),b(i),5),10)=tr;
                
            end
           % tr=tr+1;
        end
    end
end
spot_traj=spots;
end