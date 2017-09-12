function msd = get_MSD(tracks)
for i=1:length(tracks)
  temp2 = zeros(length(tracks(i).time),1);
    for n=0:length(tracks(i).time)-1 % timei ndexes
      
    N = length(tracks(i).time);

    for j=1:N-n
       temp2(n+1)=temp2(n+1)+(tracks(i).xvalues(j+n)-tracks(i).xvalues(j)).^2+(tracks(i).yvalues(j+n)-tracks(i).yvalues(j)).^2;
    end
    temp2(n+1) = temp2(n+1)/(N-n);
    end

    msd = temp2;
    
end 