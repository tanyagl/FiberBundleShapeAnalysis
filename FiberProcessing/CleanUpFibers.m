function fibers_cleaned = CleanUpFibers(fibers)
% remove the fibers that are too short or
% too long compared to the average arclength of all fibers
    for i=1:length(fibers)
        fiber = fibers{i};
        a(i) = arclength(fiber(1,:),fiber(2,:),fiber(3,:));
    end
    m = median(a);s = std(a);
    keep = ((a >= m-4*s)& (a <= m+4*s))';
%     if(sum(keep)<length(fibers))
%         stop_here = 1;
%     end
    fibers_cleaned = fibers(keep);
  