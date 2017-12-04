function fibers_zc = ZeroCenterFibers( fibers)
   centroid = GetCentroid( fibers);
   for fiber = 1:length(fibers)
       fibers_zc{fiber} = fibers{fiber} - centroid';
   end
