function is_flipped_orientation = CheckOrientationWrtRef(SuperFiber, SuperFiberRef)

curve1 = SuperFiberRef;
curve2 = SuperFiber;
curve2flipped =fliplr(SuperFiber);
[noflipdist, ~] = DiscreteFrechetDist(curve1,curve2);%tanya  - changed it since the other measure doesn't work well for ring-shaped tracts
[flipdist,~] = DiscreteFrechetDist(curve1,curve2flipped);
if flipdist<noflipdist
    is_flipped_orientation = 1;
else
    is_flipped_orientation = 0;
end

