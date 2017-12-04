function [fibers_out,SuperFiber_out] = FlipFibers(fibers_in,SuperFiber_in)
SuperFiber_out.mean = fliplr(SuperFiber_in.mean);
SuperFiber_out.spread = fliplr(SuperFiber_in.spread);

for i=1:length(fibers_in)
    fibers_out{i} = fliplr(fibers_in{i});
end
fibers_out = fibers_out';
