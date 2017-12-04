function transformed_fibers = ApplyRigidTransform(fibers, T, R)
for fiber = 1:length(fibers)
    transformed_fibers{fiber} = (R*fibers{fiber} + T);%*ones(1,length(fibers{fiber}));
end
