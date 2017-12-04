function [ d_m, d_s ] = CalculateDeformation(fg_moving,fg_moved, numNodes )
% This function registers fg_moving fiber bundle to fg_fixed fiber bundle.
% INPUTS:
%       fg_moving - the bundle before registration
%       fg_moved - the bundle after registration
%       numNodes - number of nodes to which to resample each fiber
% OUTPUTS:
%       fibers_moved_step2 - the output of step 2 of registration (wireframe
%                            points alignment)
%       fibers_moved_step1 - the output of step 2 of registration
%                            (SuperFiber points alignment)
%       fibers_moving_rt - moving bundle after rigid transformation to
%                          align with the fixed bundle
%       fibers_moving_zc - moving bundle zero-centered (and flipped if
%                          mirror_right_to_left flag is on)
%       fibers_moving_zc_no_flip -  moving bundle zero-centered (not flipped)
%       fibers_fixed_zc - fixed bundle after zero centered
%
% written by Tanya Glozman, Stanford University, 2017


% % How to calculate deformation: 
% %---------------------------------
% % we calculate 3 deformation measures: one for each step of the
% % registration and the total deformation. Since the number of nodes or the
% % order of nodes on the fibers does not change during the deformation, to
% % calculate the deformation for each node, we simply calculate the norm
% % between the node location before the deformation and its lcoation
% % thereafter.
% % To normalize the deformation we calculate the arclength of the SuperFiber
% % at each step and that is used as a normalization factor, so that a
% % % deformation is relative
% The deformation we calculate should be relative to the original bundle
% dimensions
% we calculate the deformation normalization factor as follows: we find the
% length of the SuperFiber curve using arclen and we find the average area
% of a cross section using the surfaceArea of an alphaShape surounding the
% points. Mupliplying these yields an estimate of the bundle volume. This
% is ouor normalization factor. Actually the alphaShape surface area is not
% a very good estimate because it isn't very consistent with age - because
% of the template problems etc. which probably causes us to not see all the
% fibers.
  %%% compute SuperFibers
[SuperFiber_moved,~,~] = ComputeSuperFiber(fg_moved,numNodes);
[SuperFiber_moving,~,~] = ComputeSuperFiber(fg_moving,numNodes);
[arclen_moving,~] = arclength(SuperFiber_moving.mean(1,:),SuperFiber_moving.mean(2,:),SuperFiber_moving.mean(3,:));                                
%[arclen_moved,~] = arclength(SuperFiber_moved.mean(1,:),SuperFiber_moved.mean(2,:),SuperFiber_moved.mean(3,:));
cnt = 0;
for fiber = 1:length(fg_moved)
    for point_in_fiber = 1:numNodes
        cnt = cnt+1;
        deformation_factor(cnt) = norm(fg_moving{fiber}(:,point_in_fiber) - fg_moved{fiber}(:,point_in_fiber));
 %       displacement_vector{fiber}(point_in_fiber,:) = fg_moved{fiber}(:,point_in_fiber) - fg_moving{fiber}(:,point_in_fiber);
    end
end
d = reshape(deformation_factor,numNodes,length(deformation_factor)/numNodes);
d_m = mean(d(:))/arclen_moving; %normalizing factor - notice this is reasonable for a rather long bundle, not necessarily for a short bundle
d_s = std(d(:))/arclen_moving;
end

