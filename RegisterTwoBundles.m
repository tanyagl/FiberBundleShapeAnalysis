function [ fibers_moved_step2, fibers_moved_step1, fibers_moving_rt, fibers_moving_zc, fibers_moving_zc_no_flip, fibers_fixed_zc ] = RegisterTwoBundles( fg_moving, fg_fixed, numNodes, mirror_right_to_left )
% This function registers fg_moving fiber bundle to fg_fixed fiber bundle.
% INPUTS:
%       fg_moving - the bundle to be aligned
%       fg_fixed - the reference bundle we're alighning to
%       numNodes - number of nodes to which to resample each fiber
%       mirror_right_to_left - when the registration is between left and
%       right homologous fiber bundles (e.g. for a brain lateralization
%       study), we first mirror-flip the right bundle.
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

    [resampled_fibers_ref, ~, ~, ~] = dtiReorientFibers(fg_fixed, numNodes);
    [resampled_fibers_moving, ~, ~, ~] = dtiReorientFibers(fg_moving, numNodes);

    %%% clean up the bundle
    resampled_fibers_ref_cleaned = CleanUpFibers(resampled_fibers_ref.fibers);
    resampled_fibers_moving_cleaned = CleanUpFibers(resampled_fibers_moving.fibers);
    %%% zerocenter all fibers
    fibers_fixed_zc = ZeroCenterFibers( resampled_fibers_ref_cleaned);
    fibers_moving_zc = ZeroCenterFibers( resampled_fibers_moving_cleaned);
    %%% compute SuperFibers
    [SuperFiber_ref_zc,~,~] = ComputeSuperFiber(fibers_fixed_zc,numNodes);
    [SuperFiber_moving_zc,~,~] = ComputeSuperFiber(fibers_moving_zc,numNodes);
    %%%flip if necessary
    is_flipped_orientation = CheckOrientationWrtRef(SuperFiber_moving_zc.mean,SuperFiber_ref_zc.mean);
    if(is_flipped_orientation)
        [fibers_moving_zc,SuperFiber_moving_zc] = FlipFibers(fibers_moving_zc,SuperFiber_moving_zc);
    end
    
    if( mirror_right_to_left)
        fibers_moving_zc_no_flip = fibers_moving_zc;
        for fiber = 1:length(fibers_right_zc)
            fibers_moving_zc{fiber} = (fibers_moving_zc{fiber}.*[-1;1;1]);
        end
    end
    %%% calculate clustering
    RefData = [];
    [~, seed_points_ref_zc, RefData] = SimpleClustering(fibers_fixed_zc,SuperFiber_ref_zc,RefData,1,3,1);
    %    [~, seed_points2_zc{subj}, ~] = SimpleClustering(fibers2_zc{subj},SuperFiber2_zc{subj},RefData{subj},0,3,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% DEFORMATION ESTIMATION based on registration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% STEP 1 - ALIGN THE SUPERFIBERS %%%%%%%%%%%%%
    % TBD: limit rotation within some angle range - otherwise for some short
    % bundles the result is >180 rotation
    [R,T] = rigid_transform_3D(SuperFiber_moving_zc.mean',SuperFiber_ref_zc.mean');
    fibers_moving_rt = (ApplyRigidTransform(fibers_moving_zc,T, R))';

    %seed_points2_rt{subj} = R{subj}* seed_points2_zc{subj}'+T{subj};
    %%% need to recompute SuperFiber mean again for the zc tract
    SuperFiber_moving_rt.mean = R*SuperFiber_moving_zc.mean+T;
    %%% 'p' stands for positive. The moving points have to be all
    %%% in the positive coordinates since the grid is calculated on
    %%% the positive coordinates. The grid has to be tight around
    %%% the boundaries of the moving points.
    xMoving_p = SuperFiber_moving_rt.mean' + abs(min(SuperFiber_moving_rt.mean'));
    xStatic_p = SuperFiber_ref_zc.mean' + abs(min(SuperFiber_ref_zc.mean'));

    [p,h] = size(fibers_moving_rt); 
    all_points_rt = [];
    for fiber = 1:max(p,h)
        all_points_rt = [all_points_rt, fibers_moving_rt{fiber}];
    end
    cuboid_dim_s1 = abs( max(all_points_rt') - min(all_points_rt'));
    sizeI = ceil(abs( (cuboid_dim_s1)));
    option=struct; options.MaxRef=6;
    [O_trans1,Spacing,xReg] = point_registration(sizeI,xMoving_p,xStatic_p,options);
    aligned_fiber = {};fibers_moved_step1 = {};
    for fiber = 1:size(fibers_moving_rt,1)
        translation = min(fibers_moving_rt{fiber}');
        fiber_p = fibers_moving_rt{fiber}' - translation;
        aligned_fiber{fiber} = (bspline_trans_points_double(O_trans1,Spacing,fiber_p))';% - translation';
        aligned_fiber{fiber} = dtiFiberResample(aligned_fiber{fiber}, numNodes);
        fibers_moved_step1{fiber} = (movingmean(aligned_fiber{fiber}',3)');
    end
    % now that we aligned and smoothed the fibers according to the
    % SuperFiber alignment, recompute the SuperFiber for the aligned bundle
    % and rigidly align the bundles again
    [SuperFiber_moving_s1,~,~] = ComputeSuperFiber(fibers_moved_step1,numNodes);
    [R1,T1] = rigid_transform_3D(SuperFiber_moving_s1.mean',SuperFiber_ref_zc.mean');
    fibers_moved_step1 = (ApplyRigidTransform(fibers_moved_step1,T1,R1));
    [SuperFiber_moving_s1,~,~] = ComputeSuperFiber(fibers_moved_step1,numNodes);
    %now the bundle has been coarsely moved and aligned according to the
    %SuperFiber, need to compute the clustering for wireframe alignment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% STEP 2 - ALIGN THE WIREFRAME %%%%%%%%%%%%%

    wireframepoints = 21:340;%10:2:153;%7:2:85;
    %[~, seed_points1_s1, RefData] = SimpleClustering(fibers1_zc,SuperFiber1_zc,RefData,1,3,1);
    [~, seed_points_moving_s1, ~] = SimpleClustering(fibers_moved_step1,SuperFiber_moving_s1,RefData,0,3,1);
    xMoving_p = seed_points_moving_s1(wireframepoints,:)+ abs(min(seed_points_moving_s1(wireframepoints,:)));
    xStatic_p = seed_points_ref_zc(wireframepoints,:)+ abs(min(seed_points_ref_zc(wireframepoints,:)));
    all_points_s1 = [];
    [p,h] = size(fibers_moved_step1); 

    for fiber = 1:max(p,h)
        all_points_s1 = [all_points_s1, fibers_moved_step1{fiber}];
    end
    %cuboid_dim_s2 = abs( max(all_points_s1') - min(all_points_s1'));
    cuboid_dim_s2 = max(abs( max(all_points_s1') - min(all_points_s1')),abs( max(seed_points_moving_s1(wireframepoints,:)) - min(seed_points_ref_zc(wireframepoints,:))));
    sizeI = ceil(abs( (cuboid_dim_s2)));
    option=struct; options.MaxRef=6;
    tic
    [O_trans,Spacing,xReg] = point_registration(sizeI,xMoving_p,xStatic_p,options);
    toc

    aligned_fiber_s2 = {};fibers_moved_step2 = {};
    [p,h] = size(fibers_moved_step1); 

    for fiber = 1:max(p,h)
        translation = abs(min(fibers_moved_step1{fiber}'));
        fiber_p = fibers_moved_step1{fiber}' + translation;
        %                 translation = min(aligned_fiber_resampled{fiber}');
        %                 fiber_p = aligned_fiber_resampled{fiber}' - translation;
        aligned_fiber_s2{fiber} = (bspline_trans_points_double(O_trans,Spacing,fiber_p))';% - translation';
        fibers_moved_step2{fiber} = dtiFiberResample(aligned_fiber_s2{fiber}, numNodes);
        fibers_moved_step2{fiber} = (movingmean(fibers_moved_step2{fiber}',3))';
    end

    %post processing required due to translation and to make sure
    %they are aligned
    [SuperFiber_moving_s2,~,~] = ComputeSuperFiber(fibers_moved_step2,numNodes);
    [R1,T1] = rigid_transform_3D(SuperFiber_moving_s2.mean',SuperFiber_ref_zc.mean');
    aligned_fiber_s2 = (ApplyRigidTransform(aligned_fiber_s2,T1,R1));
    fibers_moved_step2 = (ApplyRigidTransform(fibers_moved_step2,T1,R1));
    %     [SuperFiber_moving_s2,~,~] = ComputeSuperFiber(aligned_fiber_resampled_s2,numNodes);
    %     [~, seed_points_moving_s2, ~] = SimpleClustering(aligned_fiber_resampled_s2,SuperFiber_moving_s2,RefData,0,3,1);
end
