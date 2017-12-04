function [IDX, seed_points, ReferenceData] = SimpleClustering(resampled_fibers,SuperFiber,RefData,is_ref, stride, step)
%load 'C:\MyWork\Stanford\data\data.mat';
%load 'C:\MyWork\Stanford\data\RGdeformation.mat';
%note: here we assume that all fibers were resampled to 100 points. If this
%isn't so the following parameters need updating:axial, direction_point_gap
num_seeds_per_slab = 16;
ReferenceData = {};
IDX = {};
num_points = size(resampled_fibers{1},2);
%stride = 3;%ceil(num_points/10);
axial = [2:stride:num_points-1];
%step = 1;

%axial = [1,10*(1:10)]; %we're taking 11 points along the superfiber - these point are called the superfiber seeds. at each such point there will be a 'cake'
direction_point_gap = 2;%ceil(stride/2); %this is the length of the vector along the superfiber that we use to determine the direction of the fiber at each superfiberseed point
%update this parameter to get more localized direction
c = SuperFiber.mean(:,axial)'; %superfiber seed points - centerpoints along the superfiber
for i = 1:length(axial)
    seeds = [];
 
    lb = axial(i)-step;%lower bound
    if lb<1
        lb = 1;
    end
    ub = axial(i)+step;%upper bound
    if ub>length(SuperFiber.mean)
        ub = length(SuperFiber.mean);
    end
    SF_direction(i,:) = SuperFiber.mean(:,ub) - SuperFiber.mean(:,lb);
    SF_direction(i,:) = SF_direction(i,:)./norm(SF_direction(i,:));
    
    [p,p_projected] = GetPointsAround(axial,i,resampled_fibers,SuperFiber.mean,SF_direction(i,:),lb,ub);
    
    % p contains all the points in the fibers that correspond to the specific
    % location on the Superfiber. The larger eigenvectors of the covariance
    % matrix of these points points at the direction of the greatest
    % variance (major axis of the "ellipse" formed by the points). The
    % 2nd eigenvector points at the direction of the minor axis. The
    % corresponding eigenvalues will be used to scale the eigenvectors
    if(is_ref)
        not_enough_points = 0;
        empty_seeds = 0;

        if(size(p_projected,1)<3)
            empty_seeds = 1;
            if(size(p_projected,1)==2)
                d = p_projected - SuperFiber.mean(:,axial(i))';
                n = sqrt(diag(d*d')); %this gives the norm of each row
                [~,id_max] = max(n);%get the index of the maximal norm
                V1 = d(id_max,:)';
                tmp = null(V1(:).')';%get vectors perpedicular to V1
                [~,ind] = min(abs(SF_direction(i,:)*tmp'));
                V2 = tmp(ind,:)';
            elseif(size(p_projected,1)==1)
                not_enough_points = 1;
            end
        else
            [coeff,~,~] = pca(p_projected);
            V1 = coeff(:,1);
            V2 = coeff(:,2);
       end
%         [v,d] = eig(cov(p)); %v's columns are the eigenvectors, d is diagonal matrix of the eigenvalues
%         d = 2*sqrt(diag(d));
%         %one of the directions will be in (or very close to) the
%         %SF_direction. Don't use it, but use the other two
%         [val,ind] = max(abs(SF_direction(i,:)*v));
%         v(:,ind) = [];
%         d(ind) = [];
%         %now only two directions are left = these are kindof "perpendicular" to
%         %the superfiber
%         V1 = v(:,2)./norm(v(:,2));
%         V2 = v(:,1)./norm(v(:,1));
        if(not_enough_points && exist('basis','var'))
           basis(i,:) = basis(i-1,:) 
        else
            basis(i,:) = [V1',V2'];
        end
        basis12 = [basis(i,1:3); basis(i,4:6)];
        if(~empty_seeds)
            [seeds] = FindSeeds(basis12,p_projected,mean(p_projected)');%SuperFiber.mean(:,axial(i)));
        end
        ReferenceData.basis = basis;
        ReferenceData.SF_direction = SF_direction;
        ReferenceData.SuperFiber = SuperFiber.mean;
    else
        rotmat = vrrotvec2mat(vrrotvec(SF_direction(i,:),RefData.SF_direction(i,:)));

        basis1 = RefData.basis(i,1:3)*rotmat;
        basis2 = RefData.basis(i,4:6)*rotmat;
        %  basis{subj}(i,1:3) = basis1;  basis{subj}(i,4:6) = basis2;
        basis12 = [basis1;basis2];
        %[seeds] = FindSeeds(basis12,p,SuperFiber.mean(:,axial(i)));
        if (~isempty(p))
            [seeds] = FindSeeds(basis12,p_projected,mean(p_projected)');
        end

    end
   % new_seeds = (ones(length(seeds),1)*SuperFiber.mean(:,axial(i))') + seeds;
    if(~isempty(p) && ~isempty(seeds))
        new_seeds = (ones(length(seeds),1)*mean(p_projected)) + seeds;
    else
        new_seeds = ones(num_seeds_per_slab,1)*SuperFiber.mean(:,round((ub+lb)/2))';
    end

    c = [c;new_seeds];
    seed_points = c; 
end
%IDX = AssignPointsToClusters(seed_points,resampled_fibers);
%%
%
% %
% %load 'C:\MyWork\Stanford\data\data_all.mat';
% load 'C:\MyWork\Stanford\data\data_all_warped_cake_clustering3.mat'
% common_clusters = intersect(intersect(intersect(intersect(IDX{1},IDX{2}),IDX{3}),IDX{4}),IDX{5});
%
% common_clusters = intersect(intersect(intersect(intersect(intersect(IDX{1},IDX{2}),IDX{3}),IDX{4}),IDX{5}),IDX{6});
% % for subj = 1:length(IDX)
% % for i=1:length(common_clusters)
% %     IDX_common{subj} = IDX{subj}(IDX{subj} == common_clusters);
% % end
% % end
% % z=unique(randi(length(common_clusters),1,length(common_clusters)));
% cluster_indices = common_clusters;% intersect(14:3:275,common_clusters);%[1:11];%common_clusters;%(z);
% for subj = 1:length(IDX)
%     % for i=1:length(C{1})
%     %     [~,IDX{subj}(i)] = max(C{subj}(i,:));
%     % end
%     colors = distinguishable_colors(length(cluster_indices)+1);
%     colors(1,:) = [0.5 0.5 0.5];
%     fascicle_colormap = resampled_fibers{subj};
%     cnt = 1;
%     for n_fiber = 1: length(resampled_fibers{subj})
%         for point_in_fiber = 1:length(resampled_fibers{subj}{n_fiber})
%             %fascicle_colormap{n_fiber}(:,point_in_fiber) = [0.5 myVals(cnt) 0]';
%             %fascicle_colormap{n_fiber}(:,point_in_fiber) = [myVals(cnt) 0.5  0]';
%             if(ismember(IDX{subj}(cnt),cluster_indices))
%                 [~, color_index] = find(cluster_indices==IDX{subj}(cnt));
%             else
%                 color_index = 1;
%             end
%             fascicle_colormap{n_fiber}(:,point_in_fiber) = colors(color_index,:)';
%             cnt = cnt+1;
%         end
%     end
%     p.name          = sprintf('Other');
%     p.h             = figure('name',p.name, 'color','b');
%     p.colorType     = 'map';
%     p.faceAlpha     =1;
%     p.nSurfaceFaces = 19;
%     p.colorMap      = hsv(255);
%     p.radius        = .3;
%     p.nNodesMin     = 3;
%     p.color         = fascicle_colormap;
%     [fh1, lh1] = mbaDisplayFascicles(resampled_fibers{subj},p);
%     %view(0,25)
%     delete(lh1)
%     lh = light('Position',[0 0 1]);
%     % view([-90 30]);
%     view([0 100])
%     fh = p.h;
% %     str = ['C:\MyWork\Stanford\Documents\fascicle_fmaps\WarpedCakeClustering\','subject',num2str(subj),'.svg'];
% %     plot2svg(str,1)
% end
%
% %     fascicle_colormap = resampled_fibers{1};
% % cnt = 1;
% % for n_fiber = 1: length(resampled_fibers{1})
% %     for point_in_fiber = 1:length(resampled_fibers{1}{n_fiber})
% %         %fascicle_colormap{n_fiber}(:,point_in_fiber) = [0.5 myVals(cnt) 0]';
% %         %fascicle_colormap{n_fiber}(:,point_in_fiber) = [myVals(cnt) 0.5  0]';
% %         if(ismember(IDX(cnt),z))
% %             [~, color_index] = find(z==IDX(cnt));
% %         else
% %            color_index = 1;
% %         end
% %         fascicle_colormap{n_fiber}(:,point_in_fiber) = colors(color_index,:)';
% %         cnt = cnt+1;
% %     end
% % end
% %
% %
% % p.name          = sprintf('Other');
% %     p.h             = figure('name',p.name, 'color','b');
% %     p.colorType     = 'map';
% %     p.faceAlpha     =1;
% %     p.nSurfaceFaces = 19;
% %     p.colorMap      = hsv(255);
% %     p.radius        = .3;
% %     p.nNodesMin     = 3;
% %     p.color         = fascicle_colormap;
% %     [fh1, lh1] = mbaDisplayFascicles(resampled_fibers{1},p);
% %     %view(0,25)
% %     delete(lh1)
% %     lh = light('Position',[0 0 1]);
% %     % view([-90 30]);
% % view([-90 30])
% %     fh = p.h;
% %     %plot2svg('C:\MyWork\Stanford\Documents\fascicle_fmaps\FAvals\S1.svg',1)
% %
% %     saveas(p.h,['C:\MyWork\Stanford\Documents\fascicle_fmaps\FAvals\s12.png'])
% %
% figure
% i=2;
% starts = [SuperFiber{subj}.mean(:,axial(i))';SuperFiber{subj}.mean(:,axial(i))';SuperFiber{subj}.mean(:,axial(i))']';%zeros(3,3);
% ends = [SuperFiber{subj}.mean(:,axial(i)+8) - SuperFiber{subj}.mean(:,axial(i))]';
% ends = [opposite(1,:);opposite(2,:);SF_direction'];
% hold on;
% %ends = [V;V1; V2];
% quiver3(starts(1,1),starts(2,1),starts(3,1), SF_spread(1)*SF_direction(1), SF_spread(1)*SF_direction(2),SF_spread(1)*SF_direction(3),'r');hold on;
% quiver3(starts(1,2),starts(2,2),starts(3,2), SF_spread(3)*basis(1,1),SF_spread(3)*basis(1,2), SF_spread(3)*basis(1,3),  'color','g');hold on
% quiver3(starts(1,3),starts(2,3),starts(3,3),  SF_spread(1)*basis(2,1),SF_spread(1)*basis(2,2), SF_spread(1)*basis(2,3),'color','b');hold on;
% quiver3(starts(1,3),starts(2,3),starts(3,3), between(1,1), between(1,2),between(1,3),'color','k');hold on;
% quiver3(starts(1,3),starts(2,3),starts(3,3), opposite_between(1,1), opposite_between(1,2),opposite_between(1,3),'color','k');hold on;
% quiver3(starts(1,3),starts(2,3),starts(3,3), opposite(1,1), opposite(1,2),opposite(1,3),'color','g');hold on
% quiver3(starts(1,3),starts(2,3),starts(3,3), opposite(2,1), opposite(2,2),opposite(2,3),'color','b');hold on
% quiver3(starts(1,3),starts(2,3),starts(3,3), between_other(1,1), between_other(1,2),between_other(1,3),'color','k');hold on;
% quiver3(starts(1,3),starts(2,3),starts(3,3), opposite_between_other(1,1), opposite_between_other(1,2),opposite_between_other(1,3),'color','k');hold on;
% %plotcube(SF_spread,[ 0 0 0], 0.2,[0 1 1]);plotcube([-1 -1 1].*SF_spread,[ 0 0 0], 0.2,[0 1 1]);
% plotcube(2*SF_spread,[-1 -1 -1].*SF_spread, 0.2,[0 1 1]);
% figure
% hold on
% quiver3(starts(1,1), starts(2,1), starts(3,1), ends(1,1), ends(1,2), ends(1,3),'color','c');hold on;
% quiver3(starts(1,2), starts(2,2), starts(3,2), ends(2,1), ends(2,2), ends(2,3),'color','r');hold on;
% quiver3(starts(1,3), starts(2,3), starts(3,3), ends(3,1), ends(3,2), ends(3,3),'color','k');
%
% % axis equal
% figure
% k=2
% ind = [3,60:83];
% ind = 14:3:275;%[14,17,20,23,26,29,32,35]
% %ind = [1,12:35,3,60:83,9,204:227,6,132:155];
% plot3(SuperFiber{k}.mean(1,:),SuperFiber{k}.mean(2,:),SuperFiber{k}.mean(3,:),'r','LineWidth',5);hold on;
% plot3(seed_points{k}(ind,1),seed_points{k}(ind,2),seed_points{k}(ind,3),'b.');hold on
% plot3(seed_points{k}(ind,1),seed_points{k}(ind,2),seed_points{k}(ind,3),'bo');hold on
% for ang=0:90
%     view([0 ang])
%     pause(0.03)
% end
% for ang=89:-1:0
%     view([0 ang])
%     pause(0.03)
% end
% figure
% ind = 1:275;
% plot3(SuperFiber{k}.mean(1,:),SuperFiber{k}.mean(2,:),SuperFiber{k}.mean(3,:),'r','LineWidth',5);hold on;
% plot3(seed_points{k}(ind,1),seed_points{k}(ind,2),seed_points{k}(ind,3),'b.');hold on
% plot3(seed_points{k}(ind,1),seed_points{k}(ind,2),seed_points{k}(ind,3),'bo');hold on
% for ang=0:90
%     view([0 ang])
%     pause(0.03)
% end
% for ang=89:-1:0
%     view([0 ang])
%     pause(0.03)
% end