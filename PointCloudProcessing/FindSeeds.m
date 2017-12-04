function [seeds] = FindSeeds(basis12,p,centerPoint)
basis1 = basis12(1,:);
basis2 = basis12(2,:);
between = (basis12(1,:) + basis12(2,:))./2; % ((V1+V2)./2)'
between = between./norm(between);
opposite = -basis12;
opposite = opposite./norm(opposite);
opposite1 = opposite(1,:);
opposite2 = opposite(2,:);
opposite_between = -between;
opposite_between = opposite_between./norm(opposite_between);
between_other = (basis12(1,:) + opposite(2,:))./2;
between_other = between_other./norm(between_other);
opposite_between_other = -between_other;
opposite_between_other = opposite_between_other./norm(opposite_between_other);
%next calculate how far from the superfiber each direction goes
%- meaning the distance from the superfiber to the farthest
%point in each direction along the basis vectors. We're actually not taking
%the farthest point, but rather sorting the points by distance and taking
%the 90th percentile. This is to avoid situations like single outlier fiber
%very distant from the bundle.
max_basis1 = 0;
max_basis2 = 0;
max_between = 0;
max_between_other = 0; 
%To get all the points that belong to a particular slice in the slab we do
%the following:
%(1)calculate the vectors from centerPoint to each of the points in p
%(2)calculate the angles between these vectors and the 8 main directions
%(at 12,1.5,3,4.5,6,7.5,9,10.5 o'clock). This is done by calculating the
%angles for the first 4 directions and taking 180-whatever was the output
%for the opposite directions
%(3)the points that belong to a particular slice in the slab have an angle
%between 0-45 or 135-180. So in the code below, points for which
%angle_basis1 is within these values will belong to the slab determined by
%basis1 direction.
vector_to_point =   bsxfun(@minus,p,centerPoint');
v2p_normalized = normr(vector_to_point);
for i=1:size(vector_to_point,1)
    angle_basis1(i) = atan2d(norm(cross(v2p_normalized(i,:),basis1)), dot(v2p_normalized(i,:),basis1));
    angle_basis2(i) = atan2d(norm(cross(v2p_normalized(i,:),basis2)), dot(v2p_normalized(i,:),basis2));
    angle_between(i) = atan2d(norm(cross(v2p_normalized(i,:),between)), dot(v2p_normalized(i,:),between));
    angle_between_other(i) = atan2d(norm(cross(v2p_normalized(i,:),between_other)), dot(v2p_normalized(i,:),between_other));
    angle_opposite1(i) = atan2d(norm(cross(v2p_normalized(i,:),opposite1)), dot(v2p_normalized(i,:),opposite1));
    angle_opposite2(i) = atan2d(norm(cross(v2p_normalized(i,:),opposite2)), dot(v2p_normalized(i,:),opposite2));
    angle_opposite_between(i) = atan2d(norm(cross(v2p_normalized(i,:),opposite_between)), dot(v2p_normalized(i,:),opposite_between));
    angle_opposite_between_other(i) = atan2d(norm(cross(v2p_normalized(i,:),opposite_between_other)), dot(v2p_normalized(i,:),opposite_between_other));

end
if(~exist('angle_basis1'))
    stop_here = 1;
end
percentile = 90;
% basis1_dist = vector_to_point*basis1';
% basis1_dist_norm = normr(vector_to_point)*normr(basis1)';
% tmp = angle_basis1(basis1_dist_norm>-0.5 & basis1_dist_norm<0.5);
% max_basis1 = prctile(vector_to_point(tmp<45 && tmp>135),percentile);
cnt = 0;
tmp = vector_to_point(angle_basis1<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
cnt = cnt+length(tmp);
max_basis1 = nanmax(0,prctile(tmp_norm,percentile));
tmp = [];
tmp = vector_to_point(angle_basis2<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
max_basis2 = nanmax(0,prctile(tmp_norm,percentile));
cnt = cnt+length(tmp);

tmp = [];
tmp = vector_to_point(angle_between<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
max_between = nanmax(0,prctile(tmp_norm,percentile));
cnt = cnt+length(tmp);
tmp = [];
tmp = vector_to_point(angle_between_other<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
max_between_other = nanmax(0,prctile(tmp_norm,percentile));
cnt = cnt+length(tmp);
tmp = [];
tmp = vector_to_point(angle_opposite1<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
max_opposite1 = nanmax(0,prctile(tmp_norm,percentile));
cnt = cnt+length(tmp);
tmp = [];
tmp = vector_to_point(angle_opposite2<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
max_opposite2 = nanmax(0,prctile(tmp_norm,percentile));
cnt = cnt+length(tmp);
tmp = [];
tmp = vector_to_point(angle_opposite_between<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
max_opposite_between = nanmax(0,prctile(tmp_norm,percentile));
cnt = cnt+length(tmp);
tmp = [];
tmp = vector_to_point(angle_opposite_between_other<22.5,:); tmp_norm = sqrt(sum(tmp.^2,2));
max_opposite_between_other = nanmax(0,prctile(tmp_norm,percentile));
cnt = cnt+length(tmp);
tmp = [];
% col = distinguishable_colors(20);
% figure;image(reshape(col,[1 size(col)]))
% col1 = [col(6:11,:);col(13:18,:);col(20,:)];
% col = col1;
% col(2,:) = [];
% figure;myPlot3(p(angle_basis1<22.5,:),col(1,:),'.',10,'none')
% hold on;myPlot3(p(angle_between<22.5,:),col(2,:),'.',10,'none')
% hold on;myPlot3(p(angle_basis2<22.5,:),col(3,:),'.',10,'none')
% hold on;myPlot3(p(angle_opposite_between<22.5,:),col(4,:),'.',10,'none')
% hold on;myPlot3(p(angle_opposite_between_other<22.5,:),col(5,:),'.',10,'none')
% hold on;myPlot3(p(angle_opposite1<22.5,:),col(6,:),'.',10,'none')
% hold on;myPlot3(p(angle_opposite2<22.5,:),col(7,:),'.',10,'none')
% hold on;myPlot3(p(angle_between_other<22.5,:),col(8,:),'.',10,'none')
% 
% axis equal
% set(gca, 'visible', 'off') ;
% basis1_dist = vector_to_point*basis1';
% max_basis1 = prctile(basis1_dist(basis1_dist>0),percentile);
% max_basis1_opposite = prctile(abs(basis1_dist(basis1_dist<0)),percentile);
% 
% basis2_dist = vector_to_point*basis2';
% max_basis2 = prctile(basis2_dist(basis2_dist>0),percentile);
% max_basis2_opposite = prctile(abs(basis2_dist(basis2_dist<0)),percentile);
% 
% between_dist = vector_to_point*between';
% max_between = prctile(between_dist(between_dist>0),percentile);
% max_between_opposite = prctile(abs(between_dist(between_dist<0)),percentile);
% 
% between_other_dist = vector_to_point*between_other';
% max_between_other = prctile(between_other_dist(between_other_dist>0),percentile);
% max_between_other_opposite = prctile(abs(between_other_dist(between_other_dist<0)),percentile);

% % basis1_dist = sort(abs(vector_to_point*basis1'));%sorting the dot product of the direction vector (basis1) with the vector of distance between the point and the centerPoint
% % basis2_dist = sort(abs(vector_to_point*basis2'));
% % between_dist = sort(abs(vector_to_point*between'));
% % between_other_dist = sort(abs(vector_to_point*between_other'));
% % max_basis1 = prctile(basis1_dist,99);
% % max_basis2 = prctile(basis2_dist,99);
% % max_between = prctile(between_dist,99);
% % max_between_other = prctile(between_other_dist,99);
%now that we have how far from the superfiber each direction
%goes, we can create seed points in appropriate distances
num_radial_cuts = 2;%was 3
basis1_axis = (max_basis1/num_radial_cuts):(max_basis1/num_radial_cuts): max_basis1;
if(isempty(basis1_axis))
    basis1_axis = [0 0];
end
basis2_axis = (max_basis2/num_radial_cuts):(max_basis2/num_radial_cuts): max_basis2;
if(isempty(basis2_axis))
    basis2_axis = [0 0];
end

between_axis = (max_between/num_radial_cuts):(max_between/num_radial_cuts): max_between;
if(isempty(between_axis))
    between_axis = [0 0];
end

between_other_axis = (max_between_other/num_radial_cuts):(max_between_other/num_radial_cuts): max_between_other;
if(isempty(between_other_axis))
    between_other_axis = [0 0];
end

opposite_basis1_axis = (max_opposite1/num_radial_cuts):(max_opposite1/num_radial_cuts): max_opposite1;
if(isempty(opposite_basis1_axis))
    opposite_basis1_axis = [0 0];
end

opposite_basis2_axis = (max_opposite2/num_radial_cuts):(max_opposite2/num_radial_cuts): max_opposite2;
if(isempty(opposite_basis2_axis))
    opposite_basis2_axis = [0 0];
end

opposite_between_axis = (max_opposite_between/num_radial_cuts):(max_opposite_between/num_radial_cuts): max_opposite_between;
if(isempty(opposite_between_axis))
    opposite_between_axis = [0 0];
end

other_opposite_between_axis = (max_opposite_between_other/num_radial_cuts):(max_opposite_between_other/num_radial_cuts): max_opposite_between_other;
if(isempty(other_opposite_between_axis))
    other_opposite_between_axis = [0 0];
end

axis_seeds = [basis1_axis'*basis12(1,:);basis2_axis'*basis12(2,:);opposite_basis1_axis'*opposite(1,:);opposite_basis2_axis'*opposite(2,:)];
between_seeds = [between_axis'*between;between_other_axis'*between_other;opposite_between_axis'*opposite_between;other_opposite_between_axis'*opposite_between_other];
% %now calculate the seed points coords
% axis_seeds = [basis1_axis'*basis12(1,:);basis2_axis'*basis12(2,:);basis1_axis'*opposite(1,:);basis2_axis'*opposite(2,:)];
% between_seeds = [between_axis'*between;between_other_axis'*between_other;opposite_between'*opposite_between;between_other_axis'*opposite_between_other];

seeds = [axis_seeds;between_seeds];
if(length(seeds)<16)
    stop_here = 1;
end