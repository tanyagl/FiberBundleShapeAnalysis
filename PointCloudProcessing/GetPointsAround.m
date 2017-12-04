function [p_within,p_projected] = GetPointsAround(axial,i,resampled_fibers,SF_mean,SF_dir, lb,ub)
p = [];
p_within = [];
sample_around = floor(length(SF_mean)/length(axial));
%we want to first define a 'slab' by finding two planes perpendicular to
%the SF direction at a certain point. A plane is defined by a point and a
%normal, so the normals of the two sides of the slab have to face each
%other.
%Then, for every given point we need to check whether it is between the two
%planes. We can do this as follows:
%rotate the points and the normals such that the planes are parallel to
%the Z-axis, so we will only need to check whether a point is <z2 and >z1.
rotmat = vrrotvec2mat(vrrotvec(SF_dir,[0 0 1]));
lb_z = rotmat*SF_mean(:,lb);lb_z = lb_z(3);
ub_z = rotmat*SF_mean(:,ub);ub_z = ub_z(3);

ub = axial(i)+sample_around;
if ub>length(SF_mean)
    ub = length(SF_mean);
end
lb =  axial(i)-sample_around;
if lb<1
    lb = 1;
end
for fiber=1:length(resampled_fibers)
    p = [p;resampled_fibers{fiber}(:,lb:ub)'];%get all points sampled around the location in original coord system
end
p_rotated = (rotmat*p')'; %rotate the points such that their "z" axis (SF_dir) is aligned with [0 0 1], it makes choosing points that belong to a planar slab easier
p_within = p((p_rotated(:,3)<ub_z & p_rotated(:,3)>=lb_z),:); % this is in original coord system
basisPlane = null(SF_dir);%orthonormal basis for the null space of SF_dir', two vectors on the plane
pointInplane = SF_mean(:,axial(i));%point on the plane
basisCoefficients= bsxfun(@minus,p_within,pointInplane')*basisPlane ; %2d projection of the points onto the basis spread by basisPlane
p_projected=bsxfun(@plus,basisCoefficients*basisPlane.', pointInplane');
%Idea taken from here: https://www.mathworks.com/matlabcentral/newsreader/view_thread/293244


