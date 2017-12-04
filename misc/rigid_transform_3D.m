function [R,t] = rigid_transform_3D(A, B)

    assert(size(A,1) == size(B,1) && size(A,2) == size(B,2))

    centroid_A = mean(A);
    centroid_B = mean(B);

    N = size(A,1);

    H = (A - repmat(centroid_A, N, 1))' * (B - repmat(centroid_B, N, 1));

    [U,S,V] = svd(H);

    R = V*U';

    if det(R) < 0
        sprintf('Reflection detected\n');
        V(:,3)  = -1*V(:,3);
        R = V*U';
    end

    t = -R*centroid_A' + centroid_B';
end

