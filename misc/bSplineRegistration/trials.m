%# these don't all have to be the same
x = -8:2:8; y = -8:2:8; z = -8:2:8;

[X1 Y1 Z1] = meshgrid(x([1 end]),y,z);
X1 = permute(X1,[2 1 3]); Y1 = permute(Y1,[2 1 3]); Z1 = permute(Z1,[2 1 3]);
X1(end+1,:,:) = NaN; Y1(end+1,:,:) = NaN; Z1(end+1,:,:) = NaN;
[X2 Y2 Z2] = meshgrid(x,y([1 end]),z);
X2(end+1,:,:) = NaN; Y2(end+1,:,:) = NaN; Z2(end+1,:,:) = NaN;
[X3 Y3 Z3] = meshgrid(x,y,z([1 end]));
X3 = permute(X3,[3 1 2]); Y3 = permute(Y3,[3 1 2]); Z3 = permute(Z3,[3 1 2]);
X3(end+1,:,:) = NaN; Y3(end+1,:,:) = NaN; Z3(end+1,:,:) = NaN;

%#figure('Renderer','opengl')
h = line([X1(:);X2(:);X3(:)], [Y1(:);Y2(:);Y3(:)], [Z1(:);Z2(:);Z3(:)]);
set(h, 'Color',[0.5 0.5 1], 'LineWidth',1, 'LineStyle','-')

%#set(gca, 'Box','on', 'LineWidth',2, 'XTick',x, 'YTick',y, 'ZTick',z, ...
%#  'XLim',[x(1) x(end)], 'YLim',[y(1) y(end)], 'ZLim',[z(1) z(end)])
%#xlabel x, ylabel y, zlabel z
axis off
view(3), axis vis3d
camproj perspective, rotate3d on

figure;
i = 0:0.2:2;
[X Y] = meshgrid(i,i);                         
x = [X(:) X(:)]';                                
y = [Y(:) Y(:)]';
z = [repmat(i(1),1,length(x)); repmat(i(end),1,length(x))];
col = 'b';
hold on;
plot3(x,y,z,col);                                         
plot3(y,z,x,col);
plot3(z,x,y,col);

figure
for i=1:25:size(O_trans,1)
    for j=1:75:size(O_trans,2)
        for k=1:50:size(O_trans,3)
            quiver3(i, j, k, squeeze(O_trans(i,j,k,1)),squeeze(O_trans(i,j,k,2)),squeeze(O_trans(i,j,k,3)),'Color','b');
            hold on
        end
    end
end


%%
load 'C:\MyWork\Stanford\data\KYEOM\misc\pedunclesWfibers.mat'
younger = LeftICP{1}.fibers{1}';
older = LeftICP{37}.fibers{1}';
 older_resampled = dtiFiberResample(older', size(younger,1))';
 options.Verbose=true;
options.MaxRef=5;
% Fit a bspline grid to transform points Xstatic into Xmoving
[O_trans,Spacing]=point_registration([6 31 19],younger,older_resampled,options);
% Transforming some other point with the b-spline grid
younger_transformed = bspline_trans_points_double(O_trans,Spacing,younger);