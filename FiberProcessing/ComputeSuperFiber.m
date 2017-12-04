function [SuperFiber,endpoints1, endpoints2] = ComputeSuperFiber(fibers,numNodes)
SuperFiber.mean = zeros(3,numNodes);
SuperFiber.spread = zeros(3,numNodes);
%allPoints=[];
if(length(fibers)<2)
    stop_here = 1;
end
% for i=1:length(fibers)
%     fibers_resampled{i} = dtiFiberResample(cell2mat(fibers(i)), numNodes);
%     allPoints=[allPoints;fibers_resampled{i}'];
% end
%new_fibers = fibers_resampled';
%new_fibers = fibers;
empty_fibers_cnt = 0;
for node=1:numNodes
    for i=1:length(fibers)
        if(isempty(fibers{i}))
            empty_fibers_cnt = empty_fibers_cnt+1;
            sprintf('%d',i)
        else
            cur_node(i,:) = fibers{i}(:,node);
        end
    end
    C = calcCentroid(cur_node);
    D = cur_node-C;
    N = sqrt(sum(D.^2,2));
    P = prctile(N,70);
    subset_cur_node = cur_node(N<P,:);
    SuperFiber.mean(:,node) = mean(subset_cur_node)';
    SuperFiber.spread(:,node) = std(subset_cur_node)';
    if(node == 1)
        endpoints1 = cur_node;
    elseif(node == numNodes)
        endpoints2 = cur_node;
    end
    clear cur_node;
end
% normalized_coords = allPoints - repmat(calcCentroid(allPoints),[length(allPoints) 1]);
% [ PCAed_coords,~,~ ] = nipals(normalized_coords,3);
% cnt = 1;
% new_fibers = fibers;
% for i=1:s, %go back to fibers structure
%    new_fibers{i}= PCAed_coords(cnt:cnt+length(cell2mat(fibers(i)))-1,:)';
%    cnt = cnt+length(cell2mat(fibers(i)));
% end
% 
% 
% 
% 
% for i=1:length(fibers)
% %    fibers_resampled{i,1} = dtiFiberResample(cell2mat(fibers(i)), numNodes);
%     fibers_resampled{i,1} = dtiFiberResample(cell2mat(new_fibers(i)), numNodes);
% 
% end
% for node=1:numNodes
%     for i=1:length(fibers_resampled)
%         cur_node(i,:) = fibers_resampled{i}(:,node);
%     end
%     SuperFiber.mean(:,node) = mean(cur_node)';
%     SuperFiber.spread(:,node) = std(cur_node)';
% end
