function [B, Bd] = assemble_B(s, deg, p)

% serial B and Bd computation
B  = spalloc(sum(deg-1), length(p), sum( 2*(deg-1) ));
Bd = spalloc(sum(deg-1), length(p), sum( deg.*(deg-1) ));
idx = 1;
for i = 1 : length(s)
    ind = find( p == s(i) );
    count = 1;
    for j = 1 : deg(i)-1
        B(idx,ind(j)) = 1;
        B(idx,ind(j+1)) = -1;
        for k = 1 : j
            Bd(idx,ind(k)) = 1-count/deg(i);
        end
        for k = j+1 : length(ind)
            Bd(idx,ind(k)) = -count/deg(i);
        end
        count = count+1;
        idx = idx+1;
    end
end

% % distributed B and Bd computation
% BBc =cell(ndoms,1);
% BBdc=cell(ndoms,1);
% for curr_dom=1:ndoms
%     BBc{curr_dom} =spalloc(sum(deg-1),length(sep{curr_dom}),sum(deg-1));
%     BBdc{curr_dom}=spalloc(sum(deg-1),length(sep{curr_dom}),sum(deg-1));
%     idx=1;
%     for i=1:length(s)
% %         ind = find( All{curr_dom} == s(i) );
%         ind = find( sep{curr_dom} == s(i) );
%         if ( ~isempty(ind) )
%             prev_dom = curr_dom - 1;
%             pprev = p( 1 : blk(prev_dom+1)-1 );
%             ind_prev = find( pprev == s(i) );
%             len = length(ind_prev);
%             if (len == 0)
%                 BBc{curr_dom}(idx,ind) = 1;
%             elseif (len == deg(i)-1)
%                 BBc{curr_dom}(idx+deg(i)-2,ind) = -1;
%             else
%                 BBc{curr_dom}(idx+len-1,ind) = -1;
%                 BBc{curr_dom}(idx+len,ind) = 1;
%             end
%             count = 1;
%             for j=idx:(idx+len-1)
%                 BBdc{curr_dom}(j,ind) = -count/deg(i);
%                 count = count + 1;
%             end
%             for j=(idx+len):(idx+deg(i)-2)
%                 BBdc{curr_dom}(j,ind) = 1-count/deg(i);
%                 count = count + 1;
%             end
%         end
%         idx = idx + deg(i) - 1;
%     end
% end
% Bc=BBc;
% Bdc=BBdc;

end
