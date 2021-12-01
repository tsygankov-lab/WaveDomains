function output = merge_domains(WS, theta, merge_theta_th, merge_len_th)

% Merge domains based on the length of the border and median value of theta
% in each domain

%INPUT
%WS - 2d array with domains labels
%theta - 2d array with values of theta
%merge_theta_th - [DOUBLE] max difference on mean theta to merge domains
%merge_len_th - min fraction of the border in comparison to the characteristic size of the domain to merge domains

%OUTPUT
%output - 2d array with merged domains

d_ids = unique(WS(:)); %domains ids
d_ids = d_ids(2:end); %exclude zero

%characteristic lengths of domains
ch_lens = zeros(1, length(d_ids));

% cell array with indexes of intersecting borders of domains
b_inter = cell(length(d_ids), length(d_ids));

% matrix of merging domains
merge_m = zeros(length(d_ids), length(d_ids));

%before running the algorithm, compute necessary data
for i_id = 1:(length(d_ids))
    i = d_ids(i_id);
    
    %compute characteristic lenght of the domain
    ch_lens(i_id) = sqrt(2/pi*sum(WS(:) == i));
end
for i_id = 1:length(d_ids)
    i = d_ids(i_id);
    %cell outline of the domain
    CO1 = cell_outline_v2(WS == i);
    
    %find ids of the domains that are in contact with domain i
    d_ids2 = unique(WS(cell_outline_v2((WS == i) + CO1) == 1));
    d_ids2 = d_ids2(2:end); %exclude zero
    
    for j_id = 1:length(d_ids2)
        j = d_ids2(j_id);
        
        %compute indexies by which two domains interact 
        %intersect outlines of two domains and exclude pixels, which
        %intersect with outlines of other domains (it can connect domains which should be separated)
        inter = CO1 & (cell_outline_v2(WS == j)) & ~(cell_outline_v2(WS > 0 & WS ~= i & WS ~= j));
        b_inter{i,j} = find(inter);
        b_inter{j,i} = b_inter{i,j};
        
        %compute theta distance
        theta_dist = abs(circ_dist(circ_mean(theta(WS == i)), circ_mean(theta(WS == j))));
        
        %initiate elements of merging matrix if neccesary conditions are
        %met
        if ~isempty(b_inter{i,j}) && ...
                theta_dist <= merge_theta_th && ...
                length(b_inter{i,j}) >= merge_len_th*min(ch_lens(i), ch_lens(j))
            merge_m(i,j) = 1;
            merge_m(j,i) = 1;
        end
    end
end

%add connection of domains for neighbors (if one domain is connected to two 
%domains, these domais also should be connected)
%in the other words, idetify cliques and make them connected

merge_m1 = merge_m;
merge_m2 = merge_m;
check = true; % check if some non-added nodes still exist

%update matrix until covergence
while check
    check = false;
    for i_id = 1:(length(d_ids)) 
        for j_id = 1:length(d_ids)
            if merge_m1(i_id,j_id) == 1 && sum(merge_m1(i_id,:)) >= 1 && sum(merge_m1(:,j_id)) >= 1
                row = merge_m1(i_id,:); %row(j) = 0;
                col = merge_m1(:,j_id); %col(i) = 0;
                merge_m2(row == 1, j_id) = 1;
                merge_m2(i_id, col == 1) = 1;
            end
        end
    end
    if ~isequal(merge_m1, merge_m2)
        check = true;
        merge_m1 = merge_m2;
    end
end

mask = zeros(size(WS));

%merge domains
for i_id = 1:(length(d_ids)-1)
    i = d_ids(i_id);
    %for j_id = (i_id+1):length(d_ids)
    for j_id = 1:length(d_ids)
        j = d_ids(j_id);
        
        mask(WS == i) = 1;
        mask(WS == j) = 1;
        
        if merge_m2(i_id, j_id) == 1 && length(b_inter{j_id,i_id}) >= 1
            mask(b_inter{i_id,j_id}) = 1;
        end
    end
end

%filter out separated pixels from mask
%such pixels appear after merging multiple interfaces of domains
mask = ~bwareaopen(~mask,2);

output = bwlabel(mask);

end