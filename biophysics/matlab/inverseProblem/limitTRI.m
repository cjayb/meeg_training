function [newgrid, newtri] = limitTRI(grid, tri)

newgrid = grid;
newtri = tri;
cont = 1; ii = 1;
ind = flipud(sort(find(newgrid(:,3) > -0.03)));
while cont,
    
    curtris = size(newtri(:,1),1);
    for j = curtris:-1:1,
        if (newtri(j,1) == ind(ii) || newtri(j,2) == ind(ii) || newtri(j,3) == ind(ii)),
            %                 newtri = [newtri(1:(j-1),:); newtri(min(curtris, j+1):end,:)];
            newtri = [newtri(1:(j-1),:); newtri(j+1:end,:)];
            curtris = curtris - 1;
            
            %else
        end
    end
    temp = newtri(:);
    temp(temp > ind(ii)) = temp(temp > ind(ii)) - 1;
    try
        newtri = reshape(temp, curtris, 3);
    catch
        newtri = newtri;
    end
    
    all = [1:size(newgrid, 1)];
    all = all(all~=ind(ii));
    newgrid = newgrid(all,:);
    
    ind = sort(find(newgrid(:,3) > -0.03), 'descend');
    if isempty(ind),
        cont = 0;
    end
    
    %         if ~mod(length(ind), 10),
    %             fprintf(1, '%d extra points left...\n', length(ind));
    %         end
end
newgrid(:,3) = -1*newgrid(:,3);

