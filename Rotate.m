function Simplex = Rotate(Simplex,RowEl,ColEl)  
Simplex(RowEl,:) = Simplex(RowEl,:)/Simplex(RowEl,ColEl); % The element found is replaced by 1
for j=1:size(Simplex,1) % In each of the remaining rows, zero the elements on the rotation column
    if j == RowEl
        continue;
    else
        Simplex(j,:) = Simplex(j,:) - Simplex(j,ColEl) * Simplex(RowEl,:);
    end
end
end

