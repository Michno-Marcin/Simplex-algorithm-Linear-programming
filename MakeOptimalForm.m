function Simplex = MakeOptimalForm(Simplex)
fprintf('\t We will determine the optimal array through simplex rotations. \n');
DRt=0;
while max(~(Simplex(1,2:end)>=0)) % As long as there is a negative coefficient in the objective function.
    ColNr = find(Simplex(1,2:end)==min(Simplex(1,2:end)),1)+1; % Rotation column selection
    CCQuotient = Simplex(2:end,1) ./ Simplex(2:end,ColNr);
    RowNr = find(CCQuotient == min(CCQuotient(CCQuotient>0)))+1; % Rotation row selection
    MinEl = Simplex(RowNr,ColNr);
    Simplex(RowNr,:) = Simplex(RowNr,:)/MinEl; % The element found is replaced by 1.
    for j=1:size(Simplex,1) % In each of the remaining rows, zero the elements on the rotation column.
        if j == RowNr
            continue;
        else
            Simplex(j,:) = Simplex(j,:) - Simplex(j,ColNr) * Simplex(RowNr,:); 
        end
    end
    DRt=DRt+1;
    fprintf('\t The array after rotating array no.%d looks as follows: \n',DRt);
    disp(Simplex);
end
end

