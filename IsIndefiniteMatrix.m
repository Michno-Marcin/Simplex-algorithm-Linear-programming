function IsInd = IsIndefiniteMatrix(Simplex)
IsInd=false;
if ~IsContradiction(Simplex)==1
    if ~IsContradiction(Simplex)==2
        for i=2:size(Simplex,2)
            if Simplex(1,i)<0 % When there is any negative free expression in the objective function
                if Simplex(:,i)<=0 % When each element of a column is less than or equal to zero
                    IsInd = true; % Unrestricted form
                    break;
                end
            end
        end
    end
end
end

