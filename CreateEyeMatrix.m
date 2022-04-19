function Simplex = CreateEyeMatrix(Simplex)
if IsContradiction(Simplex)==0
    DRt=0;
    while IsEyeMatrix(Simplex) == false
        [col,row] = find(Simplex(2+DRt:end,2:end)'); 
        ElementValue = Simplex(row(1)+1+DRt,col(1)+1);
        Simplex(row(1)+1+DRt,:) = Simplex(row(1)+1+DRt,:)/ElementValue;
        for i=1:size(Simplex,1)
            if i == row(1)+1+DRt
                continue;
            else 
                Simplex(i,:) = Simplex(i,:) - Simplex(i,col(1)+1) * Simplex(row(1)+1+DRt,:); 
            end
        end
        DRt=DRt+1;
        fprintf('\t The array after rotation no.%d looks as follows: \n',DRt);
        disp(Simplex);
        if IsContradiction(Simplex)~=0 
            if IsContradiction(Simplex)==1
                fprintf('\t The resulting array is a contradictory form of I kind ! \n\n');
            else
                fprintf('\t The resulting array is a contradictory form of II kind ! \n\n');
            end
            break;
        end
    end
else
    if IsContradiction(Simplex)==1
        fprintf('\t This array is in contradictory form of I kind ! There is no need to create a unitary array. \n\n');
    else
        fprintf('\t This array is in contradictory form of II kind ! There is no need to create a unitary array. \n\n');
    end
end
DRm=0; % Number of rows deleted
for i=size(Simplex,1):-1:1 
        if Simplex(i,:)==0
            Simplex = [Simplex(1:i-1,:);Simplex(i+1:end,:)]; DRm=DRm+1; 
            fprintf('\t The array after removing the redundant row no.%d looks as follows: \n',DRm);
            disp(Simplex);
        end
end
end

