function Simplex = OtherResult(Simplex)
Simplex = CreateEyeMatrix(Simplex); % Check for another optimal solution of the given simplex array (we assume no radius of optimality)   
[~, UnitMatrixColumns] = IsEyeMatrix(Simplex);
if length(find(Simplex(1,2:end)==0))>length(UnitMatrixColumns) % Check for non-baseline zeros in the objective function
    NonBaseX=setdiff(find(Simplex(1,2:end)==0)+1,str2double(UnitMatrixColumns));
    % We assume no optimality radius, i.e., non-baseline columns contain at least one value greater than zero
    if (Simplex(:,NonBaseX)<=0)         
        return;
    else
        j=0; Stop = false; Text=true;
        while ~Stop
            if j+1>length(NonBaseX) % If there is no rotation row, there will be a certain radius of optimality that we omit
                Stop = true; Text=false;
            else
                CCQuotient = Simplex(2:end,1) ./ Simplex(2:end,NonBaseX(j+1));
                MinQuo = min(CCQuotient(CCQuotient>0));
                if isempty(MinQuo)
                    j=j+1;
                    continue;
                else
                    RowEl = find(CCQuotient == MinQuo)+1; % Rotation row selection
                    Simplex = Rotate(Simplex,RowEl,NonBaseX(j+1)); % Array rotation
                    j=j+1; Stop = ~isempty(RowEl); % Condition to stop the loop when a rotation row is found
                end
            end
        end
        if Text
            Vec=[1:length(OptimalValues(Simplex)); cell2mat(OptimalValues(Simplex))]; Vec=Vec(:)';
            fprintf('\t The array with another optimal solution looks as follows: \n'); disp(Simplex);
            fprintf('\t The values of the resulting variables are: \n'); fprintf('\t x%d: %d\n',Vec); 
        end
    end
end
end

