function Simplex = MakeBaseForm(Simplex)
Simplex2=Simplex;
if IsContradiction(Simplex)==0
    fprintf('\t We will now create an array in base form. \n');
    IsBaseForm = false; DRt = 0;
    while ~IsBaseForm
        if (IsContradiction(Simplex)==0)
            IsBaseForm = true;
            for i=2:size(Simplex,1) 
                if Simplex(i,1)<0
                    IsBaseForm = false;
                    [Value,MinCol] = min(Simplex(i,2:end)); MinCol=MinCol+1;
                    Quotients = [Simplex(1:i-1,1);nan;Simplex(i+1:end,1)] ./ [Simplex(1:i-1,MinCol);nan;Simplex(i+1:end,MinCol)];
                    Quotients(Quotients<=0) = nan; 
                    MinQuotient = min(Quotients);
                    MinRow = find(Quotients==MinQuotient,1);
                    Min = Simplex(MinRow, MinCol);
                    if (isnan(Min) || MinQuotient == Inf)
                        Min = Value; 
                        MinRow = i;
                    end
                    Simplex(i,:) = Simplex(i,:)/Min;
                    for j=1:size(Simplex,1) 
                        if j == MinRow
                            continue;
                        else
                            Simplex(j,:) = Simplex(j,:) - Simplex(j,MinCol) * Simplex(MinRow,:); 
                        end
                    end
                    DRt=DRt+1;
                    fprintf('\t The array after rotating array no.%d looks as follows: \n',DRt);
                    disp(Simplex);
                end
            end
        else
            if IsContradiction(Simplex)==1
                fprintf('\t The resulting array is in contradictory form of I kind !\n\n'); break;
            else
                fprintf('\t The resulting array is in contradictory form of II kind !\n\n'); break;
            end
        end
    end
end
if IsContradiction(Simplex)==0
    if(Simplex2==Simplex)
        fprintf('\t But the array above is a simplex array in base form with columns of a unitary array, so no actions needed. \n');
    else
        fprintf('\t Above, a simplex array was created in base form with columns of a unitary array. \n');
    end
end
end


