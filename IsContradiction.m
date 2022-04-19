function IsContr = IsContradiction(Simplex)
IsContr=0;
for i=2:size(Simplex,1)
    if (Simplex(i,1)~= 0 && ~any(Simplex(i,2:end)~=0))
        IsContr = 1;
    end
end
for i=2:size(Simplex,1)
    if ((Simplex(i,1)> 0 && ~any(Simplex(i,2:end)>0)) || (Simplex(i,1)< 0 && ~any(Simplex(i,2:end)<0)))
        IsContr=2;
    end
end
end

