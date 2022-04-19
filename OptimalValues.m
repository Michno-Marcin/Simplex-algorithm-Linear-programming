function XValues = OptimalValues(Simplex)
I = [zeros(1,size(Simplex,1)-1);eye(size(Simplex,1)-1)]; 
UnitMatrixNr=int8.empty;
for i=1:size(Simplex,1)-1 % Finding the column numbers corresponding to the unit matrix
    UnitMatrixNr(i) = find(ismember(Simplex',I(:,i)','rows'));
end
XValues = cell(1, size(Simplex,2)-1);
for i=1:length(UnitMatrixNr)
    XValues{UnitMatrixNr(i)-1} = Simplex(1+i,1); % Value assignment to items in unit columns
end
for i=1:length(XValues)
    if ~isempty(XValues{i}) == 0 
        XValues{i}=0; % Value assignment: zero for all other elements
    end
end
end