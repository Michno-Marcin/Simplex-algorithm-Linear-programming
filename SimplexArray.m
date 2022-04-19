function Simplex = SimplexArray(CanonicalTab,AreFreeX)
if AreFreeX
    Simplex2 = [[CanonicalTab(1,end-2);CanonicalTab(2:end,end)],CanonicalTab(:,2:end-3)];
else
    Simplex2 = [[0;CanonicalTab(2:end,end)],CanonicalTab(:,2:end-2)];
end
Simplex2=str2double(Simplex2);
fprintf('\t The resulting simplex array from the canonical form is: \n');
disp(Simplex2);
Simplex=Simplex2;
end

