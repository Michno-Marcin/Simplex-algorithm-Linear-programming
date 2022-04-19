function [IsEyeMat, UnitMatrixColumns, UnitMatrixColumns2] = IsEyeMatrix(Simplex)

if IsContradiction(Simplex)==0
     % Below is a unitary array with as many rows as there are in the simplex array with zeros in the first row
     I = [zeros(1,size(Simplex,1)-1);eye(size(Simplex,1)-1)]; UnitMatrixNo = cell(1, size(Simplex,1)-1);
     for i=1:size(Simplex,1)-1 % Finding the column numbers corresponding to the unit matrix
         UnitMatrixNo{i} = find(ismember(Simplex',I(:,i)','rows'));
     end
     UnitMatrixColumns2=string.empty; i=1;
     while i<=size(Simplex,1)-1
         if isempty(UnitMatrixNo{i})
             IsEyeMat = false;
             break;
         else
             UnitMatrixColumns2(i) =UnitMatrixNo{i};
         end
         if i == size(Simplex,1)-1
             IsEyeMat = true;
         end
         i=i+1;
     end
     UnitMatrixColumns = sort(UnitMatrixColumns2);
 else
     if IsContradiction(Simplex)==1
         fprintf('\t This array is in contradictory form of I kind ! There is no need to look for a unitary array. \n\n');
     else
         fprintf('\t This array is in contradictory form of II kind ! There is no need to look for a unitary array. \n\n');
     end 
end
end
