function CanonicalArr = CanonicalArray(NoOfX,Z,a,b,TypesBon,TypesX)
% Calling a ready-made function for specified parameters
% SimplexArray(Z,a,b,TypesBon,TypesX)
% Z - vector of objective function weights (coefficient values)
% a - bounding matrix (values of the coefficients for the bounding conditions)
% AreFreeX - auxiliary variable indicating whether there is a free expression 
%   in the coefficients of the objective function
AreFreeX = (NoOfX~=length(Z));
% b - vector of free expressions
% Assume that: 1 means ">=", 0 means "=", -1 means "<=":
% TypesBon - inequality types for bounding conditions
% TypesX - types of boundings for the decision variables of the problem

if AreFreeX
    LeftMatrix=[-Z;a,zeros(NoOfX,1)]; % Converted to a minimization problem 
    % (reversing the sign of the objective function value)
else
    LeftMatrix=[-Z;a];
end
% Removing inequalities - creating free variables with the help of a "adjugate" matrix
% Numbers: -1,1 denote here the added free variables for the inequality.
% They do not affect the result of the objective function (values of its last coefficients equal zero). 
% The first row of the "adjugate" matrix indicates the coefficients of the objective function.
a2 = zeros(size(a,1)+1,size(a,1)); % First we create a null matrix
for i=1:length(TypesBon) % Assigning row/column values to added free variables
    if(TypesBon(i)==1 || TypesBon(i)==-1)
        a2(i+1,i) = - TypesBon(i);
    end
end

% Searching the columns of the resulting matrix to find those without added free variables and remove them
i=size(a2,2);    
while i>=1
    if (a2(i+1,i)==0)
        a2(:,i) = [];
        i=i-1;
    end
    i=i-1;
end

% Creating text fields for matrix with structure and simplex array information
NewColForX = cell(1, size(LeftMatrix,1)-1);
for i=1:size(LeftMatrix,1)-1
     NewColForX{i}=convertCharsToStrings(['RestrictiveCon.',num2str(i),':']);
end 
NewColForX=["ObjectiveFunction:";NewColForX'];
Comp=string.empty;
for i=1:size(LeftMatrix,1)-1
     Comp(i)="=";
end 
Comp=Comp';
Comp(i+1)="-> Min";

% Displaying the simplex array after the first creation step
% In the figure, all columns except the first and *three 
% (two when there is no free expression in the objective function) last are assigned to the corresponding decision variables.
% The last column is the free expressions, and the third from the end are the last coefficients of the variables
if AreFreeX
    CanonicalArr=[NewColForX,[LeftMatrix(:,1:end-1),a2,LeftMatrix(:,end),[Comp,['FreeExpressions'; string(b)]]]];
else
    CanonicalArr=[NewColForX,[LeftMatrix,a2,[Comp,['FreeExpressions'; string(b)]]]];
end
fprintf('\t The resulting matrix after the first substitution step into canonical form: \n');
disp(CanonicalArr); % The resulting figure after the first substitution step

%% Creating the canonical form for the maximization problem - second stage of the conversion
% Including: obtaining non-negative decision variables

% Convertsion of negative decision variables to positive ones
for i=1:length(TypesX) 
    if(TypesX(i)==-1)
        LeftMatrix(:,i) = -LeftMatrix(:,i);
    end
end

% Conversion of any decision variables values to positive 
for i=length(TypesX):-1:1 
    if(TypesX(i)==0)
        LeftMatrix=[LeftMatrix(:,1:i),-LeftMatrix(:,i),LeftMatrix(:,i+1:size(LeftMatrix,2))];
    end
end

% Displaying the created canonical array
if AreFreeX
    CanonicalArr=[NewColForX,[LeftMatrix(:,1:end-1),a2,LeftMatrix(:,end),[Comp,['FreeExpressions'; string(b)]]]];
else
    CanonicalArr=[NewColForX,[LeftMatrix,a2,[Comp,['FreeExpressions'; string(b)]]]];
end
fprintf('\t The resulting matrix after the second substitution step into canonical form is: \n');
disp(CanonicalArr); % The resulting figure after the second substitution step
end