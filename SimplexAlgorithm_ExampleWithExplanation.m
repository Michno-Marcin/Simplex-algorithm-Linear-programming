%% The simplex method of linear programming - implementation in Matlab language
% Including: the canonical form of the maximization problem, the form of the simplex array and its optimization

%% Preparation for script execution, pre-cleaning console, variables, open windows
clear; close all; clc; 

%% Sample parameters

NoOfX = 3; % Number of decision variables in the problem
Z = [-1 3 2 6]; % Vector of objective function weights (coefficient values)
a = [3 -1 2; -2 4 0;-4 3 8]; % Bounding matrix (values of the coefficients for the bounding conditions)
b = [7; 12; 10]; % Vector of free expressions
% Types of inequalities for the boundary conditions of the problem
% Assume that: 1 means ">=", 0 means "=", -1 means "<=".
TypesCon = [0 -1 1];
% Types of boundings on the decision variables of the problem
% Assume that: 1 means ">= 0", 0 means "any", -1 means "<= 0".
TypesX = [1 0 -1];

AreFreeX = (NoOfX~=length(Z)); % Auxiliary information whether there is a free expression in the coefficients of the objective function

%% Creation of the canonical form for the maximization problem - the first stage of substitution
if AreFreeX
    LeftMatrix=[-Z;a,zeros(NoOfX,1)]; % Convertion to a minimization problem (reversing the sign of the value of the objective function)
else
    LeftMatrix=[-Z;a];
end
% Removing inequalities - creating free variables with the help of a "adjugate" matrix
% The numbers: -1,1 indicate here the added free variables for the inequality.
% They have no effect on the result of the objective function (values of its last coefficients equal zero). 
% The first row of the "adjugate" matrix indicates the coefficients of the objective function.
a2 = zeros(size(a,1)+1,size(a,1)); % First we create a null matrix
for i=1:length(TypesCon) % Assigning row/column values to added free variables
    if(TypesCon(i)==1 || TypesCon(i)==-1)
        a2(i+1,i) = - TypesCon(i);
    end
end

% We will search the columns of the resulting matrix to find those without added free variables and remove them
i=size(a2,2);    
while i>=1
    if (a2(i+1,i)==0)
        a2(:,i) = [];
        i=i-1;
    end
    i=i-1;
end

% Create text fields for matrix with structure and simplex array information
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

% Conversion of negative decision variables to positive ones
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
fprintf('\t The resulting canonical matrix (after the second substitution step) is: \n');
disp(CanonicalArr); % The resulting figure after the second substitution step
% 
% % Calling a ready-made function for specified parameters
% % SimplexArray(Z,a,b,TypesBon,TypesX)
% % Z - vector of objective function weights (coefficient values)
% % a - bounding matrix (values of the coefficients for the bounding conditions)
% % b - vector of free expressions
% % Assume that: 1 means ">=", 0 means "=", -1 means "<=":
% % TypesBon - inequality types for bounding conditions
% % TypesX - types of boundings for the decision variables of the problem
% CanonicalArr = CanonicalArray(3,[-1 3 2 6],[3 -1 2; -2 4 0;-4 3 8],[7; 12; 10],[0 -1 1],[1 0 -1]);

%% Creating a simplex array from the canonical form
% In the first case: the first column indicates the free expressions of the objective and constraint functions
if AreFreeX
    Simplex = [[CanonicalArr(1,end-2);CanonicalArr(2:end,end)],CanonicalArr(:,2:end-3)];
else
    Simplex = [[0;CanonicalArr(2:end,end)],CanonicalArr(:,2:end-2)];
end
Simplex=str2double(Simplex);
fprintf('\t The resulting simplex array from the canonical form is: \n');
disp(Simplex); % The resulting figure after the first substitution step

% % Calling a ready-made function for a specified canonical form
% Simplex = SimplexArray(CanonicalArr,AreFreeX);

%% Contradiction check, creation of unitary matrix
% Conflict check of type I boundary conditions (b_i~=0 ^ a_i=0)
IsContr=0;
for i=2:size(Simplex,1)
    if (Simplex(i,1)~= 0 && ~any(Simplex(1,2:end)~=0))
        IsContr=1;
        fprintf ('\t This array is in contradictory form of the 1st kind ! (Line(s) no.:');
        disp(i); fprintf('\n No acceptable solutions for this array.\n');
    end
end

% Conflict check of type II boundary conditions (b_i<0 ^ a_i >=0 or b_i>0 ^ a_i <=0)
for i=2:size(Simplex,1)
    if ((Simplex(i,1) > 0 && ~any(Simplex(i,2:end)>0)) || (Simplex(i,1)< 0 && Simplex(i,2:end)>=0))
        IsContr=2;
        fprintf ('\t This array is in contradictory form of II kind ! (Line(s) no.:');
        disp(i); fprintf('\n No acceptable solutions for this array.\n');
    end
end
if IsContr==0
    fprintf ('\t This array is not in contradictory form of the I, or II kind. \n\n');
end
  
% % Calling a ready-made function for the specified simplex array to check for contradictions
% IsContr = IsContradiction(Simplex);
% if IsContr==1
%     fprintf ('\t This array is in contradictory form of type I !')
%     fprintf('\t No solutions allowed for this array !');
% else 
%     if IsContr==2
%         fprintf ('\t This array is of contradictory form II type !');
%         fprintf('\t No solutions allowed for this array !');
%     else
%         fprintf ('\t This array is not in contradictory form of the I, or II kind.');
%     end
% end

% Check for existence of unit matrix and zero coefficients of the objective function for the baseline variables
 if IsContradiction(Simplex)==0
     % Below is a unitary array with as many rows as there are in the simplex array with zeros in the first row
     I = [zeros(1,size(Simplex,1)-1);eye(size(Simplex,1)-1)]; UnitMatrixNr = cell(1, size(Simplex,1)-1);
     for i=1:size(Simplex,1)-1 % Finding the column numbers corresponding to the unit matrix
         UnitMatrixNr{i} = find(ismember(Simplex',I(:,i)','rows'));
     end
     UnitMatrix_Columns=int8.empty; i=1;
     while i<=size(Simplex,1)-1
         if isempty(UnitMatrixNr{i})
             IsEyeMat = false;
             break;
         else
             UnitMatrix_Columns(i) = UnitMatrixNr{i};
         end
         if i == size(Simplex,1)-1
             IsEyeMat = true;
         end
         i=i+1;
     end
 else
     if IsContradiction(Simplex)==1
         fprintf('\t This array is in contradictory form of I kind ! There is no need to look for a unitary array. \n\n');
     else
         fprintf('\t This array is in contradictory form of the II kind ! There is no need to look for a unitary array. \n\n');
     end 
 end

% % Call a ready-made function for the specified symplex array and column information of the unit array
% fprintf('We will now create a unit array in the array') ;   
% IsEyeMat = IsEyeMatrix(Simplex);
% if IsEyeMat==0
%     fprintf('Not enough columns of the unit array in the array');
% else
%     fprintf('The array has a unit matrix. The numbers of the baseline variables are: ');
%     disp(UnitMatrix_Columns);
% end     

% % Example of how the function works for other parameters of the simplex array
% Simplex=[0 1 3 2 -1 1 0;
%     -10 0 -1 0 0 0 1;
%     -15 1 -2 2 -2 0 1;
%     -20 1 -2 3 -3 1 0];
% IsEyeMat = IsEyeMatrix(Simplex);
% if IsEyeMat==0
%     fprintf('\t Not enough unit array columns in the array.\n\n');
% else
%     fprintf('\t The array has a unitary array. The numbers of the baseline variables are: ');
%     disp(UnitMatrix_Columns);
% end   

%% Creating a unit array in a simplex array
% We look for the first non-zero element in the following rows of the bounding conditions.
% We multiply its row by an appropriate number so that the element equals: 1.
% Then we do rotations of the array around the given row until we get a unitary array.
% I.e., zero the remaining values in the element column before adding the element rows accordingly

if IsContradiction(Simplex)==0
    fprintf('\t We will now create a unit array in an array. \n');
    DRt=0; % Done rotations - the number of elements around which rotations have already been performed
    while IsEyeMatrix(Simplex) == false % As long as there are no columns of the unit matrix
        [col,row] = find(Simplex(2+DRt:end,2:end)'); % Searching for a non-zero element
        ElementValue = Simplex(row(1)+1+DRt,col(1)+1); % Value of non-zero element
        Simplex(row(1)+1+DRt,:) = Simplex(row(1)+1+DRt,:)/ElementValue; % The element found is replaced by 1
        
        for i=1:size(Simplex,1) % To each of the other rows
            if i == row(1)+1+DRt
                continue;
            else % Add the appropriate number of times to the rotation row to reset the elements in the rotation column
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
        fprintf('\t This array is in contradictory form of I kind ! There is no need to create a unitary array.. \n\n');
    else
        fprintf('\t This array is in contradictory form of II kind ! There is no need to create a unitary array. \n\n');
    end
end
% Removing unnecessary (null) lines
DRm=0; % Number of rows deleted
for i=size(Simplex,1):-1:1 
        if Simplex(i,:)==0
            Simplex = [Simplex(1:i-1,:);Simplex(i+1:end,:)]; DRm=DRm+1; 
            fprintf('\t The array after removing the redundant row no.%d looks as follows: \n',DRm);
            disp(Simplex);
        end
end
if IsContradiction(Simplex)==0
    fprintf('\t A simplex array in "near-baseline" form with unit array columns is: \n');
    disp(Simplex);
    fprintf('\t The numbers of the "near baseline" variables are: \n');
    [~, UnitMatrixColumns]=IsEyeMatrix(Simplex);
    disp(UnitMatrixColumns);
end 

% % Calling a predefined function for a specified simplex array to create a unit array 
% Simplex = CreateEyeMatrix(Simplex);
% if IsContradiction(Simplex)==0
%    fprintf('\t A simplex array in "near-baseline" form with unit array columns is: \n');
%    disp(Simplex);
%    fprintf('\t The numbers of the "near baseline" variables are: \n');
%    [~, UnitMatrixColumns]=IsEyeMatrix(Simplex);
%    disp(UnitMatrixColumns);
% end 


%% Checking the baseline of the simplex array, eliminating free negative elements in the bounding conditions
% Including the technique of creating subproblems 

if IsContradiction(Simplex)==0
    fprintf('\t We will now create an array in base form. \n');
    % When any free element (except the objective function) is negative, perform appropriate array rotations to obtain the base array
    IsBaseForm = false; DRt = 0;
    while ~IsBaseForm
        if (IsContradiction(Simplex)==0)
            IsBaseForm = true; % Variable indicating whether array is in base form, loop end condition
            for i=2:size(Simplex,1) 
                if Simplex(i,1)<0 % When there is any negative free expression in the objective function
                    IsBaseForm = false;
                    % Looking for another negative value on its row to select the rotation column
                    [Value,MinCol] = min(Simplex(i,2:end)); MinCol=MinCol+1;
                    
                    % Next, we select a pivot row:
                    % we look for the minimum positive quotient of a column element by the corresponding free expression (we skip the given row)
                    Quotients = [Simplex(1:i-1,1);nan;Simplex(i+1:end,1)] ./ [Simplex(1:i-1,MinCol);nan;Simplex(i+1:end,MinCol)];
                    Quotients(Quotients<=0) = nan; 
                    MinQuotient = min(Quotients);
                    MinRow = find(Quotients==MinQuotient,1);
                    Min = Simplex(MinRow, MinCol);
                    if (isnan(Min) || MinQuotient == Inf)
                        Min = Value; 
                        MinRow = i;
                    end
                    
                    % In the description the situation when each free word in the bounding condition is negative was omitted.
                    % In this case, we can choose any row of the rotation, column to satisfy the quotient condition.
                    
                    % Now perform a rotation around a given row and column on the entire array
                    Simplex(i,:) = Simplex(i,:)/Min; % The element found is replaced by 1
                    for j=1:size(Simplex,1) % In each of the remaining rows, zero the elements on that column
                        if j == MinRow
                            continue;
                        else
                            % Add the appropriate number of times to the rotation row to reset the elements in the rotation column
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
    fprintf('\t The simplex array in base form without negative free expressions is: \n'); disp(Simplex);
end
% 
% % Calling a predefined function for a specified simplex array to create a base form
% Simplex = MakeBaseForm(Simplex); disp(Simplex);

%% Checking for the possible existence of an unbounded form in the array / Determining the optimal solution

% Check that the array is not in unbounded form
IsIndefinite=false;
if IsContradiction(Simplex)==1
        fprintf('\t This array is in contradictory form of I kind !\n\n');
else
    if IsContradiction(Simplex)==2
        fprintf('\t This array is in contradictory form of II kind !\n\n');
    else
        for i=2:size(Simplex,2)
            if Simplex(1,i)<0 % When there is any negative free expression in the objective function
                if Simplex(:,i)<=0 % When each element of a column is less than or equal to zero
                    IsIndefinite = true; % Unbounded form
                    fprintf('\t This array is in unbounded form ! (Check column no.%d \n\n',i);
                    break;
                end
            end
        end
        if IsIndefinite==false
            fprintf('\t This array is not in unrestricted form.\n\n');
        end
    end
end

% % Calling a predefined function for a specified simplex array to check the unbounded form in the array
% if IsIndefiniteMatrix(Simplex)
%     fprintf('\t The simplex array is in unbounded form ! \n');
% else
%     fprintf('\t The simplex array is not in unbounded form. \n');
% end

% For example:
Simplex=[-55 0 0 4 -3 5 0;
    15 0 1 -1 1 -1 0;
    10 1 0 1 -1 -1 0;
    5 0 0 -1 1 -1 1];

% Determination of the optimal array (symplectic rotation)
DRt=0;
Simplex2=Simplex;
while max(~(Simplex(1,2:end)>=0)) % As long as there is a negative coefficient in the objective function
    % For faster optimization, we search for the minimum element in the objective function
    ColNr = find(Simplex(1,2:end)==min(Simplex(1,2:end)),1)+1; % Rotation column selection
    CCQuotient = Simplex(2:end,1) ./ Simplex(2:end,ColNr);
    RowNr = find(CCQuotient == min(CCQuotient(CCQuotient>0)))+1; % Rotation row selection
    MinEl = Simplex(RowNr,ColNr);
    
    % Now perform a rotation around a given row and column on the entire array
    Simplex(RowNr,:) = Simplex(RowNr,:)/MinEl; % The element found is replaced by 1
    for j=1:size(Simplex,1) % In each of the remaining rows, zero the elements on the rotation column
        if j == RowNr
            continue;
        else
            % Dodajemy odpowiedni¹ iloœæ razy wiersz obrotu, aby wyzerowaæ elementy w kolumnie obrotu
            % Add the appropriate number of the rotation row to zero the elements in the rotation column
            Simplex(j,:) = Simplex(j,:) - Simplex(j,ColNr) * Simplex(RowNr,:); 
        end
    end
    DRt=DRt+1;
    fprintf('\t The array after rotating array no.%d looks as follows: \n',DRt);
    disp(Simplex);
end
fprintf('\t The above array is in optimal form. \n');
    
% % Calling a ready function to determine the optimal array
% Simplex = MakeOptimalForm(Simplex);
% fprintf('\t The array in its optimal form looks as follows: \n'); disp(Simplex); 

% Assign the obtained values of the optimal solution to the relevant variables
I = [zeros(1,size(Simplex,1)-1);eye(size(Simplex,1)-1)]; 
UnitMatrixNr=int16.empty;
for i=1:size(Simplex,1)-1 % Finding the column numbers corresponding to the unit matrix
    UnitMatrixNr(i) = find(ismember(Simplex',I(:,i)','rows'));
end
XValues = cell(1, size(Simplex,2)-1);
for i=1:length(UnitMatrixNr)
    XValues{UnitMatrixNr(i)-1} = Simplex(1+i,1); % Assigning values to items in unit columns
end
for i=1:length(XValues)
    if ~isempty(XValues{i}) == 0 
        XValues{i}=0; % Value assignment: zero for all other elements
    end
end
Vec=[1:length(XValues); cell2mat(XValues)]; Vec=Vec(:)';
fprintf('\t The obtained values of the optimal solution for the following variables are: \n'); fprintf('\t x%d: %d\n',Vec); 

% % Calling ready function to display optimal solution from optimal array
% Vec=[1:length(OptimalValues(Simplex)); cell2mat(OptimalValues(Simplex))]; Vec=Vec(:)';
% fprintf('\t Wartoœci otrzymanych zmiennych to: \n'); fprintf('\t x%d: %d\n',Vec); 

%% Checking another optimal solution of the given simplex array (we assume that no radius of optimality arises)

% For example, we can check the following matrix:
Simplex=[30 0 3 0 0 2 4 0;
    10 1 -1 0 0 1 2 1;
    20 0 1 0 1 -1 0 1;
    15 0 2 1 0 1 1 -1]; disp (Simplex);

% Simplex=[3 0 0 1 0 0;
%     1 1 1 3 0 -2;
%     2 0 -2 -1 1 1]; disp (Simplex);

Simplex = CreateEyeMatrix(Simplex); % We will use the previously created function to find the columns of the unit matrix
[~, ~, UnitMatrixColumns_2]=IsEyeMatrix(Simplex);
if length(find(Simplex(1,2:end)==0))>length(UnitMatrixColumns_2) % Check for non-baseline zeros in the objective function
    NonBaseX=setdiff(find(Simplex(1,2:end)==0)+1,str2double(UnitMatrixColumns_2));
    % We assume no optimality radius, i.e., non-base columns contain at least one value greater than zero
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
                    RowEl = find(CCQuotient == MinQuo)+1;  % Rotation row selection
                    Simplex = Rotate(Simplex,RowEl,NonBaseX(j+1)); % Array rotation
                    j=j+1; Stop = ~isempty(RowEl); % Condition to stop the loop when a rotation row is found
                end
            end
        end
        if Text
            Vec=[1:length(OptimalValues(Simplex)); cell2mat(OptimalValues(Simplex))]; Vec=Vec(:)';
            fprintf('\t The array with the other optimal solution looks as follows: \n'); disp(Simplex);
            fprintf('\t The values of the resulting variables are: \n'); fprintf('\t x%d: %d\n',Vec); 
        end    
    end
end

% % Call a ready-made function to check another optimal solution to a given symplectic array
% % (we assume that no optimality radius is created)
% Simplex = OtherResult(Simplex);

%% 
% There is a possibility to develop the project into a set of optimal
% solutions (straight line or geometric figure).