%% The simplex method of linear programming - implementation in Matlab language
% Including: the canonical form of the maximization problem, the form of the simplex array and its optimization

%% Preparation for script execution, pre-cleaning console, variables, open windows
clear; close all; clc; 

%% Sample symplex arrays
% Simplex=[-6     1    -3     3     2     0     0;
%      7     3    -1     1    -2     0     0;
%     12    -2     4    -4     0     1    0;
%     10    -4     3    -3    -8     0    -1];

% Below examples of symplex arrays:
% Simplex=[0 -2 -1 0 0 0;
%      14 2 2 1 0 0;
%      20 4 0 0 1 0;
%     18 2 3 0 0 1]; 
% Simplex=[15 0 0 0 -2 1;
%      -10 1 0 0 2 -3;
%      10 0 1 0 -1 2;
%      4 0 0 1 1 1]; 
% Simplex=[15 0 0 0 -2 1;
%      -10 1 0 0 -1 -2;
%      10 0 1 0 -1 2;
%      4 0 0 1 1 1];
% Simplex=[0 1 3 2 -1 1 0;
%      -10 0 -1 0 0 0 1;
%      -15 1 -2 2 -2 0 1;
%      -20 1 -2 3 -3 1 0];
% Simplex=[30 0 3 0 0 2 4 0;
%      10 1 -1 0 0 1 2 1;
%      20 0 1 0 1 -1 0 1;
%      15 0 2 1 0 1 1 -1];
% Simplex=[-5 1 0 0 0;
%      2 1 0 1 -2;
%      1 0 1 -2 0];
% Simplex=[3 0 0 1 0 0;
%      5 1 -3 1 2 0;
%      2 0 -2 -1 1 1];
% We will get the first 2 solutions obtained manually, after another use of "OtherResult" function we will find the third one.
 Simplex=[-7 0 0 1 0 0 0;
     15 1 1 1 0 0 -1;
     5 0 -1 0 0 1 0;
     10 0 0 -1 1 0 -1];
% We will get 2 of the 3 solutions because the "OtherResult" function will return the first one when used again.
% It selects the first possible element to rotate the array.

if exist('Simplex','var')
    fprintf('\t The resulting simplex array is: \n'); disp(Simplex);
else
fprintf('\t No simplex array chosen ! \n'); return;
end

%% Array optimization

% Contradiction check, creation of unitary matrix
IsContr = IsContradiction(Simplex);
if IsContr==1
    fprintf ('\t The resulting array is in contradictory form of I kind !');
    fprintf('\n No acceptable solutions for this array.\n');
else
    if IsContr==2
        fprintf ('\t The resulting array is in contradictory form of II kind !');
        fprintf('\n No acceptable solutions for this array.\n');
    else
        fprintf ('\t This array is not contradictory. \n');
        
        % Creating a unit array in a symplex array
        Simplex2 = CreateEyeMatrix(Simplex);
        if IsContradiction(Simplex2)==0
            if (Simplex2==Simplex)
                fprintf('\t Array above is a "near-baseline" simplex array. \n');
            else
                fprintf('\t Above, a "near-baseline" simplex array with unit array columns was created. \n');
            end
            Simplex=Simplex2;
            fprintf('\t The numbers of the "near baseline" variables are: \n');
            [~, UnitMatrixColumns] = IsEyeMatrix(Simplex);
            disp(UnitMatrixColumns);
            
            % Checking the baseline of the simplex array, eliminating free negative elements in the bounding conditions
            Simplex = MakeBaseForm(Simplex);
            
            % Checking an unbounded form in an array (only if it does not contradict)
            if IsContradiction(Simplex)~=0
                if IsIndefiniteMatrix(Simplex)
                    fprintf('\t The simplex array is in unbounded form ! \n');
                end
            else
                fprintf('\t The simplex array is not in unbounded form. \n');
                    
                % Determination of the optimal array (simplex rotation)
                Simplex2 = MakeOptimalForm(Simplex);
                if (Simplex2==Simplex)
                    fprintf('\t But the array above is in optimal form, so no actions needed. \n');
                else
                    fprintf('\t The above array is in optimal form. \n');
                end
   
                % Displaying the optimal solution
                Vec=[1:length(OptimalValues(Simplex)); cell2mat(OptimalValues(Simplex))]; Vec=Vec(:)';
                fprintf('\t The values of the resulting variables are: \n'); fprintf('\t x%d: %d\n',Vec); 

                % Check for another optimal solution of the given simplex array (we assume no radius of optimality)
                Simplex = OtherResult(Simplex);
            end
        end
    end
end
