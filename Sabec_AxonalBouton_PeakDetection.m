%% IMPORT EXCEL DATA AND TRANSFORM 
% Define the filename - not from Onedrive
filename = ['C:\Users\ms0405\Suite2P\Traces\Rep trace.xlsx'];

% Specify the range to import 
range = 'F2:AHU350'; % Adjust if over 350 boutons in FOV

% Extract the base filename (without path and extension) - Help track files
[~, baseFileName, ~] = fileparts(filename);

% Create a valid MATLAB variable name
validVarName = matlab.lang.makeValidName(baseFileName);

% Check if the file exists
if exist(filename, 'file') == 2
    try
        % Read data as a matrix from the specified range
        dataMatrix = readmatrix(filename, 'Range', range);
        
        % Remove rows that are entirely NaN
        dataMatrix = dataMatrix(~all(isnan(dataMatrix), 2), :);

        % Transpose the matrix so that rows correspond to original columns
        dataTransposed = dataMatrix';
        
        % Dynamically assign the transposed matrix to the variable name
        eval([validVarName ' = dataTransposed;']);
        disp(['Data imported and transposed into variable: ' validVarName]);
        disp(dataTransposed);
    catch err
        % If readmatrix fails, display the error message
        disp('Failed to import data as a matrix:');
        disp(err.message);
    end
else
    disp('File does not exist.');
end
%% DRIFT CORRECTION

% Compute the moving average for each column
windowSize = 25;  % Adjust window size as needed
numCols = size(dataTransposed, 2);  % Number of columns

% Initialize corrected data matrix
data_corrected = zeros(size(dataTransposed));

% Loop through each column to apply moving average
for col = 1:numCols
    % Compute the moving average for the current column
    smoothedData = movmean(dataTransposed(:, col), windowSize);
    
    % Correct the drift by subtracting the smoothed data
    data_corrected(:, col) = dataTransposed(:, col) - smoothedData;
end

% Plot original and corrected data for the first column as an example
figure;
plot(dataTransposed(:, 1), 'b'); % Original data
hold on;
plot(data_corrected(:, 1), 'r'); % Corrected data
legend('Original Data', 'Corrected Data');
xlabel('Frames');
ylabel('dF/F');
title('Drift Correction Using Moving Average');
hold off;

% Dynamically assign the transposed matrix to the variable name
        eval([validVarName ' = data_corrected;']);
        disp(['Data imported and transposed into variable: ' validVarName]);
        disp(data_corrected);


%% DETECT TRANSIENT PEAKS AND LOCATIONS
% Assuming the variable was created and is a matrix, apply findpeaks to each row
if exist(validVarName, 'var')
    % Extract the variable from the dynamically named matrix
    data_corrected = eval(validVarName);
    
    % Check if the imported data is numeric
    if isnumeric(data_corrected)
    
        % Get the number of columns (original columns after transposition)
        numColumns = size(data_corrected, 2);
        
        % Initialize cell arrays to store peaks, locations, and standard deviations for each column
        allPks = cell(numColumns, 1);
        allLocs = cell(numColumns, 1);
        allStdDevs = cell(numColumns, 1);
        
        % Loop through each column
        for col = 1:numColumns
            % Extract the data from the current column (row in the transposed matrix)
            columnData = data_corrected(:, col);
           
            % Compute the standard deviation of the column data
            stdDev = std(columnData, 'omitnan');
            allStdDevs{col} = stdDev;
            pop_stdDev = mean(cell2mat(allStdDevs));
        end 
        
        for col = 1:numColumns

            columnData = data_corrected(:, col);

            % Display column data for debugging
            disp(['Processing Column ' num2str(col) ':']);
            disp(columnData);
        
            % Apply findpeaks with a threshold based on the standard deviation
            [pks, locs] = findpeaks(columnData, 'MinPeakProminence',1.0, 'MinPeakWidth',1.0, 'MinPeakHeight',2.5*pop_stdDev);

            % Display peaks and locations for debugging
            disp(['Peaks in Column ' num2str(col) ':']);
            disp(pks);
            disp(['Locations in Column ' num2str(col) ':']);
            disp(locs);
            
            % Store the peaks and locations in cell arrays
            allPks{col} = pks;
            allLocs{col} = locs;
            
            % Plot the data and the peaks for the current column
            figure;
            plot(columnData);
            hold on;
            plot(locs, pks, 'rv', 'MarkerFaceColor', 'r');
            title(['Peaks in Bouton ' num2str(col)]);
            xlabel('Frames');
            ylabel('dF/F');
            hold off;
        end
        
        % Construct the MAT-file name using baseFileName
        matFileName = [baseFileName, '_results.mat'];
        
        % Save the results
        %save(matFileName, 'allPks', 'allLocs', 'allStdDevs');
 

        
     
        
%% CALCULATE AVERAGE TRANSIENT AMPLITUDE FOR EACH BOUTON 

% Get the number of columns in allPks
numColumns = length(allPks);

% Initialize an array to store the average values
averagePks = NaN(numColumns, 1); % Using NaN to handle cases with no peaks

% Loop through each cell in allPks
for col = 1:numColumns
    % Extract the data from the current cell
    pks = allPks{col};
    
    % Compute the average if the cell is not empty
    if ~isempty(pks)
        averagePks(col) = mean(pks);
    else
        averagePks(col) = NaN; % Assign NaN if no data is present
    end
end

% Display the results
disp('Average values for each cell in allPks:');
disp(averagePks);


% Construct the MAT-file name with prefix based on baseFileName
matFileName = [baseFileName, '_average_peaks.mat'];

% Save the results to the MAT-file with the dynamic name- TURN ON IF FILES REQ'D
%save(matFileName, 'averagePks');
%disp(['Average peaks saved to ' matFileName]);



%% CALCULATE NUMBER AND FREQUENCY OF TRANSIENTS DETECTED


% Define recording length in seconds
Time = 232.5; % standard length for 900 frames at 0.388Hz

% Get the number of columns in allLocs
numColumns = length(allLocs);

% Initialize an array to store the counts of values
countLocs = NaN(numColumns, 1); % Using NaN to handle cases with no locations

% Initialize an array to store the counts divided by time
Frequency = NaN(numColumns, 1);

% Loop through each cell in allLocs
for col = 1:numColumns
    % Extract the data from the current cell
    locs = allLocs{col};
    
    % Compute the count of values if the cell is not empty
    if ~isempty(locs)
        countLocs(col) = numel(locs);
        % Divide the count by time
        Frequency(col) = countLocs(col) / Time;
    else
        countLocs(col) = 0; % Assign 0 if no data is present
        Frequency(col) = 0; % Assign 0 if no data is present
    end
end

% Display the results
disp('Count of values for each cell in allLocs:');
disp(countLocs);

disp('Counts divided by total number of frames:');
disp(Frequency);

% Construct the MAT-file name with prefix based on baseFileName
matFileName = [baseFileName, '_Frequency.mat'];

% Save the results to the MAT-file with the dynamic name - TURN ON IF FILES REQ'D
%save(matFileName, 'countLocs', 'Frequency');
%disp(['Frequency saved to ' matFileName]);


%% CALCULATE AREA UNDER THE CURVE (TRANSIENT ONLY) PER FOV

% Initialize cell array to store the AUCs for each column
allAUCs = cell(numColumns, 1);

% Loop through each column
for col = 1:numColumns
    % Extract the data from the current column
    columnData = data_corrected(:, col);
    
    % Initialize array to store AUCs for the current column
    aucs = zeros(length(allLocs{col}), 1);
    
    % Loop through each detected transient to calculate AUC
    for i = 1:length(allLocs{col})
        % Get the location of the peak
        loc = allLocs{col}(i);
        
        % Define the window around the peak for integration
        % You may adjust the window size to best capture the transient
        startIdx = max(1, loc - 5);  % Start 5 frames before the peak
        endIdx = min(length(columnData), loc + 5);  % End 5 frames after the peak
        
        % Extract the transient window
        transientWindow = columnData(startIdx:endIdx);
        
        % Compute the AUC using the trapezoidal method
        aucs(i) = trapz(transientWindow);
    end
    
    % Store the AUCs in the cell array
    allAUCs{col} = aucs;
end

% Initialize array to store the sum of AUCs for each column
sumAUCs = zeros(numColumns, 1);

% Loop through each column to sum the AUCs
for col = 1:numColumns
    % Sum the AUCs for the current column
    sumAUCs(col) = sum(allAUCs{col});
end

% Calculate mean AUC and stnadard devaition across columns/boutons
avgAUC = mean(sumAUCs, 'omitnan');

stdAUC = std(sumAUCs, 'omitnan');



%% COMPILE TRANSIENT AMPLITUDE STATS_FOV

% Transient amplitude dF/F
% Check if averagePks is a numeric matrix
if isnumeric(averagePks)
    % Calculate the average for each column
    avgValuesAmp = mean(averagePks, 'omitnan');
    
    % Calculate the standard deviation for each column
    stdValuesAmp = std(averagePks, 'omitnan');
    
    % Calculate the number of samples for each column
    numSamplesAmp = sum(~isnan(averagePks)); % Count non-NaN values in each column
   
     % Calculate the total number of elements
    totalElements = numel(averagePks);
    
    % Calculate the number of non-NaN elements
    nonNaNElements = sum(~isnan(averagePks(:)));
    
    % Calculate the percentage of non-NaN values
    percentageNonNaN = (nonNaNElements / totalElements) * 100;
    
    % Display the results
    disp('Average values for each column in averagePks:');
    disp(avgValuesAmp);
    
    disp('Standard deviation for each column in averagePks:');
    disp(stdValuesAmp);
    
    disp('Number of samples for each column in averagePks:');
    disp(numSamplesAmp);
    
    disp('Percentage of non-NaN values in averagePks:');
    disp(percentageNonNaN);

    % Construct the MAT-file name with prefix based on baseFileName
    %filename = 'C:\Users\ms0405\150356_240618_ROI1_dF.xlsx';
    [~, baseFileName, ~] = fileparts(filename);
    matFileName = [baseFileName, '_statistics.mat'];
    
    % Save the results to the MAT-file with the dynamic name- TURN ON IF FILES REQ'D
    %save(matFileName, 'avgValuesAmp', 'stdValuesAmp', 'numSamplesAmp');
    %disp(['Statistics saved to ' matFileName]);
else
    disp('averagePks is not a numeric matrix.');
end

% Transient frequency stats

% Calculate the average for each column
    avgValuesFreq = mean(Frequency, 'omitnan');
    
    % Calculate the standard deviation for each column
    stdValuesFreq = std(Frequency, 'omitnan');
    
    % Calculate the number of samples for each column
    numSamplesFreq = sum(~isnan(Frequency)); % Count non-NaN values in each column
    
       
    % Display the results
    disp('Average values for each column in Frequency:');
    disp(avgValuesFreq);
    
    disp('Standard deviation for each column in averagePks:');
    disp(stdValuesFreq);
    
    disp('Number of samples for each column in averagePks:');
    disp(numSamplesFreq);
    
    
    % Construct the MAT-file name with prefix based on baseFileName
    %filename = 'C:\Users\ms0405\150356_240618_ROI1_dF.xlsx';
    [~, baseFileName, ~] = fileparts(filename);
    matFileName = [baseFileName, '_freq_statistics.mat'];
    
    % Save the results to the MAT-file with the dynamic name - TURN ON IF FILES REQ'D
    %save(matFileName, 'avgValuesFreq', 'stdValuesFreq', 'numSamplesFreq');
    %disp(['Statistics saved to ' matFileName]);

%% EXPORT ALL STATS AS EXCEL TABLE_FOV 

% Check if avgValuesAmp, stdValuesAmp, numSamplesAmp, avgValuesFreq, stdValuesFreq, numSamplesFreq, avgAUC, and stdAUC are numeric vectors
if isnumeric(avgValuesAmp) && isnumeric(stdValuesAmp) && isnumeric(numSamplesAmp) && isnumeric(avgValuesFreq) && isnumeric(stdValuesFreq) && isnumeric(numSamplesFreq) && isnumeric(avgAUC) && isnumeric(stdAUC)
    % Ensure all vectors are row vectors
    avgValuesAmp = avgValuesAmp(:)';
    stdValuesAmp = stdValuesAmp(:)';
    numSamplesAmp = numSamplesAmp(:)';
    avgValuesFreq = avgValuesFreq(:)';
    stdValuesFreq = stdValuesFreq(:)';
    numSamplesFreq = numSamplesFreq(:)';
    avgAUC = avgAUC(:)'; % Ensure avgAUC is a row vector
    stdAUC = stdAUC(:)'; % Ensure stdAUC is a row vector
    
    % Combine the vectors into a single matrix
    combinedMatrix = [avgValuesAmp; stdValuesAmp; numSamplesAmp; avgValuesFreq; stdValuesFreq; numSamplesFreq; avgAUC; stdAUC];
    
    % Display the combined matrix
    disp('Combined matrix (avgValuesAmp, stdValuesAmp, numSamplesAmp, avgValuesFreq, stdValuesFreq, numSamplesFreq, avgAUC, stdAUC):');
    disp(combinedMatrix);

    % Transpose the matrix so that rows correspond to original columns
    CombinedStatsTransposed = combinedMatrix';

    % Define the labels for each column
    columnLabels = {'Av_Amplitude_df/DF', 'SD_Amplitude_df/DF', 'n_Active', 'Av_Frequency_Hz', 'SD_Frequency_Hz', 'n_Total', 'Av_AUC', 'SD_AUC'};
    
    % Convert the matrix to a table and assign column labels
    CombinedStatsTable = array2table(CombinedStatsTransposed, 'VariableNames', columnLabels);

    % Display the table
    disp('Table with column labels:');
    disp(CombinedStatsTable);

    % Define the new column data and label
    newColumnLabel = 'File_ID';
    newColumnData = repmat({baseFileName}, height(CombinedStatsTable), 1);

    % Add the new column to the table
    CombinedStatsTable.(newColumnLabel) = newColumnData;

    % Display the updated table
    disp('Updated Table with File_ID column:');
    disp(CombinedStatsTable);

    % Define the Excel file name with the same prefix
    excelFileName = [baseFileName, '_FOV_statistics.xlsx'];

    % Export the table to the Excel file
    writetable(CombinedStatsTable, excelFileName);

    % Display a message confirming the export
    disp(['Table exported to ' excelFileName]);
else
    disp('One or more input vectors are not numeric. Please check the data types.');
end


%% EXPORT INDIVIDUAL TRANSIENTS TO EXCEL
% Initialize an empty array to collect concatenated data
concatenatedAmp = [];

% Loop through each cell in the cell array
for k = 1:numel(allPks)
    % Concatenate vertically
    concatenatedAmp = [concatenatedAmp; allPks{k}(:)];
end

% Display the result
disp('Single column matrix of allPks data:');
disp(concatenatedAmp);

% Define the labels for each column
columnLabel = {baseFileName};
    
% Convert the matrix to a table and assign column labels
AllAmpTable = array2table(concatenatedAmp, 'VariableNames', columnLabel);

% Display the table
disp('Table with column label:');
disp(AllAmpTable);

% Construct the MAT-file name with prefix based on baseFileName
%filename = 'C:\Users\ms0405\150356_240618_ROI1_dF.xlsx';
[~, baseFileName, ~] = fileparts(filename);
matFileName = [baseFileName, '_all_transients.mat'];

% Save the table to a MAT-file with the dynamic name
%save(matFileName, 'AllAmpTable');
%disp(['Table saved to ' matFileName]);

% Define the Excel file name with the same prefix
excelFileName = [baseFileName, '_all_transients.xlsx'];

% Export the table to the Excel file
writetable(AllAmpTable, excelFileName);

% Display a message confirming the export
disp(['Table exported to ' excelFileName]);

%% CALCULATE SPIKE/TRANSIENT SYNCHRONICITY

% Generate raster plot to see population firing patterns
numColumns = length(allLocs);

figure;
hold on;

for col = 1:numColumns
    locs = allLocs{col};
    for peak = 1:length(locs)
        plot([locs(peak), locs(peak)], [col-0.4, col+0.4], 'k');
    end
end

xlabel('Frames');
ylabel('Boutons');
title('Raster Plot of Ca Transients');
set(gca, 'YDir', 'reverse');
yticks(1:numColumns);
yticklabels(1:numColumns);
hold off;

% Calculate bouton-bouton correlations
maxFrames = size(data_corrected, 1);
numColumns = length(allLocs);
eventMatrix = zeros(maxFrames, numColumns);

for col = 1:numColumns
    if ~isempty(allLocs{col})
        eventMatrix(allLocs{col}, col) = 1;
    end
end

correlationMatrix = corr(eventMatrix);

% Replace NaN values with 0
correlationMatrix(isnan(correlationMatrix)) = 0;

% Only keep positive values and scale them between 0 and 1
positiveCorrelations = correlationMatrix;
positiveCorrelations(positiveCorrelations < 0) = 0;

% Rescale positive values between 0 and 1
maxValue = max(positiveCorrelations(:));
minValue = min(positiveCorrelations(:));

% Avoid division by zero in case all values are the same
if maxValue ~= minValue
    scaledCorrelationMatrix = (positiveCorrelations - minValue) / (maxValue - minValue);
else
    scaledCorrelationMatrix = positiveCorrelations; % No scaling needed
end

disp('Scaled Correlation Matrix (Positive Values Only):');
disp(scaledCorrelationMatrix);

% Plot the scaled correlation matrix
figure;
imagesc(scaledCorrelationMatrix);
colorbar;
title('Scaled Correlation Matrix (Positive Values Only, Scaled 0-1)');
xlabel('Boutons');
ylabel('Boutons');

% Calculate the average positive correlation for each bouton
meanPositiveCorrelationsPerBouton = zeros(numColumns, 1);

for i = 1:numColumns
    % Exclude diagonal element and compute mean of positive correlations for each bouton
    positiveValues = scaledCorrelationMatrix(i, :);
    positiveValues(i) = []; % Remove the self-correlation (diagonal)
    meanPositiveCorrelationsPerBouton(i) = mean(positiveValues(positiveValues > 0));
end

% Display results
disp('Average Positive Correlation per Bouton:');
disp(meanPositiveCorrelationsPerBouton);

% Create a table with each bouton’s average positive correlation
BoutonResults = table((1:numColumns)', meanPositiveCorrelationsPerBouton, ...
    'VariableNames', {'Bouton', 'AveragePositiveCorrelation'});


% Calculate mean of positive correlations
n = size(scaledCorrelationMatrix, 1);
offDiagonalElements = scaledCorrelationMatrix(~eye(n));
meanPositiveCorrelation = mean(offDiagonalElements);

disp(['Mean of positive correlation matrix: ', num2str(meanPositiveCorrelation)]);

% Create a table with the calculated values
CorrelationResults = table(meanPositiveCorrelation, ...
    'VariableNames', {'MeanPositiveCorrelation'});

% Display the table
disp('Results Table:');
disp(CorrelationResults);

% Define the new column data and label
newColumnLabel = 'File_ID';
newColumnData = [baseFileName];

% Add new column to the table with the specific label
CorrelationResults.(newColumnLabel) = newColumnData;

% Export the table to an Excel file
outputFileName = [baseFileName, '_PositiveCorrelationResults.xlsx'];
writetable(CorrelationResults, outputFileName);

disp(['Results have been exported to ', outputFileName]);

%% BOUTON LEVEL DATA

% Compile bouton averages together in one table
AVGBoutonResults = table((1:numColumns)', meanPositiveCorrelationsPerBouton, averagePks, Frequency, sumAUCs, ...
   'VariableNames', {'Bouton', 'AveragePositiveCorrelation', 'Amplitude', 'Frequency', 'AUC'});

% Remove rows where Frequency equals 0
AVGBoutonResults = AVGBoutonResults(AVGBoutonResults.Frequency ~= 0, :);

% Define file identifier and add to table
AVGBoutonResults.File_ID = repmat({baseFileName}, size(AVGBoutonResults, 1), 1);

% Export the table to an Excel file without variable names in the first row
outputFileName = [baseFileName, '_Bouton_statistics.xlsx'];
writetable(AVGBoutonResults, outputFileName, 'WriteVariableNames', false);

disp(['Bouton-wise average results have been exported to ', outputFileName]);

%clear;

    else
        disp('The imported variable is not a numeric matrix.');
    end
else
    disp('The variable was not created.');
end
