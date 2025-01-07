%% INTER Vs INTRA bouton correlation 

% Specify the folder containing the CSV files
folderPath = 'C:\Users\ms0405\Suite2P\Traces'; % Update with your folder path

% Get a list of all CSV files in the folder
fileList = dir(fullfile(folderPath, '*.csv'));

% Loop through each file in the folder
for fileIdx = 1:length(fileList)
    % Get the current filename
    filename = fullfile(folderPath, fileList(fileIdx).name);

    % Extract the base filename (without path and extension) - Helps track files
    [~, baseFileName, ~] = fileparts(filename);

    % Create a valid MATLAB variable name
    validVarName = matlab.lang.makeValidName(baseFileName);

    % Check if the file exists
    if exist(filename, 'file') == 2
        try
            % STEP 1: IMPORT DATA
            % Read the full data table
            dataTable = readtable(filename);

            % Extract specific columns (e.g., B, D, F, H, J; indices 2, 4, 6, 8, 10)
            validColumns = [2, 4, 6, 8, 10];
            validColumns = validColumns(validColumns <= width(dataTable)); % Ensure within bounds
            selectedColumns = dataTable(:, validColumns);

            % Convert table to array
            dataArray = table2array(selectedColumns);

            % Normalize data
            columnMeans = mean(dataArray, 'omitnan');
            normalizedArray = dataArray ./ columnMeans;
            normalizedDataTable = array2table(normalizedArray, 'VariableNames', selectedColumns.Properties.VariableNames);

            % STEP 2: DRIFT CORRECTION
            windowSize = 25;  % Window size for moving average
            correctedDataArray = zeros(size(dataArray));

            for col = 1:size(dataArray, 2)
                smoothedData = movmean(dataArray(:, col), windowSize);
                correctedDataArray(:, col) = dataArray(:, col) - smoothedData;
            end
            correctedDataTable = array2table(correctedDataArray, 'VariableNames', selectedColumns.Properties.VariableNames);

            % STEP 3: PEAK DETECTION
            correctedData = table2array(correctedDataTable);
            numColumns = size(correctedData, 2);
            allPks = cell(numColumns, 1);
            allLocs = cell(numColumns, 1);
            allStdDevs = zeros(numColumns, 1);

            % Compute population standard deviation for each column
            for col = 1:numColumns
                allStdDevs(col) = std(correctedData(:, col), 'omitnan');
            end
            popStdDev = mean(allStdDevs);

            % Peak detection
            for col = 1:numColumns
                columnData = correctedData(:, col);
                minPeakHeight = 2.5 * allStdDevs(col);
                [pks, locs] = findpeaks(columnData, 'MinPeakProminence', 0.0, 'MinPeakWidth', 1.0, 'MinPeakHeight', minPeakHeight);
                allPks{col} = pks;
                allLocs{col} = locs;
            end

            % STEP 4: CALCULATE TRANSIENT FREQUENCY
            Time = 232.5; % Total duration (e.g., 900 frames at 0.388Hz)
            numColumns = length(allLocs);
            countLocs = NaN(numColumns, 1);
            Frequency = NaN(numColumns, 1);

            for col = 1:numColumns
                locs = allLocs{col};
                if ~isempty(locs)
                    countLocs(col) = numel(locs);
                    Frequency(col) = countLocs(col) / Time;
                else
                    countLocs(col) = 0;
                    Frequency(col) = 0;
                end
            end

            % STEP 5: CALCULATE TRANSIENT AMPLITUDE
            averagePks = NaN(numColumns, 1);
            stdPks = NaN(numColumns, 1);
            for col = 1:numColumns
                pks = allPks{col};
                if ~isempty(pks)
                    averagePks(col) = mean(pks);
                    stdPks(col) = std(pks);
                else
                    averagePks(col) = NaN;
                    stdPks(col) = NaN;
                end
            end

            % STEP 6: SYNCHRONICITY AND CORRELATIONS
            maxFrames = size(correctedData, 1);
            eventMatrix = zeros(maxFrames, numColumns);

            for col = 1:numColumns
                if ~isempty(allLocs{col})
                    eventMatrix(allLocs{col}, col) = 1;
                end
            end
            correlationMatrix = corr(eventMatrix);
            correlationMatrix(isnan(correlationMatrix)) = 0;

            % Normalize correlation matrix
            positiveCorrelations = correlationMatrix;
            maxValue = max(positiveCorrelations(:));
            minValue = min(positiveCorrelations(:));
            if maxValue ~= minValue
                scaledCorrelationMatrix = (positiveCorrelations - minValue) / (maxValue - minValue);
            else
                scaledCorrelationMatrix = positiveCorrelations;
            end

            % Exclude diagonal correlations (set them to NaN)
            scaledCorrelationMatrixNoDiagonal = scaledCorrelationMatrix;
            for i = 1:size(scaledCorrelationMatrix, 1)
                scaledCorrelationMatrixNoDiagonal(i, i) = 1;
            end

            % Add File_ID column to the correlation matrix table
            numRows = size(scaledCorrelationMatrixNoDiagonal, 1);
            fileIDColumn = repmat({baseFileName}, numRows, 1);  % Add the File_ID column
            CorrelationMatrixTable = array2table(scaledCorrelationMatrixNoDiagonal, ...
                'VariableNames', selectedColumns.Properties.VariableNames, ...
                'RowNames', selectedColumns.Properties.VariableNames);
            CorrelationMatrixTable.File_ID = fileIDColumn;  % Add the File_ID column

            % STEP 7: EXPORT RESULTS
            % Output Excel file
            outputFileName = fullfile(folderPath, [baseFileName, '_SequentialResults.xlsx']);
            
            % Write CorrelationMatrixTable with File_ID to Excel sheet
            writetable(CorrelationMatrixTable, outputFileName, 'Sheet', 'CorrelationMatrix');
            
            disp(['Results have been exported to ', outputFileName]);

        catch ME
            disp(['Error processing file: ' baseFileName]);
            disp(ME.message);
        end
    else
        disp(['File not found: ' filename]);
    end
end