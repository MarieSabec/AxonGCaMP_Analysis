%% SETUP FOLDER AND FILES
% Specify the folder containing the Excel files
folderPath = 'C:\Users\ms0405\Suite2P\'; % Update with your folder path

% Get a list of all Excel files in the folder
fileList = dir(fullfile(folderPath, '*.xlsx'));

% Loop through each Excel file in the folder
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
            % Read data as a matrix from the specified range
            range = 'F2:AHU350';  % Adjust if needed
            dataMatrix = readmatrix(filename, 'Range', range);

            % Remove rows that are entirely NaN
            dataMatrix = dataMatrix(~all(isnan(dataMatrix), 2), :);

            % Transpose the matrix so that rows correspond to original columns
            dataTransposed = dataMatrix';

            % Dynamically assign the transposed matrix to the variable name
            eval([validVarName ' = dataTransposed;']);
            disp(['Data imported and transposed into variable: ' validVarName]);

            %% DRIFT CORRECTION
            % Drift correction using moving average
            windowSize = 25;  % Adjust window size as needed
            numCols = size(dataTransposed, 2);  % Number of columns
            data_corrected = zeros(size(dataTransposed));

            for col = 1:numCols
                smoothedData = movmean(dataTransposed(:, col), windowSize);
                data_corrected(:, col) = dataTransposed(:, col) - smoothedData;
            end

            % Dynamically assign the corrected data to the variable name
            eval([validVarName ' = data_corrected;']);
            disp(['Data corrected for drift and stored in variable: ' validVarName]);

            %% DETECT TRANSIENT PEAKS AND LOCATIONS
            if exist(validVarName, 'var')
                % Extract the data from the dynamically named variable
                data_corrected = eval(validVarName);

                numColumns = size(data_corrected, 2);
                allPks = cell(numColumns, 1);
                allLocs = cell(numColumns, 1);
                allStdDevs = cell(numColumns, 1);

                for col = 1:numColumns
                    columnData = data_corrected(:, col);
                    stdDev = std(columnData, 'omitnan');
                    allStdDevs{col} = stdDev;
                end

                pop_stdDev = mean(cell2mat(allStdDevs));

                for col = 1:numColumns
                    columnData = data_corrected(:, col);
                    [pks, locs] = findpeaks(columnData, 'MinPeakProminence', 1.0, 'MinPeakWidth', 1.0, 'MinPeakHeight', 2.5 * pop_stdDev);
                    allPks{col} = pks;
                    allLocs{col} = locs;
                end

                %% CALCULATE AUTOCORRELATION
                allAutocorrs = cell(numColumns, 1);
                maxLag = 100;

                for col = 1:numColumns
                    boutonData = data_corrected(:, col);
                    boutonData = boutonData(~isnan(boutonData));
                    autocorrResult = xcorr(boutonData, maxLag, 'coeff');
                    allAutocorrs{col} = autocorrResult;

                    % Plot autocorrelation for the first bouton
                    if col == 1
                        figure;
                        lags = -maxLag:maxLag;
                        plot(lags, autocorrResult);
                        title(['Autocorrelation for Bouton ' num2str(col)]);
                        xlabel('Lag (frames)');
                        ylabel('Autocorrelation');
                    end
                end

                %% DETECT RHYTHMIC FIRING AND OSCILLATORY PATTERNS USING FFT
                allPSDs = cell(numColumns, 1);
                peakFrequencies = zeros(numColumns, 1);
                peakPowers = zeros(numColumns, 1);

                for col = 1:numColumns
                    boutonData = data_corrected(:, col);
                    boutonData = boutonData(~isnan(boutonData));
                    Fs = 3.88;  % Sampling frequency
                    L = length(boutonData);
                    Y = fft(boutonData);

                    P2 = abs(Y / L);
                    P1 = P2(1:floor(L / 2) + 1);
                    P1(2:end-1) = 2 * P1(2:end-1);

                    f = Fs * (0:(L / 2)) / L;
                    allPSDs{col} = P1;

                    [maxPower, maxIdx] = max(P1);
                    peakFrequencies(col) = f(maxIdx);
                    peakPowers(col) = maxPower;

                    % Plot power spectrum for the first bouton
                    if col == 1
                        figure;
                        plot(f, P1);
                        title(['Power Spectrum for Bouton ' num2str(col)]);
                        xlabel('Frequency (Hz)');
                        ylabel('|P1(f)|');
                    end

                    disp(['Highest Peak Frequency for Bouton ' num2str(col) ': ' num2str(peakFrequencies(col)) ' Hz']);
                    disp(['Corresponding Power for Bouton ' num2str(col) ': ' num2str(peakPowers(col))]);
                end

                %% CREATE BOUTON RESULTS TABLE
                boutonNumbers = (1:numColumns)';
                fileIDs = repmat({baseFileName}, numColumns, 1);

                BoutonsresultsTable = table(fileIDs, boutonNumbers, peakFrequencies, peakPowers, ...
                    'VariableNames', {'FileID', 'BoutonNumber', 'PeakFrequency', 'PeakPower'});

                disp('Bouton Results Table:');
                disp(BoutonsresultsTable);

                % Save the table to an Excel file
                excelFileName = [baseFileName, '_BoutonPeakPower.xlsx'];
                writetable(BoutonsresultsTable, excelFileName);
                disp(['Results saved to Excel file: ' excelFileName]);

                %% AVERAGE PSD ACROSS ALL BOUTONS
                allPSDMatrix = zeros(size(allPSDs{1}, 1), numColumns);
                for col = 1:numColumns
                    allPSDMatrix(:, col) = allPSDs{col};
                end
                averagePSD = mean(allPSDMatrix, 2);

                f = Fs * (0:(L / 2)) / L;
                figure;
                plot(f, averagePSD, 'k', 'LineWidth', 1.5);
                title('Average Power Spectrum Density Across All Boutons');
                xlabel('Frequency (Hz)');
                ylabel('Average |P1(f)|');
                grid on;

                %% CREATE AND EXPORT AVERAGE PSD TABLE WITH FREQUENCY BINS
                binWidth = 0.05;
                fBins = 0:binWidth:max(f);

                binnedFrequencies = zeros(length(fBins)-1, 1);
                binnedPowers = zeros(length(fBins)-1, 1);

                for i = 1:length(fBins)-1
                    binIndices = (f >= fBins(i)) & (f < fBins(i+1));
                    binnedFrequencies(i) = mean(f(binIndices));
                    binnedPowers(i) = mean(averagePSD(binIndices), 'omitnan');
                end

                averagePSDTable = table(binnedFrequencies, binnedPowers, ...
                                        'VariableNames', {'Frequency_Hz', baseFileName});
                disp(averagePSDTable);

                % Save the average PSD table to an Excel file
                outputFileName = [baseFileName '_AveragePSD_Binned.xlsx'];
                writetable(averagePSDTable, outputFileName);
                disp(['Average PSD table saved to: ' outputFileName]);

            % catch err
            %     disp(['Error processing file ' filename ':']);
            %     disp(err.message);
        %     end
        % else
        %     disp(['File does not exist: ' filename]);
        end
    end
end

disp('All files processed.');
end
