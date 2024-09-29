function [weightedPeriodMean, weightedLagTimeMean, RsqlargestRegionMask, linearPeriods, ...
    withinConeOfInfluence,weightedPeriodMean_iclucoi,weightedLagTimeMean_iclucoi,...
    weightedinterpolatedAngleWxyMean_iclucoi,weightedinterpolatedAngleWxyMean,...
    weightedPeriodMeanxwt, weightedLagTimeMeanxwt, RsqlargestRegionMaskxwt, ...
    weightedPeriodMean_iclucoixwt,weightedLagTimeMean_iclucoixwt,...
    weightedinterpolatedAngleWxyMean_iclucoixwt,weightedinterpolatedAngleWxyMeanxwt] = coefficientweightingWTC(x, y, varargin)
% This function computes the weighted mean lag time and period using wavelet coherence (WTC) and wavelet cross-spectrum (XWT).
% 这个函数使用小波相干（WTC）和小波互谱（XWT）计算加权平均滞后时间和周期。

%范围设定为执行周期范围内

% Default configuration values
% 默认配置值
Dj = 1/100;
periodRangeNum = [8 16];

% Parse optional inputs
% 解析可选输入
if ~isempty(varargin)
    for i = 1:2:length(varargin)
        switch varargin{i}
            case 'Dj'
                Dj = varargin{i+1};
            case 'PeriodRange'
                periodRangeNum = varargin{i+1};
        end
    end
end

% Compute XWT and WTC
% 计算XWT和WTC
[waveletCrossSpectrum, ~, ~, ~, sig95xwt] = xwt(x, y, 'Dj', Dj);
[waveletCoherence, period, ~, coi, sig95] = wtc(x, y, 'Dj', Dj);

% Compute phase angle in rdian of the cross spectrum
% 计算互谱的弧度格式相角
angleWxy = angle(waveletCrossSpectrum);

% Define a new linear period array
% 定义一个新的线性周期数组
linearPeriods = 2:0.1:round(max(period));

% Define the time array
% 定义时间数组
timeArray = 1:length(x);

% Create grids for the original period and time
% 为原始周期和时间创建网格
[originalTimeGrid, originalPeriodGrid] = meshgrid(timeArray, period);
% Create grids for the linear period and time
% 为线性周期和时间创建网格
[linearTimeGrid, linearPeriodGrid] = meshgrid(timeArray, linearPeriods);

% Interpolate angleWxy onto the new linear period array
% 插值angleWxy到新的线性周期数组上
% The original wavelet analysis divides the period as an exponential rise
% with a base of 2, so the weighting process leads to an overdensity of low
% periods. Therefore, the periods are interpolated to be linear to avoid
% this problem.
% 原始小波分析划分的周期为以2为底的指数上升，因此加权过程会导致低周期过密。因此，将周期插值为线性以避免该问题。
interpolatedAngleWxy = interp2(originalTimeGrid, originalPeriodGrid, angleWxy, linearTimeGrid, linearPeriodGrid, 'makima');
% interpolatedAngleWxy(interpolatedAngleWxy < 0) = interpolatedAngleWxy(interpolatedAngleWxy < 0) + (2 * pi);

% Interpolate significantCoherence and
% 插值显著谱和小波相干系数谱
interpolatedSig95 = interp2(originalTimeGrid, originalPeriodGrid, sig95, linearTimeGrid, linearPeriodGrid, 'makima');
interpolatedSig95xwt = interp2(originalTimeGrid, originalPeriodGrid, sig95xwt, linearTimeGrid, linearPeriodGrid, 'makima');

interpolatedWaveletCoherence = interp2(originalTimeGrid, originalPeriodGrid, waveletCoherence, linearTimeGrid, linearPeriodGrid, 'makima');

% Define masks to identify valid values within the region of interest
% 定义掩码以识别感兴趣区域内的有效值
withinConeOfInfluence = (linearPeriods(:) * (1 ./ coi) > 1);
significantCoherence = (interpolatedSig95 >= 1);
significantCoherencexwt = (interpolatedSig95xwt >= 1);

% Specified period range
% 指定的周期范围
periodRange = (linearPeriods > periodRangeNum(1)) & (linearPeriods < periodRangeNum(2));

% Get row indices corresponding to the specified period range
% 获取与指定周期范围对应的行索引
rowsInRange = find(periodRange);

% Label connected regions in the original significantCoherence matrix
% 在原始显著相干性矩阵中标记连接区域
[labeledRegions, numRegions] = bwlabel(significantCoherence, 8);
[labeledRegionsxwt, numRegionsxwt] = bwlabel(significantCoherencexwt, 8);

% Get properties of each region, including area and pixel indices
% 获取每个区域的属性，包括面积和像素索引
regionProperties = regionprops(labeledRegions, 'Area', 'PixelIdxList');
regionPropertiesxwt = regionprops(labeledRegionsxwt, 'Area', 'PixelIdxList');

% Find the area of each connected region within the specified period range
% 查找指定周期范围内每个连接区域的面积
maxArea = 0;
maxRegionIdx = 0;
for i = 1:numRegions
    pixelIdxList = regionProperties(i).PixelIdxList;
    pixelIdxList = pixelIdxList(~withinConeOfInfluence(pixelIdxList));
    [rowIdx, ~] = ind2sub(size(significantCoherence), pixelIdxList);
    inRange = ismember(rowIdx, rowsInRange);
    if any(inRange)
        currentArea = sum(inRange);
        if currentArea > maxArea
            maxArea = currentArea;
            maxRegionIdx = i;
        end
    end
end
% Create a logical matrix of the same size as the original significantCoherence
% 创建与原始显著相干性矩阵大小相同的逻辑矩阵
largestRegionMask = false(size(significantCoherence));
if maxRegionIdx > 0
    largestRegionMask(regionProperties(maxRegionIdx).PixelIdxList) = true;
end
% largestRegionMask(~ismember(1:size(significantCoherence,1), rowsInRange),
% :) = false; % inrange verison worked

% Find the XWT area of each connected region within the specified period range
% 查找XWT指定周期范围内每个连接区域的面积
maxAreaxwt = 0;
maxRegionIdxxwt = 0;
for i = 1:numRegionsxwt
    pixelIdxListxwt = regionPropertiesxwt(i).PixelIdxList;
    pixelIdxListxwt = pixelIdxListxwt(~withinConeOfInfluence(pixelIdxListxwt));
    [rowIdxxwt, ~] = ind2sub(size(significantCoherencexwt), pixelIdxListxwt);
    inRangexwt = ismember(rowIdxxwt, rowsInRange);
    if any(inRangexwt)
        currentAreaxwt = sum(inRangexwt);
        if currentAreaxwt > maxAreaxwt
            maxAreaxwt = currentAreaxwt;
            maxRegionIdxxwt = i;
        end
    end
end

% Create a logical matrix of the same size as the original significantCoherence
% 创建与原始显著相干性矩阵大小相同的逻辑矩阵
largestRegionMaskxwt = false(size(significantCoherencexwt));
if maxRegionIdxxwt > 0
    largestRegionMaskxwt(regionPropertiesxwt(maxRegionIdxxwt).PixelIdxList) = true;
end
%largestRegionMaskxwt(~ismember(1:size(significantCoherencexwt,1),
%rowsInRange), :) = false; % inrange verison worked

% Compute the weighted average lag time and period using wavelet coherence
% coefficients as weights %including coi%
% 使用小波相干系数作为权重计算包括 %coi范围内% 的加权平均滞后时间和周期

% weight 权重
weightedCoherenceValues_iclucoi = (interpolatedWaveletCoherence .* largestRegionMask) ./ ...
    sum((interpolatedWaveletCoherence .* largestRegionMask), 1);
weightedCoherenceValues_iclucoi(weightedCoherenceValues_iclucoi == 0) = nan;

weightedCoherenceValues_iclucoixwt = (interpolatedWaveletCoherence .* largestRegionMaskxwt) ./ ...
    sum((interpolatedWaveletCoherence .* largestRegionMaskxwt), 1);
weightedCoherenceValues_iclucoixwt(weightedCoherenceValues_iclucoixwt == 0) = nan;

% period 周期
weightedPeriodMean_iclucoi = sum(linearPeriodGrid .* weightedCoherenceValues_iclucoi, 1, 'omitnan');
weightedPeriodMean_iclucoi(weightedPeriodMean_iclucoi == 0) = nan;

weightedPeriodMean_iclucoixwt = sum(linearPeriodGrid .* weightedCoherenceValues_iclucoixwt, 1, 'omitnan');
weightedPeriodMean_iclucoixwt(weightedPeriodMean_iclucoixwt == 0) = nan;

% 重塑为与withinConeOfInfluence相同的大小
% Reshape to the same size as withinConeOfInfluence
weightedPeriodMean_matrix = zeros(size(withinConeOfInfluence));
weightedPeriodMean_matrixxwt = zeros(size(withinConeOfInfluence));


for t = 1:length(weightedPeriodMean_iclucoi)
    if ~isnan(weightedPeriodMean_iclucoi(t))
        % 找到与linearPeriods最接近的位置
        % Find the closest position in linearPeriods
        [~, periodIdx] = min(abs(linearPeriods - weightedPeriodMean_iclucoi(t)));
        weightedPeriodMean_matrix(periodIdx, t) = weightedPeriodMean_iclucoi(t);
    end
end
weightedPeriodMean = sum(weightedPeriodMean_matrix.*~withinConeOfInfluence,1);
weightedPeriodMean(weightedPeriodMean == 0) = nan;
weightedPeriodMean_matrix_logical = logical(weightedPeriodMean_matrix);

for t = 1:length(weightedPeriodMean_iclucoixwt)
    if ~isnan(weightedPeriodMean_iclucoixwt(t))
        % 找到与linearPeriods最接近的位置
        % Find the closest position in linearPeriods
        [~, periodIdxxwt] = min(abs(linearPeriods - weightedPeriodMean_iclucoixwt(t)));
        weightedPeriodMean_matrixxwt(periodIdxxwt, t) = weightedPeriodMean_iclucoixwt(t);
    end
end

weightedPeriodMeanxwt = sum(weightedPeriodMean_matrixxwt.*~withinConeOfInfluence,1);
weightedPeriodMeanxwt(weightedPeriodMeanxwt == 0) = nan;
weightedPeriodMean_matrix_logicalxwt = logical(weightedPeriodMean_matrixxwt);


%Timelag 时滞
% Weighted average radian phase angle, circular format is computed nonlinearly
% 加权平均弧度相位角，圆格式的计算是非线性的，检测是否存在负值弧度
% weightedinterpolatedAngleWxyMean_iclucoi = angle(sum(weightedCoherenceValues_iclucoi.* exp(1i * interpolatedAngleWxy),1,"omitmissing"));
weightedinterpolatedAngleWxyMean_iclucoi = sum(interpolatedAngleWxy.* weightedPeriodMean_matrix_logical);
weightedinterpolatedAngleWxyMean_iclucoi(isnan(weightedPeriodMean_iclucoi)) = nan;

weightedinterpolatedAngleWxyMean_iclucoixwt = sum(interpolatedAngleWxy.* weightedPeriodMean_matrix_logicalxwt);
weightedinterpolatedAngleWxyMean_iclucoixwt(isnan(weightedPeriodMean_iclucoixwt)) = nan;

% 如果存在加权平均弧度相位角
if ~isempty(weightedinterpolatedAngleWxyMean_iclucoi)
    % 计算弧度相位角的差值
    weightedinterpolatedAngleWxyMean_iclucoi_diff = diff(weightedinterpolatedAngleWxyMean_iclucoi);
    
    % 处理相位角差值
    weightedinterpolatedAngleWxyMean_iclucoi_diff(weightedinterpolatedAngleWxyMean_iclucoi_diff > pi) = weightedinterpolatedAngleWxyMean_iclucoi_diff(weightedinterpolatedAngleWxyMean_iclucoi_diff > pi) - 2 * pi;
    weightedinterpolatedAngleWxyMean_iclucoi_diff(weightedinterpolatedAngleWxyMean_iclucoi_diff < -pi) = weightedinterpolatedAngleWxyMean_iclucoi_diff(weightedinterpolatedAngleWxyMean_iclucoi_diff < -pi) + 2 * pi;
    
    % 找到 weightedinterpolatedAngleWxyMean_iclucoi 第一个非 NaN 的值
    first_non_nan_index = find(~isnan(weightedinterpolatedAngleWxyMean_iclucoi), 1, 'first');
    first_non_nan_value = weightedinterpolatedAngleWxyMean_iclucoi(first_non_nan_index);
    
    % 找到 weightedinterpolatedAngleWxyMean_iclucoi_diff 第一个非 NaN 的位置
    first_non_nan_diff_index = find(~isnan(weightedinterpolatedAngleWxyMean_iclucoi_diff), 1, 'first');
    
    % 插入第一个非 NaN 的值到 diff 数组中
    if ~isempty(first_non_nan_diff_index)
        weightedinterpolatedAngleWxyMean_iclucoi_diff = [weightedinterpolatedAngleWxyMean_iclucoi_diff(1:first_non_nan_diff_index-1), first_non_nan_value, weightedinterpolatedAngleWxyMean_iclucoi_diff(first_non_nan_diff_index:end)];
    else
        weightedinterpolatedAngleWxyMean_iclucoi_diff = first_non_nan_value;
    end
    
    % 计算累加和
    weightedinterpolatedAngleWxyMean_iclucoi = cumsum(weightedinterpolatedAngleWxyMean_iclucoi_diff,"omitnan");
    weightedinterpolatedAngleWxyMean_iclucoi(isnan(weightedPeriodMean_iclucoi)) = nan;
end
weightedinterpolatedAngleWxyMean = weightedinterpolatedAngleWxyMean_iclucoi;
weightedinterpolatedAngleWxyMean(isnan(weightedPeriodMean)) = nan;

if ~isempty(weightedinterpolatedAngleWxyMean_iclucoixwt)
    % 计算弧度相位角的差值
    weightedinterpolatedAngleWxyMean_iclucoi_diffxwt = diff(weightedinterpolatedAngleWxyMean_iclucoixwt);
    
    % 处理相位角差值
    weightedinterpolatedAngleWxyMean_iclucoi_diffxwt(weightedinterpolatedAngleWxyMean_iclucoi_diffxwt > pi) = weightedinterpolatedAngleWxyMean_iclucoi_diffxwt(weightedinterpolatedAngleWxyMean_iclucoi_diffxwt > pi) - 2 * pi;
    weightedinterpolatedAngleWxyMean_iclucoi_diffxwt(weightedinterpolatedAngleWxyMean_iclucoi_diffxwt < -pi) = weightedinterpolatedAngleWxyMean_iclucoi_diffxwt(weightedinterpolatedAngleWxyMean_iclucoi_diffxwt < -pi) + 2 * pi;
    
    % 找到 weightedinterpolatedAngleWxyMean_iclucoi 第一个非 NaN 的值
    first_non_nan_indexxwt = find(~isnan(weightedinterpolatedAngleWxyMean_iclucoixwt), 1, 'first');
    first_non_nan_valuexwt = weightedinterpolatedAngleWxyMean_iclucoixwt(first_non_nan_indexxwt);
    
    % 找到 weightedinterpolatedAngleWxyMean_iclucoi_diff 第一个非 NaN 的位置
    first_non_nan_diff_indexxwt = find(~isnan(weightedinterpolatedAngleWxyMean_iclucoi_diffxwt), 1, 'first');
    
    % 插入第一个非 NaN 的值到 diff 数组中
    if ~isempty(first_non_nan_diff_indexxwt)
        weightedinterpolatedAngleWxyMean_iclucoi_diffxwt = [weightedinterpolatedAngleWxyMean_iclucoi_diffxwt(1:first_non_nan_diff_indexxwt-1), first_non_nan_valuexwt, weightedinterpolatedAngleWxyMean_iclucoi_diffxwt(first_non_nan_diff_indexxwt:end)];
    else
        weightedinterpolatedAngleWxyMean_iclucoi_diffxwt = first_non_nan_valuexwt;
    end
    
    % 计算累加和
    weightedinterpolatedAngleWxyMean_iclucoixwt = cumsum(weightedinterpolatedAngleWxyMean_iclucoi_diffxwt,"omitnan");
    weightedinterpolatedAngleWxyMean_iclucoixwt(isnan(weightedPeriodMean_iclucoixwt)) = nan;
end
weightedinterpolatedAngleWxyMeanxwt = weightedinterpolatedAngleWxyMean_iclucoixwt;
weightedinterpolatedAngleWxyMeanxwt(isnan(weightedPeriodMeanxwt)) = nan;

% Dealing with radians, such as the existence of a negative radian indicating advancement [-pi 0] 
% (wavelet direction facing up) is considered as the result of increasing from 0 to that position, i.e. [pi 2pi]

% % 去除包含 NaN 的部分计算标准差和均值
% validIndices = ~isnan(weightedinterpolatedAngleWxyMean);
% validAngles = weightedinterpolatedAngleWxyMean(validIndices);
% % [~,~,sigma] = anglemean(validAngles);
% 
% countn = 0;%记录是否平移
% %如果均值-标准差为负弧度,则认为存在显著负弧度，表示提前的负值弧度[-pi 0](小波方向朝上）视为从0开始顺时针增加到该位置的结果，即[pi
% %2pi],考虑到，负值较小时，如果仅处理负值+2pi，视为转到该位置，则会导致较小正值-较小负值之间存在接近2pi的差值，该情况不符合实际，因此
% %应该为整个时滞谱滞后 + 2pi，将相位角平移1个完整周期（2pi）以确保它们在适当的范围内。
% if any(~isnan(validAngles))
%     meantheta = mean(validAngles);
%     confmean95 = circ_confmean(validAngles');
%     if isnan(confmean95)
%         confmean95 = 0;
%     end
%     while (meantheta + (2*pi *countn) - confmean95 < 0) && countn == 0
%         weightedinterpolatedAngleWxyMean_iclucoi = weightedinterpolatedAngleWxyMean_iclucoi + 2 * pi;
%         weightedinterpolatedAngleWxyMean = weightedinterpolatedAngleWxyMean + 2 * pi;
%         % 去除包含 NaN 的部分计算标准差和均值
%         countn = 1;
%     end
% end

% % Calculate the lag time using the interpolated angleWxy
% % 使用插值后的angleWxy计算滞后时间
% lagTimeXWT = interpolatedAngleWxy ./ (2 * pi) .* linearPeriodGrid;
% % lagTimeXWT(lagTimeXWT == 0) = nan;
% lagTimeXWT_iclucoi = lagTimeXWT .* largestRegionMask;
% lagTimeXWT_iclucoi(~largestRegionMask) = NaN;
weightedLagTimeMean_iclucoi = weightedinterpolatedAngleWxyMean_iclucoi / (2 * pi) .* weightedPeriodMean_iclucoi;
weightedLagTimeMean_iclucoi(isnan(weightedPeriodMean_iclucoi)) = nan;
weightedLagTimeMean = weightedinterpolatedAngleWxyMean / (2 * pi) .* weightedPeriodMean;
weightedLagTimeMean(isnan(weightedPeriodMean)) = nan;

weightedLagTimeMean_iclucoixwt = weightedinterpolatedAngleWxyMean_iclucoixwt / (2 * pi) .* weightedPeriodMean_iclucoixwt;
weightedLagTimeMean_iclucoixwt(isnan(weightedPeriodMean_iclucoixwt)) = nan;
weightedLagTimeMeanxwt = weightedinterpolatedAngleWxyMeanxwt / (2 * pi) .* weightedPeriodMeanxwt;
weightedLagTimeMeanxwt(isnan(weightedPeriodMeanxwt)) = nan;

% Find the waveletCoherence of the largest significant region
% 找到最大显著区域的小波相干系数
RsqlargestRegionMask = largestRegionMask .*interpolatedWaveletCoherence;
RsqlargestRegionMask(RsqlargestRegionMask  == 0) = NaN;

RsqlargestRegionMaskxwt = largestRegionMaskxwt .*interpolatedWaveletCoherence;
RsqlargestRegionMaskxwt(RsqlargestRegionMaskxwt  == 0) = NaN;

    
    
end
