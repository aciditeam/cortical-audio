function [ features ] = processDescriptor(inFile)
    if ~exist(inFile, 'file')
        error([inFile ' does not exist']);
    end
    system(['./ircamdescriptor ' inFile ' ircamDescriptor.cfg -o/tmp/descriptors.sdif']);
    [fHandle, head] = Fsdifopen('/tmp/descriptors.sdif', 'r');
    infos = struct;
    features = struct;
    frametypes = struct;
    disp('* Reading SDIF descriptors headers');
	disp(head);
    for j = 1:length(head.TYP.FTD(1).msig)
        if head.TYP.FTD(1).msig(j, 1) ~= 'I'
            signal = ['MD_' char(head.TYP.FTD(1).msig(j, :))];
            infos.(signal) = struct;
            infos.(signal).name = char(head.TYP.FTD(1).mname{j});
            featName = char(head.TYP.FTD(1).mname{j});
            features.(featName) = struct;
            features.(featName).name = featName;
            features.(featName).value = [];
            features.(featName).times = [];
        else
            signal = char(head.TYP.FTD(1).msig(j, :));
            infos.(signal) = struct;
            infos.(signal).name = char(head.TYP.FTD(1).mname{j});
            tempName = char(head.TYP.FTD(1).mname{j});
            frameType = char(head.TYP.FTD(1).msig(j, 2:end));
            frametypes.(frameType) = struct;
            frametypes.(frameType).name = tempName(1:(end - 4));
        end
    end
    frames = Fsdifread(fHandle);
    disp('* Parsing through analysis frames');
    while ~isempty(frames)
        frameNames = fieldnames(frames.data);
        if strcmp(char(frames.fsig), '1DSC')
            for k = 1:length(frameNames)
                if isfield(infos, frameNames{k})
                    finalName = infos.(frameNames{k}).name;
                    if size(frames.data.(frameNames{k}), 1) > 1
                        features.(finalName).value = [features.(finalName).value frames.data.(frameNames{k})];
                    else
                        features.(finalName).value = [features.(finalName).value ; frames.data.(frameNames{k})];
                    end
                    features.(finalName).times = [features.(finalName).times frames.time];
                end
            end
        end
		frames = Fsdifread(fHandle);
    end
    featNames = fieldnames(features);
	Fsdifclose(fHandle);
    disp('* Modifying descriptors structure');
    for f = 1:length(featNames)
        switch size(features.(featNames{f}).value, 2)
            case 9
                features.(featNames{f}).value = features.(featNames{f}).value(:, 7:9);
            case 6
                features.(featNames{f}).value = features.(featNames{f}).value(:, 6);
            case 3
                features.(featNames{f}).value = features.(featNames{f}).value(:, 3);
            otherwise
        end
    end
    if isempty(features.EnergyEnvelope.value)
        features.EnergyEnvelope.value = features.TotalEnergy.value;
    end
end

