function [amp_m] = FtrcToMat(part_s, maxLength, col)
%function [amp_m] = FtrcToMat(part_s, maxLength, col)

    
    
if nargin < 2 | isempty(maxLength)
    maxLength = [];
end

if nargin < 3
    col = 3;
end

fieldToFind_c = {'HRM','TRC', 'MD_1HRM', 'MD_1TRC'};
fieldname_c = fieldnames(part_s);
field = intersect(fieldname_c, fieldToFind_c);
field = field{1};

if isempty(maxLength)
    maxLength = 0;
    for k = 1:length(part_s)
        if ~isempty(part_s(k).(field)(:,1))
            maxLength = max(maxLength, max(part_s(k).(field)(:,1)));
        end
    end
end
amp_m = zeros(length(part_s), maxLength);

for k = 1:length(part_s)
    if ~isempty(part_s(k).(field))
        maxIndice = min(maxLength, max(part_s(k).(field)(:,1)));
        pos_v = find(part_s(k).(field)(:,1) <= maxIndice);
        [sortval_v,sortpos_v] = sort(part_s(k).(field)(pos_v,1));
        pos_v = pos_v(sortpos_v);
        amp_m(k, part_s(k).(field)(pos_v,1)) = part_s(k).(field)(pos_v,col);
    end
end
