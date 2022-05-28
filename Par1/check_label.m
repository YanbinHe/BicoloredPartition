function [label_checked] = check_label(label)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if label(1) == 1
    label_checked = label;
else
    label_checked = 1*(~label);
end
end

