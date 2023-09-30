function [] = mkdir_if_not_exists(savepath)
%MKDIR_IF_NOT_EXISTS Creates directory if it did not exist already
    if ~exist(savepath, 'dir')
       mkdir(savepath)
    end
end

