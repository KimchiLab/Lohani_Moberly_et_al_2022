% simple function to read individual tiff files and load as a stack into data matrix (NxMxT)
% 04.14.09 Tsai-Wen Chen (modified by Sweyta lohani 2020)
% 09.20.10 Improve reading over network by first creating local copy
function data=readTifStack_modified(varargin)
%  readTifStack(filename)
%  readTifStack(filename,index)
%  readTifStack(filename,firstim,lastim)
movelocal=0;
index=[];
if nargin ==0
  [filename, tif_path] = uigetfile('*.tif','choose a file');
  if isequal(filename,0);return;end
  filename = [tif_path filename];  
else
  tiffs=varargin{2};
  filePath=varargin{1};
end

if nargin == 3 
    index=varargin{3};
end

if nargin ==4
    index=(varargin{3}:varargin{4});
end

if nargin ==5
    index=(varargin{3}:varargin{4});
    movelocal=varargin{5};
end

if nargin==6
    index=(varargin{3}:varargin{4});
    movelocal=varargin{5};
    TmpFile=varargin{6};
end 
%%
    filename=TmpFile;
    disp('create local copy');

%Initialize data
data_start_ind = 0;

%Get base parameters from first image
filename=fullfile(filePath,tiffs(1).name);
info=imfinfo(filename);

data = zeros(info(1).Height,info(1).Width,length(info)*size(tiffs, 1),'single');
%Go through list of tiff files
for file_num=1:size(tiffs, 1)
    filename=fullfile(filePath,tiffs(file_num).name);

    if movelocal %Copy locally if desired
        [pathstr, name]=fileparts(filename);
        if isempty(pathstr)
            filename=[pwd,'/',filename];
        end
        copyfile(filename, TmpFile,'f');
        filename=fullfile(TmpFile,tiffs(file_num).name);
    end

    %Get filename info to set length param
    info = imfinfo(filename);

    %Read all images of the stack to data
    for frame = 1:length(info)
        data(:,:, data_start_ind+frame) ...
            = imread(filename, "Index", frame);
    end

    %Iterate start point for next tif file
    data_start_ind = data_start_ind+length(info);
end
