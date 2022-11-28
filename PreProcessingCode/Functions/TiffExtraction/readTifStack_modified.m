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
% local=['C','D','E','F','G'];
%% get the name of the first file in the tiff sequence
filename=fullfile(filePath,tiffs(1).name);
if movelocal
    [pathstr, name]=fileparts(filename);
    if isempty(pathstr)
        filename=[pwd,'/',filename];
    end
    copyfile(filename, TmpFile,'f');
    filename=TmpFile;
    disp('create local copy');
end

info=imfinfo(filename);
nimage=length(info);
if isempty(index)
    index=1:nimage;
end
%%
nread=length(find(index<=size(tiffs,1)));
%data=zeros(info(1).Height,info(2).Width,nread);
data=zeros(info(1).Height,info(1).Width,nread,'single');
for i=1:length(index)
    if index(i)<=size(tiffs,1)
        filename=fullfile(filePath,tiffs(index(i)).name);
        if movelocal
            [pathstr, name]=fileparts(filename);
            if isempty(pathstr)
                filename=[pwd,'/',filename];
            end
            copyfile(filename, TmpFile,'f');
            filename=TmpFile;
        end
        data(:,:,i)=imread(filename);
    else
        break;
    end
end
