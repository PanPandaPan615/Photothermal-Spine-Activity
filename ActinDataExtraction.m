% This script takes multi image tif stacks containing segmentation of
% neurons and uses it ot measure the mean insity of each structure in the 
% neuron. Each frame in the tiff file contains either the soma, dendrite,
% or dendritic spines. 

clc, clear, close all

cd('\\10.4.151.15\mahadevan-jansenlab_data_storage\Projects\INS\Cells\Neurons\Microscopy\Confocal\Dendritic Spine Study\Actin_GCaMP');
addpath(genpath('\\10.4.151.15\Mahadevan-JansenLab_Data_Storage\Projects\INS\Cells\Neurons\Code\MATLAB\Functions'))

% Experimental parameters

tic

% Get list of directory contents
dircontents = dir(fullfile(cd, '**', '*.*')); 

% Experimental Folders
folderlist = dircontents([dircontents.isdir]);
tempmask = contains({folderlist.name},'Actin')';
folderlist = folderlist(tempmask);
folderlist = fullfile({folderlist.folder},{folderlist.name});

% Experimental Files
filelist = dircontents(~[dircontents.isdir]);  
filelist = fullfile({filelist.folder},{filelist.name});
[filepath,filename,ext] = fileparts(filelist); % get file parts

% Location of nd2 files of data to be measured in filelist
tempmask = contains(filepath,"Raw Data"); ndloc = strcmp(ext,".nd2"); 
ndloc = ndloc & tempmask;

% location of tif files in filelist
tempmask = contains(filepath,"Segmented"); tiffloc = strcmp(ext,".tif"); 
segloc = tiffloc & tempmask;

tempmask = contains(filepath,"GCaMP_MIP");
MIPloc = tiffloc & tempmask;

% initialize structure
ExpData = [];
[optimizer, metric] = imregconfig("multimodal");
optimizer.MaximumIterations = 1e3;
optimizer.Epsilon = 1.5e-4;

% For each tif file
%i = find(tiffloc);
for i = find(segloc)
   % Read tif file image
   Segmented = readTiffStack(filelist{i});

   % Location of all nd2 files with same dish and neuron number as tiff file
   nameloc = contains(filename,filename{i}) & ndloc; 
   
   % Location of tiff file in folderlist 
   folderloc = RecursiveContains(filelist{i}, folderlist);
   
   % Location of ndfiles in the folder with the same name as the tif files
   filelocs = contains(filelist, folderlist{folderloc}) & nameloc; 
   GCaMP_MIP = contains(filename,filename{i}) & MIPloc & contains(filelist, folderlist{folderloc});

   GCaMP_MIP = imread(filelist{GCaMP_MIP});
   GCaMP_MIP = rescale(GCaMP_MIP);

   % read the stack for processing take measurements in this loop
    k = 1;
    for j = find(filelocs)
        name = filename{j};
        ImAct = nd2read(filelist{j});
        ImAct = imgaussfilt(ImAct,1);
        ndinfo  = nd2info(filelist{j});

        % Registration      
        tform = imregtform(GCaMP_MIP,rescale(ImAct),'affine',optimizer,...
                           metric,'InitialTransformation',affinetform2d);

        tformInv = invert(tform);
        ImAct = imwarp(ImAct, tformInv, 'OutputView',  imref2d(size(ImAct)));
        
        Composite = cat(3,rescale(ImAct), GCaMP_MIP, zeros(size(ImAct)));

        B = max(Segmented,[],3); B = edge(B);
        Composite = imoverlay(Composite,B);
        imshow(Composite)
        
        % Measure Neuron Intensity
        NeuronMeasurements = MeasureNeuronAnatomy(Segmented,ImAct);

       
        Spines(:,k) = NeuronMeasurements.Spines;
        k = k+1;
    end

    NeuronMeasurements.Spines = Spines;
    clear Soma Dendrites Spines

    PathParts = strsplit(filelist{j},'\');
    Expdate = find(contains(PathParts,"2025"),1,'last');
    Expdate = PathParts{Expdate};

    NeuronMeasurements.filename = fullfile(Expdate,filename{j});
    NeuronMeasurements.ndinfo = ndinfo;
    NeuronMeasurements.Segmentation = Segmented(:,:,2:end);
    
    ExpData = [ExpData; NeuronMeasurements];
   
end

save(cd,"ExpData",'-v7.3');

toc

%% Functions
function NeuronMeasurements = MeasureNeuronAnatomy(MaskStack,ndimage)
   sz = size(MaskStack,3);
   x = (sz-2)/2;
   
   % Split into anatomical stacks
   Spinesroi = MaskStack(:,:,3+x:end);
    
   SpineParent = []; Spines = [];
    % Measure Spines Intensity
    for i = 1:size(Spinesroi,3)
        mask = Spinesroi(:,:,i);
        raw = regionprops(mask,ndimage,"MeanIntensity");
        raw = [raw.MeanIntensity];

        SpineParent = [SpineParent, i*ones(1,size(raw,2))];
        Spines = horzcat(Spines, raw);
    end
    
    NeuronMeasurements.Spines = Spines;
    NeuronMeasurements.SpineParent = SpineParent;
end
