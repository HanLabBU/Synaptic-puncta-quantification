% written by Seth Bensussan 
function Neurogenesis_Synapse_Count()
    %Code to analyze number of synapses from neurogenesis data.
    
    %Set directory in order to access "findNestedFiles" function
    cd('C:\Users\Administrator\Documents\Matlab\Han Lab\Neurogenesis\scripts-master\matlab\utilities');
    
    %save neurogensis files into a cell called filename_list
    filename_list = findNestedFiles('U:\eng_research_handata\eng_research_handata2\Seth (handata2)\Neurogenesis', '"*MAX_Green Stack.tif"');
    
 
    filename_idx=210; %manually choose the file for synaptic analysis 
    filename = filename_list{filename_idx}; %access individual file
        [pathstr, name, ext] = fileparts(filename); %get fileparts for individual file
        cd(pathstr); %go to path of file
        image = imread('MAX_Green Stack.tif');
        
        image_filt = imgaussfilt(image,1.5); %filter/smooth file
        
        %set lowest pixels to 0
        image_bg_subtract = image_filt;
        CutoffIndex=200;
        image_bg_subtract(image_bg_subtract<CutoffIndex)=0; 
        
        %sharpen image
        imageSharp = imsharpen(image_bg_subtract,'Radius',1.5,'Amount',100.5); 
        
       BW = imageSharp; 
       BW(BW<2000)=0;
       binaryImage=zeros(size(BW));
       
       %make intensities over threshold equal 1
        CutoffIndex2=6000;
        binaryImage(BW>CutoffIndex2)=1; 
        
        %find connected components
        prelimROIS = bwconncomp(binaryImage); %find connected components (pick ROIs)
        ROIdata = regionprops(prelimROIS);
        synapseFinal = ROIdata([ROIdata.Area]>10); %filter out regions smaller than certain num pixels
        
         %overlay synapses and rois
        comp = 255*double(image)./4067;
        center = cat(1,synapseFinal.Centroid);
        ov2 = insertShape(comp,'Circle',[center 3*ones(numel(synapseFinal),1)],'Color','red');
        ov5=squeeze(ov2(:,:,1));
        figure;
        imagesc(double(image)+14000*ov5);
        
end