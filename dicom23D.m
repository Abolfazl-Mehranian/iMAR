% Copyright (C) 2009 Pacific Northwest National Laboratory
%
% Authors:
%       Ariel Balter ariel.balter@pnl.gov
%       Biological Monitoring and Modeling
%       Fundamental Science and Computing Division
%       Pacific Northwest National Laboratory
%
% This software is licensed under the terms of the MIT license with
% some exceptions: The software may be freely redistributed under
% the condition that the copyright notices (including this entire
% header and the copyright notice) are not removed, and no
% compensation is received. Private, research, and institutional
% use is free. You may distribute modified versions of this code UNDER
% THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT IN THE
% SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH SOURCE
% AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR
% NOTICE IS GIVEN OF THE MODIFICATIONS. Distribution of this code as
% part of a commercial system is permissible ONLY BY DIRECT ARRANGEMENT
% WITH THE AUTHORS.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
% WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
% OTHER DEALINGS IN THE SOFTWARE.
 



function [volume_image,I, slice_data, image_meta_data] = ...
        dicom23D(dicom_directory, dicom_fields, save_format)
    
    % function [volume_image, slice_data, image_meta_data] = ...
    %         dicom23D(dicom_directory, dicom_fields, save_format)
    %
    % Usage:
    %      [V,S,I] = DICOM23D(P, F, S)  -  Reads the contents of the dicom
    %      directory P, or has user select if P=[], and outputs a 3D volume image
    %      in V.  Saves the header fields in F for each slice.  Saves
    %      geonmetric meta data abou the image, such as physical aspect ratio,
    %      in I.  All data is saved to MAT file in the dicom directory.  Future
    %      implementation will have a switch for data saving and allow
    %      alternate text based formats.
    %
    % Inputs:
    %      "dicom_directory" (optional) -  Path the the folder containing dicom
    %       files.  If none specified, user will select with dialog box.
    %
    %      "dicom_fields" (optional) -  Names of dicom header fields to store.
    %       The default values are:
    %             default_dicom_fields = {...
    %                 'Filename',...
    %                 'Rows',...
    %                 'Columns', ...
    %                 'PixelSpacing',...
    %                 'SliceThickness',...
    %                 'SliceLocation',...
    %                 'SpacingBetweenSlices'...
    %                 'ImagePositionPatient',...
    %                 'ImageOrientationPatient',...
    %                 'FrameOfReferenceUID',...
    %                 };
    %
    %      "save_format" (optional)  -  Not yet implemented.  Will save the
    %       data in a variety of text formats.  Default value is 'mat'.
    %
    % Outputs:
    %      "volume_image"  -  volume image data
    %
    %      "slice_data"  -  a vector of structs, one for each slice,
    %       containing the header data from the specified dicom fields (or
    %       the defaults) and the extra fields:
    %             extra_fields = {...
    %                 'PhysicalHeight',...  % Height (cols) of slice in mm
    %                 'PhysicalWidth',...  % Width (rows) of slice in mm
    %                 'PixelSliceLocation',... % Slice z-location in pixels
    %                 'PixelSliceThickness',... % Slice thickness in pixels
    %                 'SliceData'...  % The slice image data
    %                  }
    %
    %      "image_meta_data"  -  Struct of meta data about the image:
    %             image_meta_data = {...
    %                 'PhysicalTotalZ',... % Total extent of image in Z
    %                       direction
    %                 'NumberOfSlices',... % Number of slices
    %                 'PhysicalAspectRatio'... % Aspect ratio of image as a whole :=
    %                       [PhysicalWidth, PhysicalHeight,
    %                       SliceThickness*NumberOfSlices]
    %                 }
    
    default_dicom_fields = {...
        'Filename',...
        'Height', ...
        'Width', ...
        'Rows',...
        'Columns', ...
        'PixelSpacing',...
        'SliceThickness',...
        'SliceLocation',...
        'SpacingBetweenSlices'...
        'ImagePositionPatient',...
        'ImageOrientationPatient',...
        'FrameOfReferenceUID',...
        };
    
    % We need these checks because to calculate the "extra_fields", we
    % need to have the PixelSpacing and SliceThickness data.  If not, we
    % leave out the extra_fields.
    no_pixel_spacing = false;
    no_slice_thickness = false;
    
    extra_fields = {...
        'PhysicalHeight',...  % Height (cols) of slice in mm
        'PhysicalWidth',...  % Width (rows) of slice in mm
        'PixelSliceLocation',... % Slice z-location in pixels
        'PixelSliceThickness',... % Slice thickness in pixels
        'SliceData'...  % The slice image data
        };
    
    image_meta_data = struct(...
        'PhysicalTotalZ',[],... % Total extent of image in Z direction
        'NumberOfSlices',[],... % Number of slices
        'PhysicalAspectRatio',[]);... % Aspect ratio of image as a whole :=
        %  [PhysicalWidth, PhysicalHeight, SliceThickness*NumberOfSlices]
    
    if nargin<1
        dicom_directory = uigetdir();
    end
    
    if isempty(dicom_directory)
        dicom_directory = uigetdir();
    end
    
    if nargin<3
        save_format = 'mat';
    end
    
    if nargin<2
        dicom_fields = default_dicom_fields;
    end
    
    all_fields = [dicom_fields, extra_fields];
    
    
    %-----Pseudocode:
    %
    %     for all files in folder
    %         if file is good dicom file
    %             in an order sorted by slice z position, insert:
    %             slice_data <-- filtered header info including:
    %                 filename
    %                 pixel_z_location
    %                 physical_z_location
    %                 image slice
    %                 other requested header data
    %             image <-- image slice
    %          else report error for bad dicom file
    %     save image and image data in requested format (e.g. MAT)
    
    % (1) Loop through all files in folder
    % (2) If extension is ".dcm" then try to read header and image.  Report
    % any errors.
    % (3) Fill struct with relevant header data and image data
    % (4) Also add 'pixel_z_value', 'physical_z_value'
    % (5) Sort by 'SliceLocation'
    % (6) Deal slices to 3D image
    % (7) Save image and slice_data as MAT file
    
    warning off;
    
    % Get directory listing
    listing = dir(dicom_directory);
    % number of files
    N = numel(listing); % How many entries in the directory listing
    if (N<3)
        error('Empty folder');
        return
    end
    
    slice_data(N) = cell2struct(cell(size(all_fields)), all_fields, 2);
    
    h = waitbar(0,'Reading DICOM Files...','WindowStyle','modal');
    
    true_index = 0; % a sequential index of dicom files, that is ignoring
    % files of other types.
    
    for i = 3:length(listing) % loop through directory listing, but skip '.' and '..'
        filename = listing(i).name;
        [dummy_path, just_the_name, extension] = fileparts(filename);
        full_path = fullfile(dicom_directory, filename);
        
        goodfile = false;
        
        % Check for good dicom file
        if isdicom(full_path)
            true_index = true_index + 1;
            header = dicominfo(full_path);
            slice_image = dicomread(header);
            
            % Save selected header data into the structure slice_data
            for j = 1:numel(dicom_fields) % loop through dicom field names
                current_field = dicom_fields{j};
                % Deal with requested fields not found in header
                if isfield(header, current_field)
                    slice_data(true_index).(current_field) = header.(current_field);
                else
%                     ['header did not contain the field ' current_field];
                end %if
                
            end % loop through dicom field names
            % done saving filtered header data
            
            % Save slice data
            slice_data(true_index).SliceData = slice_image;
            
            % Save extra fields
            needed_header_tags = [...
                isfield(header, 'PixelSpacing'), ...
                isfield(header, 'SliceThickness'), ...
                isfield(header, 'SliceLocation')...
                ];
            
            if all(needed_header_tags)
                pixel_spacing = header.PixelSpacing;
                slice_data(true_index).PhysicalHeight = ...
                    double(pixel_spacing(1)*header.Columns);
                slice_data(true_index).PhysicalWidth = ...
                    double(pixel_spacing(2)*header.Rows);
                % need to double check which aspect ratio goes with cols/rows
                slice_data(true_index).PixelSliceLocation = ...
                    header.SliceLocation / mean(pixel_spacing);
                slice_data(true_index).PixelSliceThickness = ...
                    header.SliceThickness / mean(pixel_spacing);
            else
                no_pixel_spacing = true;
            end % if pixel spacing
            
        end % if isdicom
        
        waitbar(i/N,h);
    end % loop through directory listing
    % Eliminate empty structs at end.
    slice_data = slice_data(1:true_index);
        
    waitbar(1,h);
    close(h);
    warning on;
    
    % Check that some dicom slice was found
    if true_index < 1
        'No dicom slices found...returning empty';
        volume_image = [];
        slice_data = [];
        image_meta_data = [];
        return
    end
    
%     If SliceLocation is known, sort by that.  This is deemed more
%     accurate than going by filename order (or file number).
    if isfield(slice_data(1), 'SliceLocation')
        slice_data.SliceLocation
        [S,I] = sort([slice_data.SliceLocation]);
        slice_data = slice_data(I);
        
    end
    
    % Pre-allocate volume image array
    [rows, cols] = size(slice_data(1).SliceData);
    volume_image = ...
        zeros(rows, cols, length(slice_data));
    % Build volume image array
    h = waitbar(0,'Writing slice images to volume image array...','WindowStyle','modal');
    for i = 1:length(slice_data)
        waitbar(i/N,h);
        volume_image(:,:,i) = slice_data(i).SliceData;
    end
    close(h);
    
    % If SliceThickness is known, calculate the total Z extent of slices.
    image_meta_data.NumberOfSlices = length(slice_data);
    if isfield(slice_data(1), 'SliceThickness')
        image_meta_data.PhysicalTotalZ = ...
            double(slice_data(1).SliceThickness*length(slice_data));
    else
        no_slice_thickness = true;
    end
    
    % if PixelSpacing and SliceThickness is known, create
    % PhysicalAspectRatio.
    if ~no_pixel_spacing && ~no_slice_thickness
        image_meta_data.PhysicalAspectRatio = [...
            slice_data(1).PixelSpacing(1),...
            slice_data(1).PixelSpacing(2),...
            slice_data(1).SliceThickness...
            ];
    end
    
    % Save the data to the dicom directory.
%     h = waitbar(0, 'Saving mat files...', 'WindowStyle', 'modal');
%     waitbar(0,h);
%     

%     save(fullfile(dicom_directory,'VOLUME_IMAGE'),...
%     'volume_image', 'slice_data', 'image_meta_data');
%     local_directory = pwd;
%     eval(['cd ' dicom_directory]);
%     save VOLUME_IMAGE volume_image slice_data image_meta_data
%     eval(['cd ' local_directory]);
%     close(h);
clc
