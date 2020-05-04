function [trainX,trainY,testX]=read_chimera_data(path1,data_model_name)
%Data_Unmarked as Testing Data
data_path_org=[path1 data_model_name '\Training_Data\'];
data_path_org_test=[path1 data_model_name '\Testing_Data\'];
data_test_labels=[path1 data_model_name '\Testing_Data_Label\Labels.xlsx'];
%% CHIMERA DATA
data_path=strcat(data_path_org,'chimera\');
num_files=dir(data_path);
chimera_data=[];
for i=3:length(num_files)
    data_file=strcat(data_path,num_files(i).name);
    file_rd=importdata(data_file);
    chimera_data=[chimera_data,file_rd];
end
num_samples=length(num_files)-2; % Beacuse the first two are directories
chimera_label=zeros(num_samples,1);
chimera_label(:)=1;
chimera_d=[chimera_data',chimera_label];
%% COHERENT DATA
data_path=strcat(data_path_org,'coherent\');
num_files=dir(data_path);
coherent_data=[];
for i=3:length(num_files)
    data_file=strcat(data_path,num_files(i).name);
    file_rd=importdata(data_file);
    coherent_data=[coherent_data,file_rd];
end
num_samples=length(num_files)-2; % Beacuse the first two are directories
coherent_label=zeros(num_samples,1);
coherent_label(:)=2;%data label
coherent_d=[coherent_data',coherent_label];
%% INCOHERENT DATA
data_path=strcat(data_path_org,'incoherent\');
num_files=dir(data_path);
incoherent_data=[];
for i=3:length(num_files)
    data_file=strcat(data_path,num_files(i).name);
    file_rd=importdata(data_file);
    incoherent_data=[incoherent_data,file_rd];
end
num_samples=length(num_files)-2; % Beacuse the first two are directories
incoherent_label=zeros(num_samples,1);
incoherent_label(:)=3;%data label
incoherent_d=[incoherent_data',incoherent_label];

all_data=[chimera_d;coherent_d;incoherent_d];
trainX=all_data(:,1:end-1);
trainY=all_data(:,end);

clear chimera_data coherent_data incoherent_data  chimera_d file_rd
%% Load Testing Data
test_num_files = dir(fullfile(data_path_org_test,'*.dat'));
sorted_file_names = natsortfiles({test_num_files.name});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chimera_data=[];
for i=1:length(sorted_file_names)
    data_file=strcat(data_path_org_test,sorted_file_names(i));
    file_rd=importdata(char(data_file));
    chimera_data=[chimera_data,file_rd];
end
num_samples=length(sorted_file_names);
testX=chimera_data';
end