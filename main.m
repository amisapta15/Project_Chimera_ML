clear all
clc
close all
rng('default')
addpath(genpath(cd))
Data_Path='.\Data\';
ALL_Labels=cell(1,4);
for ii=1:4  
    model_names={'New_Kuramoto_March','Logistic','Hennon','FHN'};
    [trainX,trainY,testX]=read_chimera_data(Data_Path,model_names{1,ii});

    ntree=100;
    rng('default');
    %seed=RandStream('mcg16807','Seed',0);
    %RandStream.setGlobalStream(seed);
    
    mtry=round(sqrt(size(trainX,2)));
    
    [modelRaF,T12]=ObliqueRF_train(trainX,trainY,'ntrees',ntree,'nvartosample',mtry,'oblique',1);
    [YRaF,fv,Ts12,~]=ObliqueRF_predict(testX,modelRaF);
   % standard_RaF=length(find(YRaF==testY))/length(testY);
    
    [modelMPSVM_T,T12]=ObliqueRF_train(trainX,trainY,'ntrees',ntree,'nvartosample',mtry,'oblique',2);
    [YMPSVM_T,fv,Ts12,~]=ObliqueRF_predict(testX,modelMPSVM_T);
    %MPSVM_T=length(find(YMPSVM_T==testY))/length(testY);
    
    [modelMPSVM_P,T12]=ObliqueRF_train(trainX,trainY,'ntrees',ntree,'nvartosample',mtry,'oblique',3);
    [YMPSVM_P,fv,Ts12,~]=ObliqueRF_predict(testX,modelMPSVM_P);
   % MPSVM_P=length(find(YMPSVM_P==testY))/length(testY);
    
    [modelMPSVM_N,T12]=ObliqueRF_train(trainX,trainY,'ntrees',ntree,'nvartosample',mtry,'oblique',4);
    [YMPSVM_N,fv,Ts12,~]=ObliqueRF_predict(testX,modelMPSVM_N);
    %MPSVM_N=length(find(YMPSVM_N==testY))/length(testY);
    
    %% Optimize RVFL Parameters
    orginal_trainX=trainX; %% Saving for Later Use
    orginal_trainY=trainY;
    orginal_testX=testX;
    %orginal_testY=testY;
    
    dataX=trainX;
    [m,n]=size(dataX);
    dataY=trainY;
    
    %% Training RVFL-AE using 5 times 4-fold croos validation
    M=5;
    v=4;
    step=floor(m/v);
    num_count=M*v;
    TestingAccuracy=zeros(1,num_count);
    flag=0;
    max_acc=0;
    Range_C=[2^-6,2^-4,2^-2,2^0,2^2,2^4,2^6,2^8,2^10,2^12]; % 11 values
    Range_N=2:20:302; % Number of hidden Neurons
    option.Scale=1;
    option.method='RVFL_AE';
    for ic=1:length(Range_C)
        ic
        option.C=Range_C(ic);
        for jn=1:length(Range_N)
            option.N=Range_N(jn);
            count=1;
            flag=0;
            for i=1:M
                index=randperm(m);
                for j =1:v
                    if j~= v
                        flag=flag+1;
                        startpoint=(j-1)*step+1;
                        endpoint=(j)*step;
                    else
                        startpoint=(j-1)*step+1;
                        endpoint=m;
                    end
                    cv_p=startpoint:endpoint; %%%% test set position
                    
                    %%%%%%%%%%%%%% test set
                    testX=dataX(index(cv_p),:);
                    testY= dataY(index(cv_p),:);  %%%%label
                    %%%%%%%%%%%%%% training data
                    trainX=dataX;
                    trainX(index(cv_p),:)='';
                    trainY=dataY;
                    trainY(index(cv_p),:)='';
                    [model,train_time,train_accuracy,TestingAccuracy(count)]=RVFL_train_val_NEW(trainX,trainY,testX,testY,option);
                    count=count+1;
                end
            end
            if max_acc<mean(TestingAccuracy)
                max_acc=mean(TestingAccuracy);
                OptPara=option;
            end
        end
    end
    clear TestingAccuracy
    %% Testing RVFL-AE
    orginal_testY=zeros(size(orginal_testX,1),1);% Just passing it as dummy bcz we need not accuarcy but labels
    [model,train_time,train_accuracy,TestingAccuracy]=RVFL_train_val_NEW(orginal_trainX,orginal_trainY,orginal_testX,orginal_testY,option);
    ALL_Labels{1,ii}=[YRaF,YMPSVM_T,YMPSVM_P,YMPSVM_N,model.testY];
    save ALL_Labels 
end