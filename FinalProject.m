% Project by: Aurora Guidubaldi, Rafael Sêco and Henrique Oliveira

clear all

try
    load dataROI_E_N.mat
     allVarSubRoi=who('-file','dataROI_E_N');%List of the names of all the variable of dataROI_E_N.mat
     nrSubRois = size(allVarSubRoi,1);   % Number of subROI's
catch Error
      msgbox('Mat File Not Found','Mat File Error','error');
    return;
end

Table_features_1st_order_Normal = table([],[],[],[],[],[], 'VariableNames', {'Media','Variance','Symmetry','Kurtosis','Energy','Entropy'});
Table_features_1st_order_Emphysema = table([],[],[],[],[],[], 'VariableNames', {'Media','Variance','Symmetry','Kurtosis','Energy','Entropy'});
Table_CoOccurrence_Matrix_Descriptors = table([],[],[],[],[], 'VariableNames',{'Correlation','Contrast','Uniformity','Homogeneity','Entropy'});
FeaturesOfTexture_Normal = table([],[],[],[],[], 'VariableNames', {'SRE','LRE','GLNU','RLNU','RP'});
FeaturesOfTexture_Emphysema = table([],[],[],[],[], 'VariableNames', {'SRE','LRE','GLNU','RLNU','RP'});

TableEntropy = table([],'VariableNames',{'Entropy'});
Table_CoOccurenceMatrix_Normal = table([],[],[],[],[],'VariableNames',{'Correlation','Contrast','Energy','Homogeneity','Entropy'});
Table_CoOccurenceMatrix_Emphysema = table([],[],[],[],[],'VariableNames',{'Correlation','Contrast','Energy','Homogeneity','Entropy'});

Table_4thmethod_Normal = table([],[],[],[],[],'VariableNames',{'Variance','Angular Second Moment','Inverse Difference Moment','Correlation','Contrast'});
Table_4thmethod_Emphysema = table([],[],[],[],[],'VariableNames',{'Variance','Angular Second Moment','Inverse Difference Moment','Correlation','Contrast'});

GL1(1) = 3000;
GL1(2) = -3000;

for nrSubRoiCycle=1:1:nrSubRois %Cycle that goes through every subROI
   
    expEval = sprintf('cellDataSubROI = %s;',allVarSubRoi{nrSubRoiCycle,1});%Name of the Variable SubRoi analysed in each cycle 
    eval(expEval);  %Execute the instruction
    percentage = cellDataSubROI{1,1};

     if(percentage > 99)
        data = cellDataSubROI{1,3}; % Array with information of each subROI,

        if(min(data(:,3))<GL1(1))
            GL1(1)=min(data(:,3));
        end

        if(max(data(:,3))>GL1(2))
            GL1(2)=max(data(:,3));
        end
     end
end

 for nrSubRoiCycle=1:1:nrSubRois %Cycle that goes through every subROI
   
    expEval = sprintf('cellDataSubROI = %s;',allVarSubRoi{nrSubRoiCycle,1});%Name of the Variable SubRoi analysed in each cycle 
    eval(expEval);  %Execute the instruction 
    Pathology = cellDataSubROI{1,2};
    percentage = cellDataSubROI{1,1};

    if(percentage > 99)
        data = cellDataSubROI{1,3}; % Array with information of each subROI,

        xx=data(:,1);
        yy=data(:,2);
        values=data(:,3);

        xx1= xx-min(xx(:))+1;
        yy1 = yy-min(yy(:))+1;

        Organizeddata = zeros(max(xx1(:)),max(yy1(:)));
        for k=1:numel(xx1)
            lin=xx1(k);
            col=yy1(k);
            Organizeddata(lin,col) =  values(k);
        end

        % Dividing the normal people, from the pathological ones.
        if strcmp(Pathology,'Normal')

            % Histogram
            HSSubRoiNormal=data(:,3); % Hounsfield Scale of the Normal patient
            L=32;
            GSSubRoiNormal = scaleNL(HSSubRoiNormal,L,GL1);
            Histogram(GSSubRoiNormal,L,nrSubRoiCycle+1,'Normalized Histogram for a Normal patient');

            % Features First Order
            dataFirstOrderNormal = featuresFirstOrder(GSSubRoiNormal,L);
            TableSubRoi = array2table(dataFirstOrderNormal,'VariableNames',{'Media','Variance','Symmetry','Kurtosis','Energy','Entropy'});
            Table_features_1st_order_Normal = vertcat(Table_features_1st_order_Normal,TableSubRoi);

            % Co-occurrence Matrix
            GSOrganizedDataNormal = scaleNL(Organizeddata,L,GL1);
            CCMatrixNormal = graycomatrix(GSOrganizedDataNormal,'Offset',[0 1],'NumLevels',32,'GrayLimits',[],'Symmetric',true);

            % Graycoprops + Entropy
            statsNormal = graycoprops(CCMatrixNormal,{'Correlation','Contrast','Energy','Homogeneity'});
            Entropy_Normal = entropyProject(CCMatrixNormal,32);
            T1 = struct2table(statsNormal);
            T2 = array2table(Entropy_Normal,'VariableNames',{'Entropy'});
            CoOccurrence_Table_Normal = horzcat(T1,T2);
            Table_CoOccurenceMatrix_Normal = [Table_CoOccurenceMatrix_Normal;CoOccurrence_Table_Normal];


            %Primitives Matrix
            escala1 = min(GSOrganizedDataNormal(:));
            escala2 = max(GSOrganizedDataNormal(:));
            Primitives_Matrix_Normal = grayrlmatrix(uint8(GSOrganizedDataNormal),'NumLevels',32,'Offset',1,'GrayLimits',[escala1 escala2]);
            Primitives_Matrix_Normal = Primitives_Matrix_Normal{1,1};
            mask = ones(size(Primitives_Matrix_Normal));
            descriptors_Primitives_Normal = glrlm(Primitives_Matrix_Normal,32,mask);
            Table = array2table(descriptors_Primitives_Normal,'VariableNames',{'SRE','LRE','GLNU','RLNU','RP'});
            FeaturesOfTexture_Normal = [FeaturesOfTexture_Normal;Table];



            %Gray Level Difference Method
            dataGLDMNormal = GreyLevelDifferenceMethod(GSOrganizedDataNormal,L);
            TableSubRoi4th = array2table(dataGLDMNormal,'VariableNames',{'Variance','Angular Second Moment','Inverse Difference Moment','Correlation','Contrast'});
            Table_4thmethod_Normal = vertcat(Table_4thmethod_Normal,TableSubRoi4th);

        else

            % Histogram
            HSSubRoiEmphysema=data(:,3); % Hounsfield Scale of the Normal patient
            L=32;
            GSSubRoiEmphysema = scaleNL(HSSubRoiEmphysema,L,GL1);
            Histogram(GSSubRoiEmphysema,L,nrSubRoiCycle+1,'Normalized Histogram for an emphysema patient');

            % Features First Order
            dataFirstOrderEmphysema = featuresFirstOrder(GSSubRoiEmphysema,L);
            TableSubRoi = array2table(dataFirstOrderEmphysema,'VariableNames',{'Media','Variance','Symmetry','Kurtosis','Energy','Entropy'});
            Table_features_1st_order_Emphysema = vertcat(Table_features_1st_order_Emphysema,TableSubRoi);

            % Co-occurrence Matrix
            GSOrganizedDataEmphysema = scaleNL(Organizeddata,L,GL1);
            CCMatrixEmphysema = graycomatrix(GSOrganizedDataEmphysema,'Offset',[0 1],'NumLevels',32,'GrayLimits',[],'Symmetric',true);

            % Graycoprops
            statsEmphysema = graycoprops(CCMatrixEmphysema,{'Correlation','Contrast','Energy','Homogeneity'});
            Entropy_Emphysema = entropyProject(CCMatrixEmphysema,32);
            T1 = struct2table(statsEmphysema);
            T2 = array2table(Entropy_Emphysema,'VariableNames',{'Entropy'});
            CoOccurrence_Table_Emphysema = horzcat(T1,T2);
            Table_CoOccurenceMatrix_Emphysema = [Table_CoOccurenceMatrix_Emphysema;CoOccurrence_Table_Emphysema];


            %Primitives Matrix
            escala1 = min(GSOrganizedDataEmphysema(:));
            escala2 = max(GSOrganizedDataEmphysema(:));
            Primitives_Matrix_Emphysema = grayrlmatrix(uint8(GSOrganizedDataEmphysema),'NumLevels',32,'Offset',1,'GrayLimits',[escala1 escala2]);
            Primitives_Matrix_Emphysema = Primitives_Matrix_Emphysema{1,1};
            mask = ones(size(Primitives_Matrix_Emphysema));
            descriptors_Primitives_Emphysema = glrlm(Primitives_Matrix_Emphysema,32,mask);

            Table = array2table(descriptors_Primitives_Emphysema,'VariableNames',{'SRE','LRE','GLNU','RLNU','RP'});
            FeaturesOfTexture_Emphysema = [FeaturesOfTexture_Emphysema;Table];



            %Gray Level Difference Method
            dataGLDMEmphysema = GreyLevelDifferenceMethod(GSOrganizedDataEmphysema,L);
            TableSubRoi4th = array2table(dataGLDMEmphysema,'VariableNames',{'Variance','Angular Second Moment','Inverse Difference Moment','Correlation','Contrast'});
            Table_4thmethod_Emphysema = vertcat(Table_4thmethod_Emphysema,TableSubRoi4th);

        end
    end


 end
 
%  Create an excel file for Normal and Emphysema Patients
writetable(Table_features_1st_order_Normal,"Table_Features1stOrder_NormalPatients.xlsx")
writetable(Table_features_1st_order_Emphysema,"Table_Features1stOrder_EmphysemaPatients.xlsx")

writetable(Table_CoOccurenceMatrix_Normal,"CoOccurrence_Matrix_Attributes_NormalPatients.xlsx")
writetable(Table_CoOccurenceMatrix_Emphysema,"CoOccurrence_Matrix_Attributes_EmphysemaPatients.xlsx")

writetable(FeaturesOfTexture_Normal,"FeaturesOfTexture_NormalPatients.xlsx")
writetable(FeaturesOfTexture_Emphysema,"FeaturesOfTexture_EmphysemaPatients.xlsx")

writetable(Table_4thmethod_Normal,"FeaturesOfGrayLevelDifferenceMethod_NormalPatients.xlsx")
writetable(Table_4thmethod_Emphysema,"FeaturesOfGrayLevelDifferenceMethod_EmphysemaPatients.xlsx")


% Average of every feature for Normal Patients
Average_features_1st_order_Normal = varfun(@mean, Table_features_1st_order_Normal, 'InputVariables',@isnumeric);
Average_features_1st_order_Emphysema = varfun(@mean, Table_features_1st_order_Emphysema, 'InputVariables',@isnumeric);

Average_CoOccurenceMatrix_Normal = varfun(@mean,Table_CoOccurenceMatrix_Normal, 'InputVariables',@isnumeric);
Average_CoOccurenceMatrix_Emphysema = varfun(@mean, Table_CoOccurenceMatrix_Emphysema, 'InputVariables',@isnumeric);

Average_FeaturesOfTexture_Normal = varfun(@mean, FeaturesOfTexture_Normal, 'InputVariables',@isnumeric);
Average_FeaturesOfTexture_Emphysema = varfun(@mean, FeaturesOfTexture_Emphysema, 'InputVariables',@isnumeric);

Average_FeaturesGLDM_Normal = varfun(@mean, Table_4thmethod_Normal, 'InputVariables',@isnumeric);
Average_FeaturesGLDM_Emphysema = varfun(@mean, Table_4thmethod_Emphysema, 'InputVariables',@isnumeric);

% Standard Deviation
SDFFONormal =  std(Table_features_1st_order_Normal{:,:});
SDFFOEmphysema = std(Table_features_1st_order_Emphysema{:,:});
SDCoOccurenceNormal = std(Table_CoOccurenceMatrix_Normal{:,:});
SDCoOccurenceEmphysema = std(Table_CoOccurenceMatrix_Emphysema{:,:});
SDFOTNormal = std(FeaturesOfTexture_Normal{:,:});
SDFOTEmphysema = std(FeaturesOfTexture_Emphysema{:,:});
SDGLDMNormal = std(Table_4thmethod_Normal{:,:});
SDGLDMEmphysema = std(Table_4thmethod_Emphysema{:,:});

%Plotting the bar graphs:

%Features 1st Order
figure('Name','Comparison between the respective averages of the first order features from Normal patients and Emphysema patients','NumberTitle','off')
X = categorical({'Media','Variance','Symmetry','Kurtosis','Energy','Entropy'});
X = reordercats(X,{'Media','Variance','Symmetry','Kurtosis','Energy','Entropy'});
Y = [Average_features_1st_order_Normal.Variables;SDFFONormal;Average_features_1st_order_Emphysema.Variables;SDFFOEmphysema];

graph = bar(X,Y);

xtips1 = graph(1).XEndPoints;
ytips1 = graph(1).YEndPoints;
labels1 = string(graph(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = graph(2).XEndPoints;
ytips2 = graph(2).YEndPoints;
labels2 = string(graph(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = graph(3).XEndPoints;
ytips3 = graph(3).YEndPoints;
labels3= string(graph(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips4 = graph(4).XEndPoints;
ytips4 = graph(4).YEndPoints;
labels4 = string(graph(4).YData);
text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

set(graph,{'DisplayName'},{'Average for Normal Patients','Standard Deviation for Normal Patients','Average for Emphysema Patients','Standard Deviation for Emphysema Patients'}')
title('Comparison between the respective averages of the first order features from Normal patients and Emphysema patients')
legend()


%Features Co-Occurrence Matrix
figure('Name','Comparison between the respective averages of the Co-Occurrence Matrix features from Normal patients and Emphysema patients','NumberTitle','off')
X = categorical({'Correlation','Contrast','Energy','Homogeneity','Entropy'});
X = reordercats(X,{'Correlation','Contrast','Energy','Homogeneity','Entropy'});
Y = [Average_CoOccurenceMatrix_Normal.Variables ; SDCoOccurenceNormal;Average_CoOccurenceMatrix_Emphysema.Variables;SDCoOccurenceEmphysema];

graph = bar(X,Y);

xtips1 = graph(1).XEndPoints;
ytips1 = graph(1).YEndPoints;
labels1 = string(graph(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = graph(2).XEndPoints;
ytips2 = graph(2).YEndPoints;
labels2 = string(graph(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = graph(3).XEndPoints;
ytips3 = graph(3).YEndPoints;
labels3= string(graph(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips4 = graph(4).XEndPoints;
ytips4 = graph(4).YEndPoints;
labels4 = string(graph(4).YData);
text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

set(graph,{'DisplayName'},{'Normal Patient','SD Normal','Emphysema Patient','SD Emphysema'}')
title('Comparison between the respective averages of the Co-Occurrence Matrix features from Normal patients and Emphysema patients')
legend()

%Features of Texture (primitives)
figure('Name','Comparison between the respective averages of the primitive features of texture from Normal patients and Emphysema patients','NumberTitle','off')
X = categorical({'SRE','LRE','GLNU','RLNU','RP'});
X = reordercats(X,{'SRE','LRE','GLNU','RLNU','RP'});
Y = [Average_FeaturesOfTexture_Normal.Variables ; SDFOTNormal;Average_FeaturesOfTexture_Emphysema.Variables;SDFOTEmphysema];

graph = bar(X,Y);

xtips1 = graph(1).XEndPoints;
ytips1 = graph(1).YEndPoints;
labels1 = string(graph(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = graph(2).XEndPoints;
ytips2 = graph(2).YEndPoints;
labels2 = string(graph(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = graph(3).XEndPoints;
ytips3 = graph(3).YEndPoints;
labels3= string(graph(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips4 = graph(4).XEndPoints;
ytips4 = graph(4).YEndPoints;
labels4 = string(graph(4).YData);
text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

set(graph,{'DisplayName'},{'Normal Patient','SD Normal','Emphysema Patient','SD Emphysema'}')
title('Comparison between the respective averages of the primitive features of texture from Normal patients and Emphysema patients')
legend()

%Features of Gray Level Difference Method
figure('Name','Comparison between the respective averages of the Gray Level Difference Method features of texture from Normal patients and Emphysema patients','NumberTitle','off')
X = categorical({'Variance','ASM','IDM','Correlation','Contrast'});
X = reordercats(X,{'Variance','ASM','IDM','Correlation','Contrast'});
Y = [Average_FeaturesGLDM_Normal.Variables ; SDGLDMNormal;Average_FeaturesGLDM_Emphysema.Variables;SDGLDMEmphysema];

graph = bar(X,Y);

xtips1 = graph(1).XEndPoints;
ytips1 = graph(1).YEndPoints;
labels1 = string(graph(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = graph(2).XEndPoints;
ytips2 = graph(2).YEndPoints;
labels2 = string(graph(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips3 = graph(3).XEndPoints;
ytips3 = graph(3).YEndPoints;
labels3= string(graph(3).YData);
text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips4 = graph(4).XEndPoints;
ytips4 = graph(4).YEndPoints;
labels4 = string(graph(4).YData);
text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

set(graph,{'DisplayName'},{'Normal Patient','SD Normal','Emphysema Patient','SD Emphysema'}')
title('Comparison between the respective averages of the Grey Level Difference Method of texture from Normal patients and Emphysema patients')
legend()
