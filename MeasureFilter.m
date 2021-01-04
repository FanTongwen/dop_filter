classdef MeasureFilter
    %MEASUREFILTER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        SateMap;
        MeasureData;
        MeasureSate_Count;% 卫星个数
    end
    
    methods
        function obj = MeasureFilter()
            %MEASUREFILTER 构造此类的实例
            %   此处显示详细说明
            obj.SateMap = containers.Map();
            obj.MeasureData = {};
            obj.MeasureSate_Count = 0;
        end
        

        
        function MF_Struct = MeasureFilter_Struct()
            MF_Struct.Mdoppler_Group = [];
            MF_Struct.MPseudo_Group = [];
            MF_Struct.PseudoValid = [];
            MF_Struct.Filter_Mdoppler = 0;
            MF_Struct.Filter_MPseudo = 0;
            
            MF_Struct.DiffPseudo_Group = [];
            MF_Struct.Filter_DeltPseudo = 0;
            MF_Struct.Pre_P = 0;
            
            MF_Struct.DiffDop_Group = [];
            MF_Struct.Filter_DeltDop = 0;
            MF_Struct.Pre_Dop = 0;
            MF_Struct.FitNum = 0;
        end
        
        function outputArg = Filter_DeltDopFunc(mf,satenum,temDouble)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            
            temdoulist = [];
            % temMfs = MeasureFilter_Struct();
            temint = 0;
            
            temlt = mf.SateMap(satenum);
            outputArg = temDouble;
            
        end
    end
end

