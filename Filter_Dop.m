%% 变量声明
clc
clear
% 常量
CLIGHT = 299792458.0;
FREQ1 = 1.57542E9;
FREQ1_CMP = 1.561098E9;

% Branch = 1;
Branch = 3;
save('Branch.mat', 'Branch');
% 变量
RefPrn = 19;%GPS 19作为起始参考星

% 读取数据
%[txt_file,txt_path] = uigetfile('/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20210103_1/QR5*.19o');
txt_file = 'QR5gnss2411.19o';
txt_path ='/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20210103_1/';
% 粗差剔除
time_p = {};%时间信息
dop_p = {};%多普勒信息
PRN = [];%卫星编号
MeasureFilter_Init();
h_txt_file = fopen([txt_path,txt_file]);
[GPSTime, n] = ReadObsHead(h_txt_file);%读取文件头
while(1)
    [obs, rs, Pos_INS, Vel_INS, dts, error] = ReadObsEpoch(h_txt_file, n);
    if error == 1 %读取结束
        msgbox(['end in' num2str(GPSTime)], '提示');
        fclose(h_txt_file);
        break;
    end
    % 预测多普勒
    for i = 1:min(n, 64)
        if obs{i}.openLoopFlag == 1
            continue;
        end
        
        Rel_Pos = rs(6*i-5 : 6*i - 3) - Pos_INS;% 相对位置
        Rk = sqrt(Rel_Pos*Rel_Pos.');% 估计的距离
        obs{i}.Pack_P = Rk - dts(i);
        obs{i}.clkbias = dts(i);
        L_Sate = Rel_Pos/Rk;% 方向矢量
        
        Rel_Vel = rs(6*i-2 : 6*i) - Vel_INS;% 相对速度
        DeltRk = Rel_Vel*L_Sate.';% 反视线方向上与相对速度的点积
        if obs{i}.sat > 32
            DeltRk = -1*DeltRk/(CLIGHT/FREQ1_CMP);
        else
            DeltRk = -1*DeltRk/(CLIGHT/FREQ1);
        end
        obs{i}.Pack_D = DeltRk;
    end
    
    for i = 1:min(n, 64) %仅限于GPS卫星
        if obs{i}.sat == RefPrn
            Index = i;
            Indexflag = 1;
        end
    end
    
    if (Indexflag == 1)&&(obs{Index}.SNR > 32.0)%确保有星且载噪比合适
        % 基站观测量保存
        if (sum(PRN == obs{Index}.sat)) == 0
            PRN = [PRN, obs{Index}.sat];
            time_p = [time_p, GPSTime];
            dop_p = [dop_p, obs{Index}.D];
        elseif (sum(PRN == obs{Index}.sat)) == 1
            time_p{PRN == obs{Index}.sat} = [time_p{PRN == obs{Index}.sat}, GPSTime];
            dop_p{PRN == obs{Index}.sat} = [dop_p{PRN == obs{Index}.sat}, obs{Index}.D];
        end
        
        for i = 1:min(n, 64)
            if (i == Index)
                continue;
            end
            if obs{i}.openLoopFlag == 1
                continue;
            end
            
            Star_INS_DDiff = (obs{Index}.D - obs{Index}.Pack_D) - (obs{i}.D - obs{i}.Pack_D);% 多普勒双差
            Star_INS_DDiff_s = MeasureFilter_DoplorFunc(obs{i}.sat, Star_INS_DDiff);% 平滑
            
            Fix_Dop = obs{i}.D + (Star_INS_DDiff - Star_INS_DDiff_s);% 修正的多普勒
            % 保存修正的多普勒信息
            if (sum(PRN == obs{i}.sat)) == 0
                PRN = [PRN, obs{i}.sat];
                time_p = [time_p, GPSTime];
                dop_p = [dop_p, Fix_Dop];
            elseif (sum(PRN == obs{i}.sat)) == 1
                time_p{PRN == obs{i}.sat} = [time_p{PRN == obs{i}.sat}, GPSTime];
                dop_p{PRN == obs{i}.sat} = [dop_p{PRN == obs{i}.sat}, Fix_Dop];
            end
        end
    else %换星
        [RefSate, Chooseindex] = Select_RefSateExtraSPP(obs, n);
        if RefSate ~= RefPrn
            if RefSate > 0 %换星成功
                
%                 % 换星后的滤波器变量修改
%                 for i = 1:64
%                     if RefSate == i
%                         continue;
%                     end
%                     if RefPrn == i
%                         continue;
%                     end
%                     ChangeRefStarDoppler(RefSate, i);
%                 end
%                 ChangeRefStarDoppler_E(RefSate, RefPrn);
                % 基站观测量保存
                if (sum(PRN == obs{Chooseindex}.sat)) == 0
                    PRN = [PRN, obs{Chooseindex}.sat];
                    time_p = [time_p, GPSTime];
                    dop_p = [dop_p, obs{Chooseindex}.D];
                elseif (sum(PRN == obs{Chooseindex}.sat)) == 1
                    time_p{PRN == obs{Chooseindex}.sat} = [time_p{PRN == obs{Chooseindex}.sat}, GPSTime];
                    dop_p{PRN == obs{Chooseindex}.sat} = [dop_p{PRN == obs{Chooseindex}.sat}, obs{Chooseindex}.D];
                end
                %改变滤波器状态
                
                RefPrn = RefSate;
                
                for i = 1:min(n, 64)
                    if (i == Chooseindex)
                        continue;
                    end
                    if obs{i}.openLoopFlag == 1
                        continue;
                    end
                    
                    Star_INS_DDiff = (obs{Chooseindex}.D - obs{Chooseindex}.Pack_D) - (obs{i}.D - obs{i}.Pack_D);% 多普勒双差
                    Star_INS_DDiff_s = MeasureFilter_DoplorFunc(obs{i}.sat, Star_INS_DDiff);% 平滑
                    
                    Fix_Dop = obs{i}.D + (Star_INS_DDiff - Star_INS_DDiff_s);% 修正的多普勒
                    % 保存修正的多普勒信息
                    if (sum(PRN == obs{i}.sat)) == 0
                        PRN = [PRN, obs{i}.sat];
                        time_p = [time_p, GPSTime];
                        dop_p = [dop_p, Fix_Dop];
                    elseif (sum(PRN == obs{i}.sat)) == 1
                        time_p{PRN == obs{i}.sat} = [time_p{PRN == obs{i}.sat}, GPSTime];
                        dop_p{PRN == obs{i}.sat} = [dop_p{PRN == obs{i}.sat}, Fix_Dop];
                    end
                end
            end
        end
    end
    
    % 读取下一次的时间及星数
    [GPSTime, n] = ReadObsTimeInfo(h_txt_file);
end
% 保存数据
sate_N = length(dop_p);
PRN_str = cell(1, sate_N);
for i = 1:sate_N
    if PRN(i) > 32
        PRN_str{i} = ['C', sprintf('%02d', PRN(i) - 32)];
    else
        PRN_str{i} = ['G', sprintf('%02d', PRN(i))];
    end
end

Fixed_Data.PRN = PRN_str;
Fixed_Data.time = time_p;
Fixed_Data.doppler = dop_p;
file_name = ['doppler_all_' num2str(Branch) '.mat'];
save(file_name, 'Fixed_Data');
%% 画图
sate_N = length(dop_p);
for i = 1:sate_N
    if length(time_p{i})<100
        continue
    end
    myplot(time_p{i}, dop_p{i}, 'doppler', num2str(PRN(i)), i);
end

%% 滤波器测试
clear
clc
MeasureFilter_Init();
Result = [];
InputDop = [];
for i = 1:1000
    temDop = rand();
    satenum = 31;
    Result = [Result, MeasureFilter_DoplorFunc(satenum, temDop)];
    InputDop = [InputDop, temDop];
end
plot(Result);
hold on
plot(InputDop);

%% 函数定义
% 多普勒滤波函数
function Filter_Result =  MeasureFilter_DoplorFunc(satenum, temDop)
    global SateMap;
    global MeasureData;
    global MeasureSate_Count;
    temdoulist = [];
    temMfs = MeasureFilter_Struct();
    temint = 0;
    
    temlt = find(SateMap{1} == satenum);
    
    if isempty(temlt) % 未找到该卫星
        SateMap{1} = [SateMap{1}, satenum];
        SateMap{2} = [SateMap{2}, MeasureSate_Count];
        MeasureData = [MeasureData, temMfs];
        temint = MeasureSate_Count;
        MeasureSate_Count = MeasureSate_Count + 1; 
    else
        temint = SateMap{2}(temlt);
    end
    
    if(length(MeasureData{temint}.Mdoppler_Group) < 20)
        MeasureData{temint}.Mdoppler_Group = [MeasureData{temint}.Mdoppler_Group, temDop];
        MeasureData{temint}.Filter_Mdoppler = mean(MeasureData{temint}.Mdoppler_Group);
    else
        if(abs(temDop - MeasureData{temint}.Filter_Mdoppler) < 4.0)
            MeasureData{temint}.Mdoppler_Group = [MeasureData{temint}.Mdoppler_Group(2:end), temDop];
            temdoulist = sort(MeasureData{temint}.Mdoppler_Group);
            MeasureData{temint}.Filter_Mdoppler = mean(temdoulist(9:13));%中间5个值进行平均
        end
    end
    Filter_Result = MeasureData{temint}.Filter_Mdoppler;
end

% 全局变量初始化
function MeasureFilter_Init()
global SateMap;
global MeasureData;
global  MeasureSate_Count;% 卫星个数

SateMap = {[],[]};
MeasureData = {};
MeasureSate_Count = 1;
end
% 释放内存
function MeasureFilter_Dinit()
clear SateMap
clear MeasureData
clear MeasureSate_Count
end

% 创建结构体
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
% 换星
function [RefSate, Chooseindex] = Select_RefSateExtraSPP(obs, n)
temCN01 = 0;
temCN02 = 0;
temel1 = 0;
temel2 = 0;
sateNum1 = 0;
sateNum2 = 0;
temRefXB1 = 0;
temRefXB2 = 0;
    for i = 1:min(n, 64)
        if obs{i}.openLoopFlag == 1
            continue;
        end
        
        if obs{i}.SNR > temCN01
            temCN01 = obs{i}.SNR;
            temel1 = obs{i}.Sate_Eleva;
            sateNum1 = obs{i}.sat;
            temRefXB1 = i;
        end
        
        if obs{i}.Sate_Eleva > temel2
            temCN02 = obs{i}.SNR;
            temel2 = obs{i}.Sate_Eleva;
            sateNum2 = obs{i}.sat;
            temRefXB2 = i;
        end
    end
    
    if (temCN01 > 30.0)&&(sateNum1 > 0)
        RefSate = sateNum1;
        Chooseindex = temRefXB1;
        return;
    elseif sateNum2 > 0
        RefSate = sateNum2;
        Chooseindex = temRefXB2;
        return;
    else
        Chooseindex = 0;
        RefSate = 0;
        return;
    end
end

function ChangeRefStarDoppler(Refnum, satenum)
    global SateMap;
    global MeasureData;
    temltRef = SateMap{1} == Refnum;
    temltSate = SateMap{1} == satenum;
    if (sum(temltRef) == 1)&&(sum(temltSate) == 1)
        size = min(length(MeasureData{SateMap{2}(temltSate)}.Mdoppler_Group), ...
            (length(MeasureData{temltRef}.Mdoppler_Group)));
        if size == 0
            return;
        end
        MeasureData{SateMap{2}(temltSate)}.Mdoppler_Group(1:size) = ...
            MeasureData{SateMap{2}(temltSate)}.Mdoppler_Group(1:size) - ...
            MeasureData{SateMap{2}(temltRef)}.Mdoppler_Group(1:size);
        
        MeasureData{SateMap{2}(temltSate)}.Filter_Mdoppler = ...
            MeasureData{SateMap{2}(temltSate)}.Filter_Mdoppler - ...
            MeasureData{SateMap{2}(temltRef)}.Filter_Mdoppler;
    end
end

function ChangeRefStarDoppler_E(Refnum, satenum)
    global SateMap;
    global MeasureData;
    temltRef = SateMap{1} == Refnum;
    temltSate = SateMap{1} == satenum;
    if (sum(temltRef) == 1)&&(sum(temltSate) == 1)
        size = min(length(MeasureData{SateMap{2}(temltSate)}.Mdoppler_Group), ...
            (length(MeasureData{temltRef}.Mdoppler_Group)));
        if size == 0
            return;
        end
        MeasureData{SateMap{2}(temltSate)}.Mdoppler_Group(1:size) = ...
            -MeasureData{SateMap{2}(temltRef)}.Mdoppler_Group(1:size);
        
        MeasureData{SateMap{2}(temltSate)}.Filter_Mdoppler = ...
            -MeasureData{SateMap{2}(temltRef)}.Filter_Mdoppler;
    end
end
% 读取数据相关函数
function [GPSTime, n] = ReadObsHead(h_txt_file)
C = textscan(h_txt_file,...
    "%s %f %f %f %f %f %f %f %f",...
    1,...
    'HeaderLine', 15);

n = C{9};
GPSTime = C{5}*3600 + C{6}*60 + C{7} +345600;
end
% 读取一个历元内的数据
function [obs, rs, Pos_INS, Vel_INS, dts, error] = ReadObsEpoch(h_txt_file, n)
error = 0;
if n == 0 % 无数据时
    obs = 0;
    rs = 0;
    Pos_INS = 0;
    Vel_INS = 0;
    dts = 0;
    error = 0;
    return;
else 
   C = textscan(h_txt_file,...
    "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
    n,...
    'Delimiter', ','); 
end

% 检测是否读完数据
for i = [2:6, 8:21] %信噪比可能是nan 不能作为判断标准
    if(sum(isnan(C{i})) > 0) %数据缺失,读取结束
        error = 1;
        obs = 0;
        rs = 0;
        Pos_INS = 0;
        Vel_INS = 0;
        dts = 0;
        return;
    end
end
for i = 1:21 %判断有无缺失行
    if(length(C{i}) < n) %数据缺失,读取结束
        error = 1;
        obs = 0;
        rs = 0;
        Pos_INS = 0;
        Vel_INS = 0;
        dts = 0;
        return;
    end
end
obs = cell(1, n);
rs = zeros(1, 6*n);
dts = zeros(1, n);
% Pos_INS
Pos_INS = [C{12}(1), C{13}(1), C{14}(1)];
% Vel_INS
Vel_INS = [C{18}(1), C{19}(1), C{20}(1)];

for i = 1:n
    % obs
    obs{i}.openLoopFlag = C{2}(i);
    obs{i}.sat = C{3}(i);
    obs{i}.P = C{4}(i);
    obs{i}.D = C{5}(i);
    obs{i}.L = C{6}(i);
    
    if isnan(C{7}(i))
        obs{i}.SNR = 0;
    else
        obs{i}.SNR = C{7}(i);
    end
    
    obs{i}.Sate_Eleva = C{8}(i);
    
    % rs
    rs(6*i-5 : 6*i) = [C{9}(i), C{10}(i), C{11}(i), C{15}(i), C{16}(i), C{17}(i)];

    % dts
    dts(i) = C{21}(i);
end
    

end
% 读取时间 星数信息
function [GPSTime, n] = ReadObsTimeInfo(h_txt_file)
C = textscan(h_txt_file,...
    "%s %f %f %f %f %f %f %f %f",...
    1);

n = C{9};
GPSTime = C{5}*3600 + C{6}*60 + C{7} +345600;
end

% 画图
function myplot1(X, Y, ylabel_s, title_str, n, color)
%PLOTWAV 此处显示有关此函数的摘要
%   此处显示详细说明

figure(n);
plot(X, Y, 'Color', color);
axis xy
axis tight
ylabel(ylabel_s);
xlabel('Time (secs)');
title(title_str);
end