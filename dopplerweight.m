%% 读取文件*Dopplerdiff.txt PseudoInfo.txt
clc
clear
[txt_file,txt_path] = uigetfile...
    ('/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/*Dopplerdiff.txt');
if txt_file ~= 0
    h_txt_file = fopen([txt_path,txt_file]);
    % 提取24行数据
    Dopplerdiff = textscan(h_txt_file, "%f %f", 'HeaderLine', 1, 'Delimiter', ',');
    fclose(h_txt_file);
end
% 
[txt_file,txt_path] = uigetfile...
    ('/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201229_PseAidCodeLoop_Timeratio/*PseudoInfo.txt');
if txt_file ~= 0
    h_txt_file = fopen([txt_path,txt_file]);
    % 提取24行数据
    PseudoInfo = textscan(h_txt_file, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", 'HeaderLine', 0, 'Delimiter', ',');
    fclose(h_txt_file);
end
% plot
myplot(Dopplerdiff{1}, Dopplerdiff{2}, 'Hz', '', 1);

myplot(PseudoInfo{1}, PseudoInfo{5}, 'tempCha[3]', 'StarDiff', 2);
myplot(PseudoInfo{1}, PseudoInfo{6}, 'tempCha[4]', 'Filter', 3);
myplot(PseudoInfo{1}, PseudoInfo{7}, 'tempCha[5]', 'fabs', 4);
myplot(PseudoInfo{1}, PseudoInfo{14}, 'tempCha[12]', 'obs[i].Pack\_D', 5);
myplot(PseudoInfo{1}, PseudoInfo{15}, 'tempCha[13]', 'obs[Index].Pack\_D', 6);
myplot(PseudoInfo{1}, PseudoInfo{16}, 'tempCha[14]', 'obs[i].D[0]', 7);
myplot(PseudoInfo{1}, PseudoInfo{17}, 'tempCha[15]', 'obs[Index].D[0]', 8);

%% 读取单差文件并保存为 SDTest_Data.mat
clear
SingleDiffFile = ...
    dir('/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/DopplerdoublediffQR3/*singlediff*Dopplerdiff.txt');

for i = 1:length(SingleDiffFile)
    h_txt_file = fopen([SingleDiffFile(i).folder filesep SingleDiffFile(i).name]);
    
    Dopplerdiff = textscan(h_txt_file, "%f %f", 'HeaderLine', 1, 'Delimiter', ',');
    if SingleDiffFile(i).name(1) == 'G'
        Celltmp = textscan(SingleDiffFile(i).name,'GPSsinglediff %d Dopplerdiff.txt');
        PRNtmp{i} = sprintf('G%02d', Celltmp{1});
        timetmp{i} = Dopplerdiff{1};
        doptmp{i} = Dopplerdiff{2};
    elseif SingleDiffFile(i).name(1) == 'B'
        Celltmp = textscan(SingleDiffFile(i).name,'BDsinglediff %d Dopplerdiff.txt');
        PRNtmp{i} = sprintf('C%02d', Celltmp{1});
        timetmp{i} = Dopplerdiff{1};
        doptmp{i} = Dopplerdiff{2};
    end
    fclose(h_txt_file);
end
SDTest_Data = GetDataStruct(timetmp, doptmp, PRNtmp);
save('SDTest_Data.mat', 'SDTest_Data');
%% 读取双差文件并保存为 DDTtst_Data.mat
clear
DoubleDiffFile = ...
    dir('/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/DopplerdoublediffQR3/*_*Dopplerdiff.txt');
for i = 1:length(DoubleDiffFile)
    h_txt_file = fopen([DoubleDiffFile(i).folder filesep DoubleDiffFile(i).name]);
    
    Dopplerdiff = textscan(h_txt_file, "%f %f", 'HeaderLine', 1, 'Delimiter', ',');
    if DoubleDiffFile(i).name(1) == 'G'
        Celltmp = textscan(DoubleDiffFile(i).name,'GPS %d _ %d Dopplerdiff.txt');
        PRNtmp{i} = sprintf('G%02d_%02d', Celltmp{1}, Celltmp{2});
        timetmp{i} = Dopplerdiff{1}.';
        doptmp{i} = Dopplerdiff{2}.';
    elseif DoubleDiffFile(i).name(1) == 'B'
        Celltmp = textscan(DoubleDiffFile(i).name,'BD %d _ %d Dopplerdiff.txt');
        PRNtmp{i} = sprintf('C%02d_%02d', Celltmp{1}, Celltmp{2});
        timetmp{i} = Dopplerdiff{1}.';
        doptmp{i} = Dopplerdiff{2}.';
    end
    fclose(h_txt_file);
end
DDTest_Data = GetDataStruct(timetmp, doptmp, PRNtmp);
save('DDTest_Data.mat', 'DDTest_Data');
%% base station data

%% 读取obs文件并保存为mat格式
clc
clear
% QR3
%{
[txt_file,txt_path] = uigetfile('/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/QR3*.19o');
if txt_file ~= 0
    h_txt_file = fopen([txt_path,txt_file]);
    C0 = textscan(h_txt_file,...
        "%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
        'HeaderLine', 15);
    fclose(h_txt_file);
end
%}
txt_file = 'QR3gnss2411.19o';
txt_path = '/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/';
h_txt_file = fopen([txt_path,txt_file]);
C0 = textscan(h_txt_file,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
        'HeaderLine', 15);
% QR0
%{
[txt_file,txt_path] = uigetfile('/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/QR0*.19o');
if txt_file ~= 0
    h_txt_file = fopen([txt_path,txt_file]);
    C1 = textscan(h_txt_file,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
        'HeaderLine', 15);
    fclose(h_txt_file);
end
%}
txt_file = 'QR0gnss2411.19o';
txt_path = '/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/';
h_txt_file = fopen([txt_path,txt_file]);
C1 = textscan(h_txt_file,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
        'HeaderLine', 15);
% Ture Value
txt_file = 'QR3gnss2411.19o';
txt_path = '/media/ftw/diske/GNSSDATA/0829Second/IFData/Result_20201127_PseAidCodeLoop_Timeratio/PseudoAnalyze_Result/';
h_txt_file = fopen([txt_path,txt_file]);
C2 = textscan(h_txt_file,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",...
        'HeaderLine', 15);


[time_p_0, dop_p_0, PRN_0] =  GetObsDop(C0);
[time_p_1, dop_p_1, PRN_1] =  GetObsDop(C1);
[time_p_2, dop_p_2, PRN_2] =  GetObsDop(C2);
%% 保存
QR3_Data = GetDataStruct(time_p_0, dop_p_0, PRN_0);
QR0_Data = GetDataStruct(time_p_1, dop_p_1, PRN_1);
True_Data = GetDataStruct(time_p_2, dop_p_2, PRN_2);
save('doppler_01.mat','QR3_Data', 'QR0_Data', 'True_Data');
%% 画图
sate_N = length(dop_p_0);
for i = 1:sate_N
    if length(time_p_0{i})<100
        continue
    end
    myplot1(time_p_0{i}, dop_p_0{i}, 'doppler', PRN_0{i}, i, 'r');
    hold on
    myplot1(time_p_1{i}, dop_p_1{i}, 'doppler', PRN_1{i}, i, 'blue');
    hold on
    legend('剔除误差后','剔除误差前');
end

%% 函数定义
% 得到obs文件内的数据 [时间序列,多普勒,卫星编号]
function [time_p, dop_p, PRN] =  GetObsDop(C)
Date_Index=strcmp(C{1}, '>');
Data_Index1 = find(Date_Index);
Time_P = length(Data_Index1);
PRN = {};
time_p = {};%时间信息
dop_p = {};%多普勒信息
for i = 1:Time_P
    time_in_sec = C{5}(Data_Index1(i))*3600 + C{6}(Data_Index1(i))*60 + C{7}(Data_Index1(i)) + 345600;%GPS time
    sate_num = C{9}(Data_Index1(i));%卫星数量
    if sate_num ~= 0
        for j = 1:sate_num
            if (Data_Index1(i) + j) > (length(C{1}) - 1) break; end
            temp1 = strcmp(PRN, C{1}{Data_Index1(i) + j});
            if sum(temp1) == 0 %没有该PRN卫星的信息
                PRN = [PRN, C{1}{Data_Index1(i) + j}];
                time_p = [time_p, time_in_sec];
                
                if isnan(C{5}(Data_Index1(i) + j))
                    dop_p = [dop_p, C{3}(Data_Index1(i) + j)];%无载波相位
                else
                    dop_p = [dop_p, C{4}(Data_Index1(i) + j)];
                end
            elseif sum(temp1) == 1%已经记录过该PRN的信息
                time_p{temp1} = [time_p{temp1}, time_in_sec];
                if isnan(C{5}(Data_Index1(i) + j))
                    dop_p{temp1} = [dop_p{temp1}, C{3}(Data_Index1(i) + j)];
                else
                    dop_p{temp1} = [dop_p{temp1}, C{4}(Data_Index1(i) + j)];
                end
            end
        end
    end
end
end
function myplot1(X, Y, ylabel_s, title_str, n, color)
%MYPLOT1 画图函数
%   X: x轴数据
%   Y: y轴数据
%   y_label_s: y轴单位
%   title_str: 图标题
%   n: 图编号
%   color: 线颜色

figure(n);
plot(X, Y, 'Color', color);
axis xy
axis tight
ylabel(ylabel_s);
xlabel('Time (secs)');
title(title_str);
end
function [Data] = GetDataStruct(time, dop, PRN)
Data.time = time;
Data.PRN =PRN;
Data.doppler = dop;
end