%% 比较剔除误差前后的多普勒
clear
[txt_file,txt_path] = uigetfile(".\*.19o");
if txt_file ~= 0
    h_txt_file = fopen([txt_path,txt_file]);
    % 提取24行数据
    C = textscan(h_txt_file,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f");
    fclose(h_txt_file);
end

%% 提取信息
Date_Index=strcmp(C{1}, '>');
Data_Index1 = find(Date_Index);
Time_P = length(Data_Index1);
PRN = {};
time_p = {};%时间信息
dop_p = {};%多普勒信息
for i = 1:Time_P
    time_in_sec = C{5}(Data_Index1(i))*3600 + C{6}(Data_Index1(i))*60 + C{7}(Data_Index1(i));
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
%% 画图
sate_N = length(dop_p);
for i = 1:sate_N
    myplot(time_p{i}, dop_p{i}, 'doppler', PRN{i}, i);
end

%% function

%% 
clear
clc
[txt_file,txt_path] = uigetfile('E:\GNSSDATA\0829Second\IFData\QR3*.19o');
if txt_file ~= 0
    h_txt_file = fopen([txt_path,txt_file]);
    % 提取24行数据
    i = 0;
    while (true)
        head = textscan(h_txt_file,'%s',1);
        i = i+1;
        if i >200 || strcmp(head{1},'HEADER')
            break
        end
    end
    if i == 201
        errordlg("未找到文件头");
        return
    end
    C0 = textscan(h_txt_file,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f");
    fclose(h_txt_file);
end

[txt_file,txt_path] = uigetfile("E:\GNSSDATA\0829Second\IFData\QR4*.19o");
if txt_file ~= 0
    h_txt_file = fopen([txt_path,txt_file]);
    % 提取24行数据
    i = 0;
    while (true)
        head = textscan(h_txt_file,'%s',1);
        i = i+1;
        if i >200 || strcmp(head{1},'HEADER')
            break
        end
    end
    if i == 201
         errordlg("未找到文件头");
        return
    end
    C1 = textscan(h_txt_file,"%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f");
    fclose(h_txt_file);
end
%%
[time_p_0, dop_p_0, PRN_0] =  GetObsDop(C0);
[time_p_1, dop_p_1, PRN_1] =  GetObsDop(C1);
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

%% 
function [time_p, dop_p, PRN] =  GetObsDop(C)
    Date_Index=strcmp(C{1}, '>');
Data_Index1 = find(Date_Index);
Time_P = length(Data_Index1);
PRN = {};
time_p = {};%时间信息
dop_p = {};%多普勒信息
for i = 1:Time_P
    time_in_sec = C{5}(Data_Index1(i))*3600 + C{6}(Data_Index1(i))*60 + C{7}(Data_Index1(i));
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

