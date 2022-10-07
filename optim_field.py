import geatpy as ga
import numpy as np
import pandas as pd
import time
import pymysql
from sqlalchemy import create_engine
import os

MINVALUE = -2000  # 若不满足该约束条件，则设置为该值
# 染色体数须大于等于2

def filter_function(Intialchrom, diction, AMN, row_number, nind, maxleng, standard_speed_difference,big_row, lind, varlen, safe_distance, digraph, diction2, num_property):
    # 根据初始种群和 diction 将每个位置料堆对应的属性整合成一个list
    chrom_property = np.array(diction[int(Intialchrom[0][0] - 1), :])
    for i in range(len(Intialchrom) - 1):
        chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
    for j in range(len(Intialchrom[1]) - 1):
        temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
        for i in range(len(Intialchrom) - 1):
            temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
        chrom_property = np.column_stack((chrom_property, temporary))

    # 计算固定搭配组合的适应度
    d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    # 计算是否又因为固定搭配组合放到一个料堆导致的放不下的情况， 若有将固定搭配组合中最长的料堆拆成两个料堆再进行计算
    while max(d) == MINVALUE:
        temporary = diction[:, 4]
        for i,ii in enumerate(temporary):
            if diction[i][5] == 0:
                temporary[i] = 0
        temporary = np.argmax(temporary)   # 有固定搭配的最长物料的序号
        relate_temp = int(diction[temporary][5])  # 对应的固定搭配物料的序号
        replace_ind = [3,4,9,10,12,13,15,16,18,19,21,22]
        for ind_item in replace_ind:
            diction[temporary][ind_item] = int(diction[temporary][ind_item] / 2)
            diction[relate_temp][ind_item] = int(diction[relate_temp][ind_item] / 2)
        diction[temporary][-1] = len(diction)
        diction = np.row_stack((diction, diction[temporary]))
        diction[relate_temp][-1] = len(diction)
        diction = np.row_stack((diction, diction[temporary]))

        # 调整
        lind += 2
        big_row += 2
        varlen = len(diction)
        if big_row >= AMN:
            row_number += 1
            big_row %= AMN

        # 判断分裂的堆数太多并且无解情况
        lea_space = safe_distance*(AMN+big_row+row_number*AMN)
        if sum(maxleng)< lea_space:
            print("Please try reset parameters.")
            raise Exception ("No answer!")

        Intialchrom = ga.crtpp(nind, lind, varlen)
        chrom_property = np.array(diction[int(Intialchrom[0][0] - 1), :])
        for i in range(len(Intialchrom) - 1):
            chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
        for j in range(len(Intialchrom[1]) - 1):
            temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
            for i in range(len(Intialchrom) - 1):
                temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
            chrom_property = np.column_stack((chrom_property, temporary))
        d = (chrom_property, diction, AMN, row_number, nind, big_row, num_property)

    # 计算料堆容量约束
    a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)
    # print(max(a))
    #________________________
    # df = pd.DataFrame([Intialchrom[label_2],label_list_2])
    # df.to_csv("C:/Users/Dorothy.S/Desktop/the_result.csv",encoding='gbk')
    #————————————————————————————

    # 检查是否存在某个料堆过长导致现有的平均堆数放到每个料条的情况难以实现的情况，若出现则将其看作两个料堆
    while max(a) == MINVALUE:
        temporary = diction[:, 4]
        temporary = np.argmax(temporary)
        replace_ind = [3,4,9,10,12,13,15,16,18,19,21,22]
        for ind_item in replace_ind:
            diction[temporary][ind_item] = int(diction[temporary][ind_item] / 2)
        diction[temporary][-1] = len(diction)
        diction = np.row_stack((diction, diction[temporary]))
        lind += 1
        big_row += 1
        varlen = len(diction)
        if big_row == AMN:
            row_number += 1
            big_row = 0

        # 判断分裂的堆数太多并且无解情况
        lea_space = safe_distance*(AMN+big_row+row_number*AMN)
        if sum(maxleng)< lea_space:
            print("Please try reset parameters.")
            raise Exception ("No answer!")

        Intialchrom = ga.crtpp(nind, lind, varlen)
        chrom_property = np.array(diction[int(Intialchrom[0][0] - 1), :])
        for i in range(len(Intialchrom) - 1):
            chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
        for j in range(len(Intialchrom[1]) - 1):
            temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
            for i in range(len(Intialchrom) - 1):
                temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
            chrom_property = np.column_stack((chrom_property, temporary))
        a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)

    # 计算分散搭配约束的适应度
    c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    # 计算位置偏好约束和产品相容性约束的适应度
    e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)
    f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    # 合并五条适应度的结果，选取同时满足五条适应度的方案放入新的初始种群
    IntialObjV = a + c + d + e + f
    chrom, ObjV = [], []
    for i,ii in enumerate(IntialObjV):
        if ii == 0:
            chrom.append(list(Intialchrom[i]))
            ObjV.append(list(IntialObjV[i]))

    # 重复生成初始种群并计算五条适应度直到新的初始种群的大小不小于原来的大小
    # Intialchrom = ga.crtpp(10*nind, lind, varlen)
    # chrom_property = list(diction[int(Intialchrom[0][0] - 1), :])
    # for i in range(len(Intialchrom) - 1):
    #     chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
    # for j in range(len(Intialchrom[1]) - 1):
    #     temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
    #     for i in range(len(Intialchrom) - 1):
    #         temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
    #     chrom_property = np.column_stack((chrom_property, temporary))

    # a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)
    # c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
    # d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
    # e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)
    # f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)
    # IntialObjV = a + c + d + e + f
    # for i, ii in enumerate(IntialObjV):
    #     temp_list = [float(a[i]),float(c[i]),float(d[i]),float(e[i]),float(f[i])]
    #     if ii == 0 and min(temp_list)>=0:
    #         chrom.append(list(Intialchrom[i]))
    #         ObjV.append(list(ii))
    # if len(chrom)>nind:
    #     chrom = chrom[:(nind+1),:]
    
    k = 0   # 循环次数
    while len(chrom) < nind and k < 10*nind:
        Intialchrom = ga.crtpp(nind, lind, varlen)
        chrom_property = list(diction[int(Intialchrom[0][0] - 1), :])
        for i in range(len(Intialchrom) - 1):
            chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
        for j in range(len(Intialchrom[1]) - 1):
            temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
            for i in range(len(Intialchrom) - 1):
                temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
            chrom_property = np.column_stack((chrom_property, temporary))

        a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)
        c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        IntialObjV = a + c + d + e + f
        for i, ii in enumerate(IntialObjV):
            temp_list = [float(a[i]),float(c[i]),float(d[i]),float(e[i]),float(f[i])]
            if ii == 0 and min(temp_list)>=0:
                chrom.append(list(Intialchrom[i]))
                ObjV.append(list(ii))
        k += 1
    
    nind = len(chrom)

    # chrom_property = list(diction[int(chrom[0][0] - 1), :])
    # for i in range(len(chrom) - 1):
    #     chrom_property = np.row_stack((chrom_property, list(diction[int(chrom[i + 1][0] - 1)])))
    # for j in range(len(chrom[1]) - 1):
    #     temporary = list(diction[int(chrom[0][j + 1] - 1)])
    #     for i in range(len(chrom) - 1):
    #         temporary = np.row_stack((temporary, list(diction[int(chrom[i + 1][j + 1] - 1)])))
    #     chrom_property = np.column_stack((chrom_property, temporary))
    
    return chrom, nind, varlen, lind, diction, big_row

    # # 计算相近原则的标准
    # # 初始化该约束的适应度值
    # b = np.zeros((nind, 1))

    # # 遍历所有个体
    # for i in range(nind):
    #     temporary = 0
    #     # 对每个料条进行循环
    #     for m in range(AMN):
    #         # 判断分堆时剩下的料堆是否以依次插入前几个料堆
    #         if m + 1 <= big_row:
    #             # 遍历每一行的每个料堆
    #             for n in range(int(row_number)):
    #                 # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
    #                 # 则该料堆所对应的相近原则的适应度值加一
    #                 Usedifference = abs(
    #                     (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(
    #                         chrom_property[i, (temporary + n) * num_property + 13]) + float(
    #                         chrom_property[i, (temporary + n) * num_property + 16])
    #                      + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(
    #                                 chrom_property[i, (temporary + n) * num_property + 22])) -
    #                     (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(
    #                         chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
    #                         chrom_property[i, (temporary + n + 1) * num_property + 16])
    #                      + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(
    #                                 chrom_property[i, (temporary + n + 1) * num_property + 22]))) * 24 * 30
    #                 if Usedifference < standard_speed_difference:
    #                     b[i] += 1
    #             temporary += row_number + 1
    #         else:
    #             # 遍历每一行的每个料堆
    #             for n in range(int(row_number) - 1):
    #                 # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
    #                 # 则该料堆所对应的相近原则的适应度值加一
    #                 Usedifference = abs(
    #                     (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(
    #                         chrom_property[i, (temporary + n) * num_property + 13]) + float(
    #                         chrom_property[i, (temporary + n) * num_property + 16])
    #                      + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(
    #                                 chrom_property[i, (temporary + n) * num_property + 22])) -
    #                     (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(
    #                         chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
    #                         chrom_property[i, (temporary + n + 1) * num_property + 16])
    #                      + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(
    #                                 chrom_property[i, (temporary + n + 1) * num_property + 22]))) * 24 * 30
    #                 if Usedifference < standard_speed_difference:
    #                     b[i] += 1
    #             temporary += row_number
    # b_standard = max(b)

    # 计算料条利用水平的衡量标准
    # # 初始化该约束的适应度值
    # g = np.zeros((nind, 1))

    # # 遍历所有个体
    # for m in range(nind):
    #     count = 0
    #     # 对每个料条进行循环
    #     for i in range(AMN):
    #         Length = int(chrom_property[m, count * num_property + 4])
    #         # 计算第一料条的长度作为被减数
    #         if i == 0:
    #             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
    #             if i + 1 <= big_row:
    #                 # 遍历每一行的每个料堆
    #                 for j in range(int(row_number)):
    #                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
    #                 count += row_number + 1
    #             else:
    #                 # 遍历每一行的每个料堆
    #                 for j in range(int(row_number) - 1):
    #                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
    #                 count += row_number
    #             pre_row = Length
    #         # 计算其余料条长度
    #         else:
    #             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
    #             if i + 1 <= big_row:
    #                 # 遍历每一行的每个料堆
    #                 for j in range(int(row_number)):
    #                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
    #                 count += row_number + 1
    #             else:
    #                 # 遍历每一行的每个料堆
    #                 for j in range(int(row_number) - 1):
    #                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
    #                 count += row_number
    #             # 令其余料条长度与第一料条长度做差
    #             g[m] += abs(Length-pre_row)
    # g_standard = min(g)


    # # 计算堆取料机运行距离约束的标准值
    # # 初始化该约束的适应度值
    # h = np.zeros((nind, 1))

    # # 遍历所有个体
    # for m in range(nind):
    #     count = 0
    #     # 对每个料条进行循环
    #     for i in range(AMN):
    #         Length = int(chrom_property[m, count * num_property + 4])
    #         # 判断分堆时剩下的料堆是否以依次插入前几个料堆
    #         if i + 1 <= big_row:
    #             # 遍历每一行的每个料堆
    #             # 计算堆取料机假定位置
    #             for j in range(int(row_number)):
    #                 Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
    #             stacker_position = Length/2

    #             # 遍历每一行的每个料堆
    #             # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
    #             temporary = 0
    #             for j in range(int(row_number)+1):
    #                 temporary_1 = temporary
    #                 temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
    #                 temporary_2 = temporary - safe_distance
    #                 ave_temp = (temporary_1+temporary_2)/2
    #                 h[m] += abs(stacker_position - ave_temp)*(int(chrom_property[m, (count + j) * num_property + 9])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 12])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 15])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 18])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 21]))

    #                 # if temporary_2 <= stacker_position:
    #                 #     h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
    #                 # else:
    #                 #     h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
    #                 #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
    #             count += row_number + 1
    #         else:
    #             # 遍历每一行的每个料堆
    #             # 计算堆取料机假定位置
    #             for j in range(int(row_number) - 1):
    #                 Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
    #             stacker_position = Length / 2

    #             # 遍历每一行的每个料堆
    #             # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
    #             temporary = 0
    #             for j in range(int(row_number)):
    #                 temporary_1 = temporary
    #                 temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
    #                 temporary_2 = temporary - safe_distance
    #                 ave_temp = (temporary_1+temporary_2)/2
    #                 h[m] += abs(stacker_position - ave_temp)*(int(chrom_property[m, (count + j) * num_property + 9])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 12])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 15])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 18])
    #                                                                  + int(chrom_property[m, (count + j) * num_property + 21]))

    #             #     if temporary_2 <= stacker_position:
    #             #         h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 12])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 15])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 18])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
    #             #     else:
    #             #         h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 12])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 15])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 18])
    #             #                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
    #             count += row_number
    # h_standard = min(h)

    # # 计算原料对公里数运行距离约束
    # # 初始化该约束的适应度值
    # i = np.zeros((nind, 1))

    # # 遍历所有个体
    # for m in range(nind):
    #     count = 0
    #     # 对每个料条进行循环
    #     for n in range(AMN):
    #         # 判断分堆时剩下的料堆是否以依次插入前几个料堆
    #         if n + 1 <= big_row:
    #             # 遍历每一行的每个料堆
    #             for j in range(int(row_number)):
    #                 # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
    #                 if chrom_property[m, (count + j) * num_property + 8] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 10])
    #                 if chrom_property[m, (count + j) * num_property + 11] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 13])
    #                 if chrom_property[m, (count + j) * num_property + 14] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 16])
    #                 if chrom_property[m, (count + j) * num_property + 17] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 19])
    #                 if chrom_property[m, (count + j) * num_property + 20] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 22])
    #             count += row_number + 1
    #         else:
    #             # 遍历每一行的每个料堆
    #             for j in range(int(row_number) - 1):
    #                 # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
    #                 if chrom_property[m, (count + j) * num_property + 8] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 10])
    #                 if chrom_property[m, (count + j) * num_property + 11] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 13])
    #                 if chrom_property[m, (count + j) * num_property + 14] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 16])
    #                 if chrom_property[m, (count + j) * num_property + 17] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 19])
    #                 if chrom_property[m, (count + j) * num_property + 20] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 22])
    #             count += row_number
    # i_standard = min(i)

    # 输出结果
    # return chrom, b_standard, g_standard, i_standard, h_standard, nind, varlen, lind, diction


def aimfuc(chrom, diction, AMN, weight, row_number,nind, maxleng, standard_speed_difference,big_row, safe_distance, b_standard , g_standard, h_standard, imax, imin, digraph, diction2, num_property):
    # 生成记录适应度值的初始函数
    ObjV= np.zeros((nind, 1))

    # 根据现有种群生成记录种群属性的 list
    chrom_property = list(diction[int(chrom[0][0] - 1), :])
    for i in range(len(chrom) - 1):
        chrom_property = np.row_stack((chrom_property, list(diction[int(chrom[i + 1][0] - 1)])))
    for j in range(len(chrom[1]) - 1):
        temporary = list(diction[int(chrom[0][j + 1] - 1)])
        for i in range(len(chrom) - 1):
            temporary = np.row_stack((temporary, list(diction[int(chrom[i + 1][j + 1] - 1)])))
        chrom_property = np.column_stack((chrom_property, temporary))

    # 分别计算每个约束属性放回的适应度，其中
    # a 料条容量约束； b 相近原则约束
    # c 分散堆放约束； d 固定搭配组合约束
    # e 位置偏好原则； f 产品相容性约束
    # g 料条利用水平约束； h 堆取料机运行距离
    # i 吨公里数运行距离最短约束
    re_pa = []    # 记录所有的约束条件值

    if len(weight)< 10:
        raise Exception ("Not enough parameters")

    if weight[1] != 0:
        a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)
        # re_pa.append(max(a))
        re_pa.append(a)
    else:
        a = 0

    if weight[6] != 0:
        b = bfunc(chrom_property, AMN, row_number, nind, standard_speed_difference, big_row, b_standard, num_property)
        re_pa.append(b)
        # re_pa.append(max(b))
    else:
        b = 0

    if weight[2] != 0:
        c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        re_pa.append(c)
        # re_pa.append(max(c))
    else:
        c = 0

    if weight[3] != 0:
        d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        re_pa.append(d)
        # re_pa.append(max(d))
    else:
        d = 0

    if weight[5] != 0:
        e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        re_pa.append(e)
        # re_pa.append(max(e))
    else:
        e = 0

    if weight[4] != 0:
        f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)
        re_pa.append(f)
        # re_pa.append(max(f))
    else:
        f = 0

    if weight[7] != 0:
        g = gfunc(chrom_property, AMN, row_number,nind, big_row, g_standard, num_property)
        re_pa.append(g)
        # re_pa.append(max(g))
    else:
        g = 0

    if weight[8] != 0:
        h = hfunc(chrom_property, AMN, row_number, nind, big_row, h_standard, safe_distance, num_property,maxleng)
        # h_temp = list(h)
        # while MINVALUE in h_temp:
        #     h_temp.remove(MINVALUE)
        # h_standard = min(h_temp)
        # if h_standard > 0:
        #     h = 2*h_standard/(h+h_standard)*100
        # else:
        #     h = [MINVALUE for _ in range(len(h))]
        re_pa.append(h)
        # re_pa.append(max(h))
    else:
        h = 0

    if weight[9] != 0:
        i = ifunc(chrom_property, AMN, row_number,nind,  big_row, imax, imin, digraph, diction2, num_property)
        re_pa.append(i)
        # re_pa.append(max(i))
    else:
        i = 0

    # 根据权重系数进行计算
    ObjV += weight[1] * a + weight[6] * b + weight[2] * c + weight[3] * d + weight[5] * e + weight[4] * f + weight[7] * g + weight[8] * h + weight[9] * i

    # 将适应度为负的个体的适应度调为0方便之后计算
    for m, i in enumerate(ObjV):
        if i < 0:
            ObjV[m] = 0

    return ObjV,re_pa

#  料条容量约束
def afunc(chrom_property,AMN,row_number,nind,maxleng,big_row, safe_distance, num_property):
    # 初始化该约束的适应度值
    a = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        list_length = []
        # 对每个料条进行循环
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if i + 1 <= big_row:
                # 遍历每一行的每个料堆
                for j in range(int(row_number)):
                    # 计算每个料条的长度
                    list_length.append(int(chrom_property[m, (count + j + 1) * num_property + 4]))
                    Length = Length + int(chrom_property[m, (count + j + 1) * num_property + 4])
                count += row_number + 1
                # 判断每个料条的长度是否超出该料条所能容纳的最大长度
                # if Length >= (maxleng[i] - (1+row_number) * safe_distance):  # 不考虑起始点的安全距离
                if Length > (maxleng[i] - (2+row_number) * safe_distance):  # 考虑起始点的安全距离
                    a[m] = MINVALUE  # 若该料条长度和超过 maxlength 则对应个体的 a 值为 MINVALUE
            else:
                # 遍历每一行的每个料堆
                for j in range(int(row_number) - 1):
                    list_length.append(int(chrom_property[m, (count + j + 1) * num_property + 4]))
                    Length = Length + int(chrom_property[m, (count + j + 1) * num_property + 4])
                count += row_number
                # 判断每个料条的长度是否超出该料条所能容纳的最大长度
                # if Length >= (maxleng[i] - row_number * safe_distance): # 不考虑起始点的安全距离
                if Length > (maxleng[i] - (row_number+1) * safe_distance):  # 考虑起始点的安全距离
                    a[m] = MINVALUE  # 若该料条长度和超过 maxlength 则对应个体的 a 值为 MINVALUE
            # ##-------------------------------
            # if a[m]==0:
            #     label_1 = m
            #     label_list = list_length
            #     print(Length,maxleng[i] - (row_number-1) * safe_distance)
            # ##---------------------------
    return a

def bfunc(chrom_property, AMN, row_number, nind, standard_speed_difference, big_row, b_standard, num_property):
    # 初始化该约束的适应度值
    b = np.zeros((nind, 1))

    # 遍历所有个体
    for i in range(nind):
        temporary = 0
        # 对每个料条进行循环
        for m in range(AMN):
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if m + 1 <= big_row:
                # 遍历每一行的每个料堆
                for n in range(int(row_number)):
                    # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
                    # 则该料堆所对应的相近原则的适应度值加一
                    Usedifference = abs(
                        (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(chrom_property[i, (temporary + n) * num_property + 13]) + float(
                            chrom_property[i, (temporary + n) * num_property + 16])
                         + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(chrom_property[i, (temporary + n) * num_property + 22])) -
                        (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
                            chrom_property[i, (temporary + n + 1) * num_property + 16])
                         + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(chrom_property[i, (temporary + n + 1) * num_property + 22])))*24*30
                    if Usedifference < standard_speed_difference:
                        b[i] += 1
                temporary += row_number + 1
            else:
                # 遍历每一行的每个料堆
                for n in range(int(row_number) - 1):
                    # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
                    # 则该料堆所对应的相近原则的适应度值加一
                    Usedifference = abs(
                        (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(
                            chrom_property[i, (temporary + n) * num_property + 13]) + float(
                            chrom_property[i, (temporary + n) * num_property + 16])
                         + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(
                                    chrom_property[i, (temporary + n) * num_property + 22])) -
                        (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(
                            chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
                            chrom_property[i, (temporary + n + 1) * num_property + 16])
                         + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(
                                    chrom_property[i, (temporary + n + 1) * num_property + 22]))) * 24 * 30
                    if Usedifference < standard_speed_difference:
                        b[i] += 1
                    # if Usedifference > standard_speed_difference:
                    #     b[i] = MINVALUE
                temporary += row_number

    b_standard = row_number*AMN+big_row-AMN
    b = b/b_standard * 100

    # if b_standard != 0:
    #     b = b/b_standard
    return b

def cfunc(chrom_property, AMN, row_number, nind, big_row, num_property):
    # 初始化该约束的适应度值
    c = np.zeros((nind, 1))

    # 遍历所有个体
    for i in range(nind):
        count = 0
        # 对每个料条进行循环
        for m in range(AMN):
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if m+1 <= big_row:
                temp_list = []
                # 遍历每一行的每个料堆
                for n in range(int(row_number)+1):
                    # 对每个料堆判断所处料条是否存在与之重名的料堆，即是否以分散堆放
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 2]
                    if temporaryvalue != '0' and temporaryvalue not in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
                        temp_list.append(temporaryvalue)
                    elif temporaryvalue != '0' and temporaryvalue in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
                        c[i] = MINVALUE
                count += row_number+1
            else:
                temp_list = []
                # 遍历每一行的每个料堆
                for n in range(int(row_number)):
                    # 对每个料堆判断所处料条是否存在与之重名的料堆，即是否以分散堆放
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 2]
                    if temporaryvalue != '0' and temporaryvalue not in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
                        temp_list.append(temporaryvalue)
                    elif temporaryvalue != '0' and temporaryvalue in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
                        c[i] = MINVALUE
                count += row_number
    return c

def dfunc (chrom_property, AMN, row_number, nind, big_row, num_property):
    # 初始化该约束的适应度值
    d = np.zeros((nind, 1))

    # 遍历所有个体
    for i in range(nind):
        count = 0
        # 对每个料条进行循环
        for m in range(AMN):
            temp_list = [[],[]]  # temp_list[0]保存固定搭配值，temp_list[1]是对应的个数
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if m+1 <= big_row:
                # 遍历每一行的每个料堆
                for n in range(int(row_number)+1):
                    # 对每个料堆验证其是否存在固定搭配组合位于该料条
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
                    if temporaryvalue != '0':
                        if temporaryvalue not in temp_list[0]:
                            temp_list[0].append(temporaryvalue)
                            temp_list[1].append(1)
                        else:
                            temp_index = temp_list[0].index(temporaryvalue)
                            temp_list[1][temp_index] += 1
                if 1 in temp_list[1]:
                    d[i] = MINVALUE
                count += row_number+1
            else:
                for n in range(int(row_number)):
                    # 对每个料堆验证其是否存在固定搭配组合位于该料条
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
                    if temporaryvalue != '0':
                        if temporaryvalue not in temp_list[0]:
                            temp_list[0].append(temporaryvalue)
                            temp_list[1].append(1)
                        else:
                            temp_index = temp_list[0].index(temporaryvalue)
                            temp_list[1][temp_index] += 1
                if 1 in temp_list[1]:
                    d[i] = MINVALUE
                count += row_number

    # # 遍历所有个体
    # for i in range(nind):
    #     count = 0
    #     # 对每个料条进行循环
    #     for m in range(AMN):
    #         # 判断分堆时剩下的料堆是否以依次插入前几个料堆
    #         if m + 1 <= big_row:
    #             # 遍历每一行的每个料堆
    #             for n in range(int(row_number) - 1):
    #                 # 对每个料堆判断其是否被分堆，若被分堆观察其被分堆的组合是否位于相邻位置
    #                 if chrom_property[i, (count + n + 1) *num_property + 23] != 0 and chrom_property[i, (count + n + 1) *num_property + 23] != 0.0:
    #                     if chrom_property[i, (count + n) * num_property + 23] != chrom_property[i, (count + n + 1) *
    #                                                                                                num_property + 23] \
    #                             and chrom_property[i, (count + n + 1) * num_property + 23] != chrom_property[i,
    #                                                                                                          (
    #                                                                                                                  count + n + 2) * num_property + 23]:
    #                         d[i] = MINVALUE
    #             count += row_number + 1
    #         else:
    #             # 遍历每一行的每个料堆
    #             for n in range(int(row_number) - 2):
    #                 # 对每个料堆判断其是否被分堆，若被分堆观察其被分堆的组合是否位于相邻位置
    #                 if chrom_property[i, (count + n + 1) *num_property + 23] != 0 and chrom_property[i, (count + n + 1) *num_property + 23] != 0.0:
    #                     if chrom_property[i, (count + n) * num_property + 23] != chrom_property[i, (count + n + 1) *
    #                                                                                                num_property + 23] \
    #                             and chrom_property[i, (count + n + 1) * num_property + 23] != chrom_property[i,
    #                                                                                                          (
    #                                                                                                                  count + n + 2) * num_property + 23]:
    #                         d[i] = MINVALUE
    #             count += row_number
    # print(max(d))
    return d

def efunc (chrom_property, AMN, row_number, nind, big_row, num_property):
    # 初始化该约束的适应度值
    e = np.zeros((nind, 1))

    # 遍历所有个体
    for i in range(nind):
        count = 0
        # 对每个料条进行循环
        for m in range(AMN):
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if m + 1 <= big_row:
                # 遍历每一行的每个料堆
                for n in range(int(row_number)+1):
                    # 判断该料堆是否存在偏好位置并注明偏好位置所在料条
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 6]
                    if int(temporaryvalue) != 0 and int(temporaryvalue) != m+1:
                        e[i] = MINVALUE
                count += row_number+1
            else:
                # 遍历每一行的每个料堆
                for n in range(int(row_number)):
                    # 判断该料堆是否存在偏好位置并注明偏好位置所在料条
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 6]
                    if int(temporaryvalue) != 0 and int(temporaryvalue) != m+1:
                        e[i] = MINVALUE
                count += row_number
    return e

def ffunc(chrom_property, AMN, row_number, nind, big_row, num_property):
    # 初始化该约束的适应度值
    f = np.zeros((nind, 1))

    # 遍历所有个体
    for i in range(nind):
        count = 0
        # 对每个料条进行循环
        for m in range(AMN):
                # 判断分堆时剩下的料堆是否以依次插入前几个料堆
                if m + 1 <= big_row:
                    # 遍历每一行的每个料堆
                    for n in range(int(row_number)):
                        temporaryvalue = chrom_property[i, (count + n) * num_property + 7]
                        # 判断该料堆是否存在不相容的组合且位于同一料条
                        if int(temporaryvalue) != 0 and temporaryvalue == chrom_property[i, (count + n + 1) * num_property + 7]:
                            f[i] = MINVALUE
                    count += row_number+1
                else:
                    # 遍历每一行的每个料堆
                    for n in range(int(row_number)-1):
                        temporaryvalue = chrom_property[i, (count + n) * num_property + 7]
                        # 判断该料堆是否存在不相容的组合且位于同一料条
                        if int(temporaryvalue) != 0 and temporaryvalue == chrom_property[i, (count + n + 1) * num_property + 7]:
                            f[i] = MINVALUE
                    count += row_number
    return f

def gfunc(chrom_property, AMN, row_number,nind, big_row, g_standard, num_property):
    # 初始化该约束的适应度值
    g = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        length_list = []   # 记录个料条的长度
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if i + 1 <= big_row:
                # 遍历每一行的每个料堆
                for j in range(int(row_number)):
                    Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                for j in range(int(row_number) - 1):
                    Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])
                count += row_number
            length_list.append(Length)

        # 令其余料条长度与平均料条长度做差
        ave_row = sum(length_list)/AMN
        for i in range(len(length_list)):
            length_list[i] = abs(length_list[i]-ave_row)
        g[m] = sum(length_list)

            # 计算第一料条的长度作为被减数
            # if i == 0:
            #     # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            #     if i + 1 <= big_row:
            #         # 遍历每一行的每个料堆
            #         for j in range(int(row_number)):
            #             Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])
            #         count += row_number + 1
            #     else:
            #         # 遍历每一行的每个料堆
            #         for j in range(int(row_number) - 1):
            #             Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])
            #         count += row_number
            #     pre_row = Length
            # # 计算其余料条长度
            # else:
            #     # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            #     if i + 1 <= big_row:
            #         # 遍历每一行的每个料堆
            #         for j in range(int(row_number)):
            #             Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])
            #         count += row_number + 1
            #     else:
            #         # 遍历每一行的每个料堆
            #         for j in range(int(row_number) - 1):
            #             Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])
            #         count += row_number
            #     # 令其余料条长度与第一料条长度做差
            #     g[m] += abs(Length-pre_row)
    g = 100*(g_standard-g)/g_standard
    return g

def cal_g(diction, maxleng, safe_distance,AMN):
    # 计算一个接近最差情况的g作为g_standard
    # diction 是dataframe格式的物料信息表，maxleng是list格式的料条最长长度,safe_distance是安全距离
    item_len = list(diction.iloc[:,4])
    item_len = sorted(item_len,reverse=True)   # 料堆长短排序
    slice_len = sorted(maxleng,reverse=True)   # 料条长短排序

    # 尽量排满最长的几排，并计算出g
    i = 0
    length_list = [] # 记录每一排的长度
    for slen in slice_len:
        if i>=len(item_len):
            break
        rest = slen-safe_distance   # 考虑开始的安全距离,计算剩余距离
        Length = item_len[i]
        i += 1
        if i>=len(item_len):
            break
        while item_len[i]+safe_distance < rest:
            Length += item_len[i]
            rest = rest - item_len[i] - safe_distance
            i += 1
            if i>=len(item_len):
                break

        # 如果最长的料条放不进最长的料条，则无解
        if i == 1:
            pirnt("Please reset parameters.")
            raise Exception ("Pile exceeds length!")

        length_list.append(Length)

    # 令其余料条长度与平均料条长度做差
    ave_row = sum(length_list)/AMN
    for i in range(len(length_list)):
        length_list[i] =abs(length_list[i]-ave_row) 
    g_standard = sum(length_list)

    return g_standard

def hfunc(chrom_property, AMN, row_number, nind, big_row, h_standard, safe_distance, num_property,maxleng):
    # 初始化该约束的适应度值
    h = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])+2*safe_distance
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if i + 1 <= big_row:
                # 遍历每一行的每个料堆
                # 计算堆取料机假定位置
                for j in range(int(row_number)):
                    Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])+safe_distance
                stacker_position = Length/2

                # 遍历每一行的每个料堆
                temporary = safe_distance
                for j in range(int(row_number)+1):
                    # 计算料堆的中点位置到堆取料机的距离
                    item_length = int(chrom_property[m, (count + j) * num_property + 4])
                    rea_position = temporary + item_length/2
                    temporary += item_length + safe_distance
                    if temporary > maxleng[i]:
                        h[m] = MINVALUE
                        break
                    elif h[m] != MINVALUE:
                        h[m] += abs(stacker_position - rea_position)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                    # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
                    # temporary_1 = temporary
                    # temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
                    # temporary_2 = temporary - safe_distance
                    # if temporary_2 <= stacker_position:
                    #     h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
                    # else:
                    #     h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                # 计算堆取料机假定位置
                for j in range(int(row_number) - 1):
                    Length = Length + int(chrom_property[m, (count + j+1) * num_property + 4])+safe_distance
                stacker_position = Length / 2

                # 遍历每一行的每个料堆
                temporary = 0
                for j in range(int(row_number)):
                    # 计算料堆的中点位置到堆取料机的距离
                    item_length = int(chrom_property[m, (count + j) * num_property + 4])
                    rea_position = temporary + item_length/2
                    temporary += item_length + safe_distance
                    if temporary>maxleng[i]:
                        h[m] = MINVALUE
                        break
                    elif h[m] != MINVALUE:
                        h[m] += abs(stacker_position - rea_position)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))

                    # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
                    # temporary_1 = temporary
                    # temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
                    # temporary_2 = temporary - safe_distance
                    # if temporary_2 <= stacker_position:
                    #     h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
                    # else:
                    #     h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
                    #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
                count += row_number
    # h = 2*h_standard/(h+h_standard)*100
    # for ind in range(len(h)):
    #     if h[ind]>0:
    temp = 5*h/h_standard
    h = 100 / (1 + np.exp(-(temp)))
  
    return h

def cal_h(maxleng, in_material):
    max_distance = max(maxleng)/2
    # 计算频率
    h_ord = [9,12,15,18,21]   # 频率的位置
    col_name = list(in_material.columns)
    name= []
    for item in h_ord:
        name.append(col_name[item])
    name.append("sum")
    new_df = in_material.iloc[:,h_ord].copy()
    new_df[name[-1]] = new_df[name[0]]
    for ind in range(1,len(name)):
        new_df[name[-1]] += new_df[name[ind]]

    sum_list = new_df[name[-1]]
    max_fre = max(sum_list)

    return max_distance*max_fre

def ifunc ( chrom_property, AMN, row_number,nind,  big_row, maxi, mini, digraph, diction2, num_property):
    # 初始化该约束的适应度值
    i = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        for n in range(AMN):
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if n + 1 <= big_row:
                # 遍历每一行的每个料堆
                for j in range(int(row_number)+1):
                    # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
                    if chrom_property[m, (count + j) * num_property + 8] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
                            chrom_property[m, (count + j) * num_property + 10])
                    if chrom_property[m, (count + j) * num_property + 11] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
                            chrom_property[m, (count + j) * num_property + 13])
                    if chrom_property[m, (count + j) * num_property + 14] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
                            chrom_property[m, (count + j) * num_property + 16])
                    if chrom_property[m, (count + j) * num_property + 17] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
                            chrom_property[m, (count + j) * num_property + 19])
                    if chrom_property[m, (count + j) * num_property + 20] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
                            chrom_property[m, (count + j) * num_property + 22])
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                for j in range(int(row_number)):
                    # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
                    if chrom_property[m, (count + j) * num_property + 8] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
                            chrom_property[m, (count + j) * num_property + 10])
                    if chrom_property[m, (count + j) * num_property + 11] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
                            chrom_property[m, (count + j) * num_property + 13])
                    if chrom_property[m, (count + j) * num_property + 14] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
                            chrom_property[m, (count + j) * num_property + 16])
                    if chrom_property[m, (count + j) * num_property + 17] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
                            chrom_property[m, (count + j) * num_property + 19])
                    if chrom_property[m, (count + j) * num_property + 20] != '0':
                        i[m] += digraph[n][
                                    diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
                            chrom_property[m, (count + j) * num_property + 22])
                count += row_number


    # # 遍历所有个体
    # for m in range(nind):
    #     count = 0
    #     # 对每个料条进行循环
    #     for n in range(AMN):
    #         # 判断分堆时剩下的料堆是否以依次插入前几个料堆
    #         if n + 1 <= big_row:
    #             # 遍历每一行的每个料堆
    #             for j in range(int(row_number)+1):
    #                 # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
    #                 if chrom_property[m, (count + j) * num_property + 8] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 10])
    #                 if chrom_property[m, (count + j) * num_property + 11] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 13])
    #                 if chrom_property[m, (count + j) * num_property + 14] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 16])
    #                 if chrom_property[m, (count + j) * num_property + 17] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 19])
    #                 if chrom_property[m, (count + j) * num_property + 20] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 22])
    #             count += row_number + 1
    #         else:
    #             # 遍历每一行的每个料堆
    #             for j in range(int(row_number)):
    #                 # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
    #                 if chrom_property[m, (count + j) * num_property + 8] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 10])
    #                 if chrom_property[m, (count + j) * num_property + 11] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 13])
    #                 if chrom_property[m, (count + j) * num_property + 14] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 16])
    #                 if chrom_property[m, (count + j) * num_property + 17] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 19])
    #                 if chrom_property[m, (count + j) * num_property + 20] != '0':
    #                     i[m] += digraph[diction2[str(n + 1)]][
    #                                 diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
    #                         chrom_property[m, (count + j) * num_property + 22])
    #             count += row_number
    # i = 2*i_standard/(i+i_standard)*100
    # i = 100*(maxi-i)/(maxi-mini)
    # for ind in range(len(i)):
    #     if i[ind]>0:
    temp = 4*(maxi-i)/(maxi-mini)
    i = 100*( 2/ (1 + np.exp(-(temp)))-1)
 
    # print(min(i),max(i))
    return i

def cal_i(dmin,udic,AMN,slicelist,userlist,max_ves,min_ves):
    # 计算得到i_min,i_max
    # (最近矩阵，最远矩阵，料条和用户对应数字编号，将插入的原始数据表，
    # 料条数, 矩阵中的所有料条序号,矩阵中的所有用户序号,最大、最小频率值)

    # i_min = dmin[udic[slicelist[0]]][udic[userlist[0]]]
    # i_max = dmin[udic[slicelist[0]]][udic[userlist[0]]]
    i_min = dmin[0][udic[userlist[0]]]
    i_max = dmin[0][udic[userlist[0]]]
    for i in range(len(slicelist)):
        for j in range(len(userlist)):
            if dmin[i][udic[userlist[j]]] < 0:
                print("No path exists.From slice ",slicelist[i]," to user ",userlist[j])
                raise Exception ("No path exists!")
            dis = dmin[i][udic[userlist[j]]]
            i_min = min(i_min,dis)
            dis = dmin[i][udic[userlist[j]]]
            i_max = max(i_max,dis)

    # i_min = dmin[udic[slicelist[0]]][udic[userlist[0]]]*velist[0]
    # i_max = dmin[udic[slicelist[0]]][udic[userlist[0]]]*velist[0]
    # for i in range(len(slicelist)):
    #     for j in range(len(userlist)):
    #         if dmin[udic[slicelist[i]]][udic[userlist[j]]] < 0:
    #             print("No path exists.From slice ",slicelist[i]," to user ",userlist[j])
    #             raise Exception ("No path exists!")
    #         dis = dmin[udic[slicelist[i]]][udic[userlist[j]]]
    #         i_min = min(i_min,dis*velist[j])
    #         dis = dmin[udic[slicelist[i]]][udic[userlist[j]]]
    #         i_max = max(i_max,dis*velist[j])

    # 对应五个用户的值
    i_min *= min_ves*5
    i_max *= max_ves*5

    return i_max,i_min

def users_velocity(indf,AMN):
    # 计算距离的最近、最远值
    c_user = [8,11,14,17,20]

    # indf.iloc[:,[8,11,14,17,20]]
    userdf = indf.iloc[:,[c_user[0],c_user[0]+2]]
    USERS = "用户"
    VES = "速率"
    userdf.columns = [USERS,VES]
    max_ves = max(userdf[VES])
    min_ves = min(userdf[VES])
    for i in range(1,len(c_user)):
        temp = indf.iloc[:,[c_user[i],c_user[i]+2]]
        temp.columns = userdf.columns
        userdf = pd.concat([userdf,temp],axis=0,join='outer',ignore_index=True)
    # print(userdf)

    userdf.drop_duplicates(subset=USERS,keep='first', inplace=True)  # 去重
    userdf = userdf.loc[(userdf[USERS]!=0)&(userdf[USERS]!='0'),:]

    # 记录所有用户
    # slicelist = [str(x) for x in range(1,AMN+1)]  # 料条
    userlist = list(userdf.iloc[:,0])
    # velist = list(userdf.iloc[:,1])    

    return userlist,max_ves,min_ves


# def filter_function(Intialchrom, diction, AMN, row_number, nind, maxleng, standard_speed_difference,big_row, lind, varlen, safe_distance, digraph, diction2, num_property):

#     # 根据初始种群和 diction 将每个位置料堆对应的属性整合成一个list
#     chrom_property = np.array(diction[int(Intialchrom[0][0] - 1), :])
#     for i in range(len(Intialchrom) - 1):
#         chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
#     for j in range(len(Intialchrom[1]) - 1):
#         temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
#         for i in range(len(Intialchrom) - 1):
#             temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
#         chrom_property = np.column_stack((chrom_property, temporary))

#     # 计算料堆容量约束
#     a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)

#     # 检查是否存在某个料堆过长导致现有的平均堆数放到每个料条的情况难以实现的情况，若出现则将其看作两个料堆
#     while max(a) == -100:
#         temporary = diction[:, 4]
#         temporary = np.argmax(temporary)
#         diction[temporary][3] = int(diction[temporary][3] / 2)
#         diction[temporary][4] = int(diction[temporary][4] / 2)
#         diction[temporary][9] = int(diction[temporary][9] / 2)
#         diction[temporary][10] = int(diction[temporary][10] / 2)
#         diction[temporary][12] = int(diction[temporary][12] / 2)
#         diction[temporary][13] = int(diction[temporary][13] / 2)
#         diction[temporary][15] = int(diction[temporary][15] / 2)
#         diction[temporary][16] = int(diction[temporary][16] / 2)
#         diction[temporary][18] = int(diction[temporary][18] / 2)
#         diction[temporary][19] = int(diction[temporary][19] / 2)
#         diction[temporary][21] = int(diction[temporary][21] / 2)
#         diction[temporary][22] = int(diction[temporary][22] / 2)
#         diction[temporary][-1] = len(diction)
#         diction = np.row_stack((diction, diction[temporary]))
#         lind += 1
#         big_row += 1
#         varlen = len(diction)
#         if big_row == AMN:
#             row_number += 1
#             big_row = 0
#         Intialchrom = ga.crtpp(nind, lind, varlen)
#         chrom_property = np.array(diction[int(Intialchrom[0][0] - 1), :])
#         for i in range(len(Intialchrom) - 1):
#             chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
#         for j in range(len(Intialchrom[1]) - 1):
#             temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
#             for i in range(len(Intialchrom) - 1):
#                 temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
#             chrom_property = np.column_stack((chrom_property, temporary))
#         a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)

#     # 计算分散搭配约束的适应度
#     c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)

#     # 计算固定搭配组合的适应度
#     d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)

#     # 计算是否又因为固定搭配组合放到一个料堆导致的放不下的情况， 若有将固定搭配组合中最长的料堆拆成两个料堆再进行计算
#     while max(d) == -100:
#         temporary = diction[:, 4]
#         for i,ii in enumerate(temporary):
#             if diction[i][5] == 0:
#                 temporary[i] = 0
#         temporary = np.argmax(temporary)
#         diction[temporary][3] = int(diction[temporary][3] / 2)
#         diction[temporary][4] = int(diction[temporary][4] / 2)
#         diction[temporary][9] = int(diction[temporary][9] / 2)
#         diction[temporary][10] = int(diction[temporary][10] / 2)
#         diction[temporary][12] = int(diction[temporary][12] / 2)
#         diction[temporary][13] = int(diction[temporary][13] / 2)
#         diction[temporary][15] = int(diction[temporary][15] / 2)
#         diction[temporary][16] = int(diction[temporary][16] / 2)
#         diction[temporary][18] = int(diction[temporary][18] / 2)
#         diction[temporary][19] = int(diction[temporary][19] / 2)
#         diction[temporary][21] = int(diction[temporary][21] / 2)
#         diction[temporary][22] = int(diction[temporary][22] / 2)
#         diction[temporary][-1] = len(diction)
#         diction = np.row_stack((diction, diction[temporary]))
#         lind += 1
#         big_row += 1
#         varlen = len(diction)
#         if big_row == AMN:
#             row_number += 1
#             big_row = 0
#         Intialchrom = ga.crtpp(nind, lind, varlen)
#         chrom_property = np.array(diction[int(Intialchrom[0][0] - 1), :])
#         for i in range(len(Intialchrom) - 1):
#             chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
#         for j in range(len(Intialchrom[1]) - 1):
#             temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
#             for i in range(len(Intialchrom) - 1):
#                 temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
#             chrom_property = np.column_stack((chrom_property, temporary))
#         d = dfunc(chrom_property, diction, AMN, row_number, nind, big_row, num_property)

#     # 计算位置偏好约束和产品相容性约束的适应度
#     e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#     f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)

#     # 合并五条适应度的结果，选取同时满足五条适应度的方案放入新的初始种群
#     IntialObjV = a + c + d + e + f
#     chrom, ObjV = [], []
#     for i,ii in enumerate(IntialObjV):
#         if ii == 0:
#             chrom.append(list(Intialchrom[i]))
#             ObjV.append(list(IntialObjV[i]))

#     # 重复生成初始种群并计算五条适应度直到新的初始种群的大小不小于原来的大小
#     while len(chrom) < nind:
#         Intialchrom = ga.crtpp(nind, lind, varlen)
#         chrom_property = list(diction[int(Intialchrom[0][0] - 1), :])
#         for i in range(len(Intialchrom) - 1):
#             chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
#         for j in range(len(Intialchrom[1]) - 1):
#             temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
#             for i in range(len(Intialchrom) - 1):
#                 temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
#             chrom_property = np.column_stack((chrom_property, temporary))

#         a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)
#         c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#         d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#         e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#         f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#         IntialObjV = a + c + d + e + f
#         for i, ii in enumerate(IntialObjV):
#             if ii == 0:
#                 chrom.append(list(Intialchrom[i]))
#                 ObjV.append(list(ii))
#     nind = len(chrom)
#     chrom_property = list(diction[int(chrom[0][0] - 1), :])
#     for i in range(len(chrom) - 1):
#         chrom_property = np.row_stack((chrom_property, list(diction[int(chrom[i + 1][0] - 1)])))
#     for j in range(len(chrom[1]) - 1):
#         temporary = list(diction[int(chrom[0][j + 1] - 1)])
#         for i in range(len(chrom) - 1):
#             temporary = np.row_stack((temporary, list(diction[int(chrom[i + 1][j + 1] - 1)])))
#         chrom_property = np.column_stack((chrom_property, temporary))

#     # 计算相近原则的标准
#     # 初始化该约束的适应度值
#     b = np.zeros((nind, 1))

#     # 遍历所有个体
#     for i in range(nind):
#         temporary = 0
#         # 对每个料条进行循环
#         for m in range(AMN):
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if m + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number)):
#                     # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
#                     # 则该料堆所对应的相近原则的适应度值加一
#                     Usedifference = abs(
#                         (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(
#                             chrom_property[i, (temporary + n) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(
#                                     chrom_property[i, (temporary + n) * num_property + 22])) -
#                         (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(
#                             chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n + 1) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(
#                                     chrom_property[i, (temporary + n + 1) * num_property + 22]))) * 24 * 30
#                     if Usedifference < standard_speed_difference:
#                         b[i] += 1
#                 temporary += row_number + 1
#             else:
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number) - 1):
#                     # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
#                     # 则该料堆所对应的相近原则的适应度值加一
#                     Usedifference = abs(
#                         (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(
#                             chrom_property[i, (temporary + n) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(
#                                     chrom_property[i, (temporary + n) * num_property + 22])) -
#                         (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(
#                             chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n + 1) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(
#                                     chrom_property[i, (temporary + n + 1) * num_property + 22]))) * 24 * 30
#                     if Usedifference < standard_speed_difference:
#                         b[i] += 1
#                 temporary += row_number
#     b_standard = max(b)

#     # 计算料条利用水平的衡量标准
#     # 初始化该约束的适应度值
#     g = np.zeros((nind, 1))

#     # 遍历所有个体
#     for m in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for i in range(AMN):
#             Length = int(chrom_property[m, count * num_property + 4])
#             # 计算第一料条的长度作为被减数
#             if i == 0:
#                 # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#                 if i + 1 <= big_row:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number)):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number + 1
#                 else:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number) - 1):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number
#                 pre_row = Length
#             # 计算其余料条长度
#             else:
#                 # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#                 if i + 1 <= big_row:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number)):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number + 1
#                 else:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number) - 1):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number
#                 # 令其余料条长度与第一料条长度做差
#                 g[m] += abs(Length-pre_row)
#     g_standard = min(g)

#     # 计算堆取料机运行距离约束的标准值
#     # 初始化该约束的适应度值
#     h = np.zeros((nind, 1))

#     # 遍历所有个体
#     for m in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for i in range(AMN):
#             Length = int(chrom_property[m, count * num_property + 4])
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if i + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 # 计算堆取料机假定位置
#                 for j in range(int(row_number)):
#                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                 stacker_position = Length/2

#                 # 遍历每一行的每个料堆
#                 # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
#                 temporary = 0
#                 for j in range(int(row_number)+1):
#                     temporary_1 = temporary
#                     temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
#                     temporary_2 = temporary - safe_distance
#                     ave_temp = (temporary_1+temporary_2)/2
#                     h[m] += abs(stacker_position - ave_temp)*(int(chrom_property[m, (count + j) * num_property + 9])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
                    
#                 #     if temporary_2 <= stacker_position:
#                 #         h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
#                 #     else:
#                 #         h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
#                 count += row_number + 1
#             else:
#                 # 遍历每一行的每个料堆
#                 # 计算堆取料机假定位置
#                 for j in range(int(row_number) - 1):
#                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                 stacker_position = Length / 2

#                 # 遍历每一行的每个料堆
#                 # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
#                 temporary = 0
#                 for j in range(int(row_number)):
#                     temporary_1 = temporary
#                     temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
#                     temporary_2 = temporary - safe_distance
#                     ave_temp = (temporary_1+temporary_2)/2
#                     h[m] += abs(stacker_position - ave_temp)*(int(chrom_property[m, (count + j) * num_property + 9])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
                    
#                 #     if temporary_2 <= stacker_position:
#                 #         h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
#                 #     else:
#                 #         h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                 #                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
#                 count += row_number
#     h_standard = min(h)

#     # 计算原料对公里数运行距离约束
#     # 初始化该约束的适应度值
#     i = np.zeros((nind, 1))

#     # 遍历所有个体
#     for m in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for n in range(AMN):
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if n + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 for j in range(int(row_number)):
#                     # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
#                     if chrom_property[m, (count + j) * num_property + 8] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
#                             chrom_property[m, (count + j) * num_property + 10])
#                     if chrom_property[m, (count + j) * num_property + 11] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
#                             chrom_property[m, (count + j) * num_property + 13])
#                     if chrom_property[m, (count + j) * num_property + 14] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
#                             chrom_property[m, (count + j) * num_property + 16])
#                     if chrom_property[m, (count + j) * num_property + 17] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
#                             chrom_property[m, (count + j) * num_property + 19])
#                     if chrom_property[m, (count + j) * num_property + 20] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
#                             chrom_property[m, (count + j) * num_property + 22])
#                 count += row_number + 1
#             else:
#                 # 遍历每一行的每个料堆
#                 for j in range(int(row_number) - 1):
#                     # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
#                     if chrom_property[m, (count + j) * num_property + 8] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
#                             chrom_property[m, (count + j) * num_property + 10])
#                     if chrom_property[m, (count + j) * num_property + 11] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
#                             chrom_property[m, (count + j) * num_property + 13])
#                     if chrom_property[m, (count + j) * num_property + 14] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
#                             chrom_property[m, (count + j) * num_property + 16])
#                     if chrom_property[m, (count + j) * num_property + 17] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
#                             chrom_property[m, (count + j) * num_property + 19])
#                     if chrom_property[m, (count + j) * num_property + 20] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
#                             chrom_property[m, (count + j) * num_property + 22])
#                 count += row_number
#     i_standard = min(i)

#     # 输出结果
#     return chrom, b_standard, g_standard, i_standard, h_standard, nind, varlen, lind, diction


# def aimfuc(chrom, diction, AMN, weight, row_number,nind, maxleng, standard_speed_difference,big_row, safe_distance, b_standard , g_standard, h_standard,  i_standard, digraph, diction2, num_property):
#     # 生成记录适应度值的初始函数
#     ObjV= np.zeros((nind, 1))

#     # 根据现有种群生成记录种群属性的 list
#     chrom_property = list(diction[int(chrom[0][0] - 1), :])
#     for i in range(len(chrom) - 1):
#         chrom_property = np.row_stack((chrom_property, list(diction[int(chrom[i + 1][0] - 1)])))
#     for j in range(len(chrom[1]) - 1):
#         temporary = list(diction[int(chrom[0][j + 1] - 1)])
#         for i in range(len(chrom) - 1):
#             temporary = np.row_stack((temporary, list(diction[int(chrom[i + 1][j + 1] - 1)])))
#         chrom_property = np.column_stack((chrom_property, temporary))

#     # 分别计算每个约束属性放回的适应度，其中
#     # a 料条容量约束； b 相近原则约束
#     # c 分散堆放约束； d 固定搭配组合约束
#     # e 位置偏好原则； f 产品相容性约束
#     # g 料条利用水平约束； h 堆取料机运行距离
#     # i 吨公里数运行距离最短约束
#     if weight[1] != 0:
#         a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)
#     else:
#         a = 0

#     if weight[2] != 0:
#         b = bfunc(chrom_property, AMN, row_number, nind, standard_speed_difference, big_row, b_standard, num_property)
#     else:
#         b = 0

#     if weight[3] != 0:
#         c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#     else:
#         c = 0

#     if weight[4] != 0:
#         d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#     else:
#         d = 0

#     if weight[5] != 0:
#         e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#     else:
#         e = 0

#     if weight[6] != 0:
#         f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)
#     else:
#         f = 0

#     if weight[7] != 0:
#         g = gfunc(chrom_property, AMN, row_number,nind, big_row, g_standard, num_property)
#     else:
#         g = 0

#     if weight[8] != 0:
#         h = hfunc(chrom_property, AMN, row_number, nind, big_row, h_standard, safe_distance, num_property)
#     else:
#         h = 0

#     if weight[9] != 0:
#         i = ifunc(chrom_property, AMN, row_number,nind,  big_row, i_standard, digraph, diction2, num_property)
#     else:
#         i = 0

#     # 根据权重系数进行计算 + weight[8] * h
#     ObjV += weight[1] * a + weight[6] * b + weight[2] * c + weight[3] * d + weight[5] * e + weight[4] * f + weight[7] * g + weight[8] * h + weight[9] * i

#     # 将适应度为负的个体的适应度调为0方便之后计算
#     for m, i in enumerate(ObjV):
#         if i < 0:
#             ObjV[m] = 0
#     return ObjV

# def afunc(chrom_property,AMN,row_number,nind,maxleng,big_row, safe_distance, num_property):
#     # 初始化该约束的适应度值
#     a = np.zeros((nind, 1))

#     # 遍历所有个体
#     for m in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for i in range(AMN):
#             Length = int(chrom_property[m, count * num_property + 4])
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if i + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 for j in range(int(row_number)):
#                     # 计算每个料条的长度
#                     Length = Length + int(chrom_property[m, (count + j + 1) * num_property + 4])
#                 count += row_number + 1
#                 # 判断每个料条的长度是否超出该料条所能容纳的最大长度
#                 if Length >= (maxleng[i] - row_number * safe_distance):
#                         a[m] = -100  # 若该料条长度和超过 maxlength 则对应个体的 a 值为 -100
#             else:
#                 # 遍历每一行的每个料堆
#                 for j in range(int(row_number) - 1):
#                     Length = Length + int(chrom_property[m, (count + j + 1) * num_property + 4])
#                 count += row_number
#                 # 判断每个料条的长度是否超出该料条所能容纳的最大长度
#                 if Length >= (maxleng[i] - (row_number-1) * safe_distance):
#                         a[m] = -100  # 若该料条长度和超过 maxlength 则对应个体的 a 值为 -100
#     return a

# def bfunc(chrom_property, AMN, row_number, nind, standard_speed_difference, big_row, b_standard, num_property):
#     # 初始化该约束的适应度值
#     b = np.zeros((nind, 1))

#     # 遍历所有个体
#     for i in range(nind):
#         temporary = 0
#         # 对每个料条进行循环
#         for m in range(AMN):
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if m + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number)):
#                     # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
#                     # 则该料堆所对应的相近原则的适应度值加一
#                     Usedifference = abs(
#                         (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(chrom_property[i, (temporary + n) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(chrom_property[i, (temporary + n) * num_property + 22])) -
#                         (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n + 1) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(chrom_property[i, (temporary + n + 1) * num_property + 22])))*24*30
#                     if Usedifference < standard_speed_difference:
#                         b[i] += 1
#                 temporary += row_number + 1
#             else:
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number) - 1):
#                     # 对每个料堆判断其每月的总消耗量并于相邻的下一料堆的该值做差，若插值小于 standard_speed_difference
#                     # 则该料堆所对应的相近原则的适应度值加一
#                     Usedifference = abs(
#                         (float(chrom_property[i, (temporary + n) * num_property + 10]) + float(
#                             chrom_property[i, (temporary + n) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n) * num_property + 19]) + float(
#                                     chrom_property[i, (temporary + n) * num_property + 22])) -
#                         (float(chrom_property[i, (temporary + n + 1) * num_property + 10]) + float(
#                             chrom_property[i, (temporary + n + 1) * num_property + 13]) + float(
#                             chrom_property[i, (temporary + n + 1) * num_property + 16])
#                          + float(chrom_property[i, (temporary + n + 1) * num_property + 19]) + float(
#                                     chrom_property[i, (temporary + n + 1) * num_property + 22]))) * 24 * 30
#                     if Usedifference < standard_speed_difference:
#                         b[i] += 1
#                 temporary += row_number
#     if b_standard != 0:
#         b = b/b_standard*100
#     return b

# def cfunc(chrom_property, AMN, row_number, nind, big_row, num_property):
#     # 初始化该约束的适应度值
#     c = np.zeros((nind, 1))

#     # 遍历所有个体
#     for i in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for m in range(AMN):
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if m+1 <= big_row:
#                 temp_list = []
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number)+1):
#                     # 对每个料堆判断所处料条是否存在与之重名的料堆，即是否以分散堆放
#                     temporaryvalue = chrom_property[i, (count + n) * num_property + 2]
#                     if temporaryvalue != '0' and temporaryvalue not in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
#                         temp_list.append(temporaryvalue)
#                     elif temporaryvalue != '0' and temporaryvalue in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
#                         c[i] = -100
#                 count += row_number+1
#             else:
#                 temp_list = []
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number)):
#                     # 对每个料堆判断所处料条是否存在与之重名的料堆，即是否以分散堆放
#                     temporaryvalue = chrom_property[i, (count + n) * num_property + 2]
#                     if temporaryvalue != '0' and temporaryvalue not in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
#                         temp_list.append(temporaryvalue)
#                     elif temporaryvalue != '0' and temporaryvalue in temp_list and chrom_property[i, (count + n) * num_property + 23] == '0.0':
#                         c[i] = -100
#                 count += row_number
#     return c

# def dfunc (chrom_property, AMN, row_number, nind, big_row, num_property):
#     # 初始化该约束的适应度值
#     d = np.zeros((nind, 1))

#     # 遍历所有个体
#     for i in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for m in range(AMN):
#             temp_list = [[],[]]  # temp_list[0]保存固定搭配值，temp_list[1]是对应的个数
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if m+1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number)+1):
#                     # 对每个料堆验证其是否存在固定搭配组合位于该料条
#                     temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
#                     if temporaryvalue != '0':
#                         if temporaryvalue not in temp_list[0]:
#                             temp_list[0].append(temporaryvalue)
#                             temp_list[1].append(1)
#                         else:
#                             temp_index = temp_list[0].index(temporaryvalue)
#                             temp_list[1][temp_index] += 1
#                 if 1 in temp_list[1]:
#                     d[i] = -100
#                 count += row_number+1
#             else:
#                 for n in range(int(row_number)):
#                     # 对每个料堆验证其是否存在固定搭配组合位于该料条
#                     temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
#                     if temporaryvalue != '0':
#                         if temporaryvalue not in temp_list[0]:
#                             temp_list[0].append(temporaryvalue)
#                             temp_list[1].append(1)
#                         else:
#                             temp_index = temp_list[0].index(temporaryvalue)
#                             temp_list[1][temp_index] += 1
#                 if 1 in temp_list[1]:
#                     d[i] = -100
#                 count += row_number
    
#     # # 遍历所有个体
#     # for i in range(nind):
#     #     count = 0
#     #     # 对每个料条进行循环
#     #     for m in range(AMN):
#     #         # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#     #         if m+1 <= big_row:
#     #             temp_list = []
#     #             # 遍历每一行的每个料堆
#     #             for n in range(int(row_number)+1):
#     #                 # 对每个料堆验证其是否存在固定搭配组合位于该料条
#     #                 temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
#     #                 if temporaryvalue != '0' and temporaryvalue in temp_list:
#     #                     temp_list.remove(temporaryvalue)
#     #                 elif temporaryvalue != '0' and temporaryvalue not in temp_list:
#     #                     temp_list.append(temporaryvalue)
#     #             if len(temp_list) != 0:
#     #                 d[i] = -100
#     #             count += row_number+1
#     #         else:
#     #             for n in range(int(row_number)):
#     #                 # 对每个料堆验证其是否存在固定搭配组合位于该料条
#     #                 temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
#     #                 if temporaryvalue != '0' and temporaryvalue in temp_list:
#     #                     temp_list.remove(temporaryvalue)
#     #                 elif temporaryvalue != '0' and temporaryvalue not in temp_list:
#     #                     temp_list.append(temporaryvalue)
#     #             if len(temp_list) != 0:
#     #                 d[i] = -100
#     #             count += row_number

#     return d

# def efunc (chrom_property, AMN, row_number, nind, big_row, num_property):
#     # 初始化该约束的适应度值
#     e = np.zeros((nind, 1))

#     # 遍历所有个体
#     for i in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for m in range(AMN):
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if m + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number)+1):
#                     # 判断该料堆是否存在偏好位置并注明偏好位置所在料条
#                     temporaryvalue = chrom_property[i, (count + n) * num_property + 6]
#                     if int(temporaryvalue) != 0 and int(temporaryvalue) != m+1:
#                         e[i] = -100
#                 count += row_number+1
#             else:
#                 # 遍历每一行的每个料堆
#                 for n in range(int(row_number)):
#                     # 判断该料堆是否存在偏好位置并注明偏好位置所在料条
#                     temporaryvalue = chrom_property[i, (count + n) * num_property + 6]
#                     if int(temporaryvalue) != 0 and int(temporaryvalue) != m+1:
#                         e[i] = -100
#                 count += row_number
#     return e

# def ffunc(chrom_property, AMN, row_number, nind, big_row, num_property):
#     # 初始化该约束的适应度值
#     f = np.zeros((nind, 1))

#     # 遍历所有个体
#     for i in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for m in range(AMN):
#                 # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#                 if m + 1 <= big_row:
#                     # 遍历每一行的每个料堆
#                     for n in range(int(row_number)):
#                         temporaryvalue = chrom_property[i, (count + n) * num_property + 7]
#                         # 判断该料堆是否存在不相容的组合且位于同一料条
#                         if int(temporaryvalue) != 0 and temporaryvalue == chrom_property[i, (count + n + 1) * num_property + 7]:
#                             f[i] = -100
#                     count += row_number+1
#                 else:
#                     # 遍历每一行的每个料堆
#                     for n in range(int(row_number) - 1):
#                         temporaryvalue = chrom_property[i, (count + n) * num_property + 7]
#                         # 判断该料堆是否存在不相容的组合且位于同一料条
#                         if int(temporaryvalue) != 0 and temporaryvalue == chrom_property[i, (count + n + 1) * num_property + 7]:
#                             f[i] = -100
#                     count += row_number
#     return f

# def gfunc(chrom_property, AMN, row_number,nind, big_row, g_standard, num_property):
#     # 初始化该约束的适应度值
#     g = np.zeros((nind, 1))

#     # 遍历所有个体
#     for m in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for i in range(AMN):
#             Length = int(chrom_property[m, count * num_property + 4])
#             # 计算第一料条的长度作为被减数
#             if i == 0:
#                 # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#                 if i + 1 <= big_row:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number)):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number + 1
#                 else:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number) - 1):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number
#                 pre_row = Length
#             # 计算其余料条长度
#             else:
#                 # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#                 if i + 1 <= big_row:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number)):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number + 1
#                 else:
#                     # 遍历每一行的每个料堆
#                     for j in range(int(row_number) - 1):
#                         Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                     count += row_number
#                 # 令其余料条长度与第一料条长度做差
#                 g[m] += abs(Length-pre_row)
#     g = 2*g_standard/(g+g_standard)*100
#     return g

# def hfunc(chrom_property, AMN, row_number, nind, big_row, h_standard, safe_distance, num_property):
#     # 初始化该约束的适应度值
#     h = np.zeros((nind, 1))

#     # 遍历所有个体
#     for m in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for i in range(AMN):
#             Length = int(chrom_property[m, count * num_property + 4])
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if i + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 # 计算堆取料机假定位置
#                 for j in range(int(row_number)):
#                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                 stacker_position = Length/2

#                 # 遍历每一行的每个料堆
#                 temporary = 0
#                 for j in range(int(row_number)+1):
#                     # 计算料堆的中点位置到堆取料机的距离
#                     item_length = int(chrom_property[m, (count + j) * num_property + 4])
#                     rea_position = temporary + item_length/2
#                     temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
#                     h[m] += abs(stacker_position - rea_position)*(int(chrom_property[m, (count + j) * num_property + 9])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
#                     # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
#                     # temporary_1 = temporary
#                     # temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
#                     # temporary_2 = temporary - safe_distance
#                     # if temporary_2 <= stacker_position:
#                     #     h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
#                     # else:
#                     #     h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
#                 count += row_number + 1
#             else:
#                 # 遍历每一行的每个料堆
#                 # 计算堆取料机假定位置
#                 for j in range(int(row_number) - 1):
#                     Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
#                 stacker_position = Length / 2

#                 # 遍历每一行的每个料堆
#                 temporary = 0
#                 for j in range(int(row_number)):
#                     # 计算料堆的中点位置到堆取料机的距离
#                     item_length = int(chrom_property[m, (count + j) * num_property + 4])
#                     rea_position = temporary + item_length/2
#                     temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
#                     h[m] += abs(stacker_position - rea_position)*(int(chrom_property[m, (count + j) * num_property + 9])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 12])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 15])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 18])
#                                                                      + int(chrom_property[m, (count + j) * num_property + 21]))
                    
#                     # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
#                     # temporary_1 = temporary
#                     # temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
#                     # temporary_2 = temporary - safe_distance
#                     # if temporary_2 <= stacker_position:
#                     #     h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
#                     # else:
#                     #     h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 12])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 15])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 18])
#                     #                                                  + int(chrom_property[m, (count + j) * num_property + 21]))
#                 count += row_number
#     h = 2*h_standard/(h+h_standard)*100
#     return h

# def ifunc ( chrom_property, AMN, row_number,nind,  big_row, i_standard, digraph, diction2, num_property):
#     # 初始化该约束的适应度值
#     i = np.zeros((nind, 1))

#     # 遍历所有个体
#     for m in range(nind):
#         count = 0
#         # 对每个料条进行循环
#         for n in range(AMN):
#             # 判断分堆时剩下的料堆是否以依次插入前几个料堆
#             if n + 1 <= big_row:
#                 # 遍历每一行的每个料堆
#                 for j in range(int(row_number)):
#                     # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
#                     if chrom_property[m, (count + j) * num_property + 8] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
#                             chrom_property[m, (count + j) * num_property + 10])
#                     if chrom_property[m, (count + j) * num_property + 11] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
#                             chrom_property[m, (count + j) * num_property + 13])
#                     if chrom_property[m, (count + j) * num_property + 14] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
#                             chrom_property[m, (count + j) * num_property + 16])
#                     if chrom_property[m, (count + j) * num_property + 17] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
#                             chrom_property[m, (count + j) * num_property + 19])
#                     if chrom_property[m, (count + j) * num_property + 20] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
#                             chrom_property[m, (count + j) * num_property + 22])
#                 count += row_number + 1
#             else:
#                 # 遍历每一行的每个料堆
#                 for j in range(int(row_number) - 1):
#                     # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
#                     if chrom_property[m, (count + j) * num_property + 8] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
#                             chrom_property[m, (count + j) * num_property + 10])
#                     if chrom_property[m, (count + j) * num_property + 11] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
#                             chrom_property[m, (count + j) * num_property + 13])
#                     if chrom_property[m, (count + j) * num_property + 14] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
#                             chrom_property[m, (count + j) * num_property + 16])
#                     if chrom_property[m, (count + j) * num_property + 17] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
#                             chrom_property[m, (count + j) * num_property + 19])
#                     if chrom_property[m, (count + j) * num_property + 20] != '0':
#                         i[m] += digraph[diction2[str(n + 1)]][
#                                     diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
#                             chrom_property[m, (count + j) * num_property + 22])
#                 count += row_number
#     i = 2*i_standard/(i+i_standard)*100
#     return i

