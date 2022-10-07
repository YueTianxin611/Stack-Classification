import geatpy as ga
import numpy as np
import pandas as pd
import time
import pymysql
from sqlalchemy import create_engine
import os


def GA():
    # 读取 ‘料堆属性表.xls’ 文件
    # 设置数据库链接路径
    conn = database_connect("10.75.4.20", "intern", "123456", "financial_dept")
    df = get_data_all(conn, "料堆属性表")
    #df = pd.read_excel('料堆属性表.xls')

    # diction 储存料堆属性，格式ndarray
    diction = np.array(df)

    # row 储存料堆数
    row = len(diction)

    # 为料堆过长的特殊情况，增加一列属性
    temporary = np.zeros((row, 1))
    diction = np.column_stack((diction, temporary))

    # 记录料堆的属性数量
    num_property = len(diction[0])

    # 读取 'weight.xls' excel 文件， 生成 weight 储存权重系数的 ndarray
    tempvalue = get_data_all(conn, "weight")
    weight = tempvalue.ix[0].values

    # 读取 'globalparameter.xls' excel 文件

    tempvalue = get_data_all(conn, "globalparameter")

    # AMN 储存料堆堆数
    AMN = tempvalue.ix[0, 0]

    # safe_distance 储存料堆间安全距离
    safe_distance = tempvalue.ix[0, 1]

    # maxleng 储存每个料条的最大长度
    maxleng = tempvalue.ix[0, 3:]

    # standard_speed_difference 储存相近约束原则的标准
    standard_speed_difference = tempvalue.ix[0, 2]

    # 用料堆数处以行数并向下取整，为一般的料条的料堆个数，并计算余下的料堆数
    if len(diction) % AMN != 0:
        row_number = row // AMN
    else:
        row_number = row // AMN
    big_row = len(diction) % AMN

    """=======================有向图初始化配置========================="""
    # 读取 'digraph.csv' excel 文件， 生成 ndarray 格式的 path 储存所有路径
    edges = get_data_all(conn, "digraph")
    #edges = pd.read_csv(r'digraph.csv', encoding='gb2312')
    path = np.array(edges)

    # 利用 Floyd 算法生成有向图，生成名为 digraph 的数据形式为 list 的有向图， 以及名为 diction2 的字典检索 digraph 用
    temporary, diction2 = [], {}
    for i, ii in enumerate(path):
        if ii[0] not in temporary:
            temporary.append(ii[0])
            diction2[ii[0]] = len(temporary) - 1
        if ii[1] not in temporary:
            temporary.append(ii[1])
            diction2[ii[1]] = len(temporary) - 1
    temporary2 = np.zeros((len(temporary), len(temporary)))
    for i in range(len(temporary2)):
        for j in range(len(temporary2[0])):
            temporary2[i][j] = 100000
    for i in path:
        temporary2[diction2[i[0]]][diction2[i[1]]] = i[2]
        temporary2[diction2[i[1]]][diction2[i[0]]] = i[2]
    digraph = list(np.zeros((len(temporary), len(temporary))))
    for i, ii in enumerate(digraph):
        digraph[i] = list(digraph[i])

    def Floyd(G, D):
        t = 0
        for u in range(0, len(G)):
            for s in range(0, len(G)):
                D[u][s] = G[u][s]
        for k in range(0, len(G)):
            for v in range(0, len(G)):
                for w in range(0, len(G)):
                    if D[v][w] > D[v][k] + D[k][w]:
                        t = t + 1
                        D[v][w] = D[v][k] + D[k][w]

    Floyd(temporary2, digraph)
    for i in range(len(digraph)):
        for j in range(len(digraph[0])):
            if digraph[i][j] == 100000:
                digraph[i][j] = 0

    """==========================初始化配置==========================="""
    # 最大遗传代数
    MAXGEN = 1000
    # 遗传算法的选择方式设为"sus"——轮盘赌算法
    selectStyle = 'sus'
    # 遗传算法的重组方式，设为两点交叉
    recombinStyle = 'xovpm'
    # 染色体数
    nind = 100
    # 染色体长度
    lind = len(diction)
    # 排列集合大小
    varlen = len(diction)
    '''===========================种群生成============================='''
    # 生成初始种群
    Intialchrom = ga.crtpp(nind, lind, varlen)

    # 调用 filter_function 函数，生成足够多满足要求的初始种群，同时返回
    # chrom 新生成的适应度全部不为零的种群， b_standard 相近原则约束的标准
    # g_standard 料条利用水平的约束标准， h_standard 堆取料机运行距离的标准
    # i_standard 堆取料机约束原则的标准, nind 为更新后的种群中个体大小，即料堆数
    # varlen 更新后的排列集合的大小， lind 更新后的染色体长度
    # diction 更新后的料堆属性
    [chrom, b_standard, g_standard, i_standard, h_standard, nind, varlen, lind, diction] = filter_function(Intialchrom,
                                                                                                           diction, AMN,
                                                                                                           row_number,
                                                                                                           nind,
                                                                                                           maxleng,
                                                                                                           standard_speed_difference,
                                                                                                           big_row,
                                                                                                           lind, varlen,
                                                                                                           safe_distance,
                                                                                                           digraph,
                                                                                                           diction2,
                                                                                                           num_property)

    chrom = np.array(chrom)

    # 记录每一代的最优个体
    pop_trace = np.zeros((MAXGEN, 1))

    # 交叉、变异操作前后父代大小于子代大小关系
    SUBPOP = 1

    # 选择操作前后父代大小于子代大小关系
    SUBPOP2 = 0.5

    # 开始计时
    start_time = time.time()

    # 变异可能
    pm = 0.4

    # 交叉可能
    RecOpt = 0.4

    # 开始遗传算法
    for gen in range(MAXGEN):

        # recombine 为交叉操作
        # 'xovpm' 可实现随即位置的互换，并自动为对应的位置增补实现去重
        # chrom 为用于重组的旧种群，每行对应着一个个体
        # RecOpt(可选参数) 指明了种群发生重组的概率，若为缺省值或 None，则默认概率 为 0.7。
        # SUBPOP(可选参数) 为子种群的数量，若为缺省值或 None，则默认为 1
        SelCh = ga.recombin(recombinStyle, chrom, RecOpt, SUBPOP)

        # 突变率设为0.1，进行变异
        SelCh = ga.mutpp(SelCh, varlen, pm)

        # 将子代种群与父代种群合并
        chrom = np.row_stack((SelCh, chrom))
        nind = len(chrom)

        # 计算适应度
        ObjV = aimfuc(chrom, diction, AMN, weight, row_number, nind, maxleng, standard_speed_difference, big_row,
                      safe_distance, b_standard, g_standard, h_standard, i_standard, digraph, diction2, num_property)

        # 记录每一代表现最好的个体，同时移除每一代表现最差的两个个体
        temporary = np.argmax(list(ObjV))
        Best_trace = list(chrom[temporary, :])
        temporary = np.argmin(list(ObjV))
        chrom = np.delete(chrom, temporary, 0)
        ObjV = np.delete(ObjV, temporary, 0)
        temporary = np.argmin(list(ObjV))
        chrom = np.delete(chrom, temporary, 0)
        ObjV = np.delete(ObjV, temporary, 0)

        # selecting 为选择操作
        # SEL_F 为包含低级选择函数名称的字符串如‘rws’
        # chrom 为包含当前种群的染色体矩阵，每一行对应一个个体的一条染色体
        # FitnV 为包含种群 chrom 中个体的适应度值的列向量。
        # SelCh 为返回的子代种群
        # SUBPOP2 为子代种群数大小为父代种群数的一半
        chrom = ga.selecting(selectStyle, chrom, ObjV, SUBPOP2)

        # 记录当代目标函数的最优值
        pop_trace[gen] = np.max(ObjV)

        # ax = ga.sgaplot(pop_trace, '种群最优个体目标函数值', False, ax, gen)

        # 将本代最好的个体插入种群
        chrom = np.row_stack((chrom, Best_trace))

        # 设置提前退出循环的条件
        if gen > 40 and pop_trace[gen] == pop_trace[gen - 40]:
            ga.sgaplot(pop_trace, '种群最优个体目标函数值', True)
            break
        if gen == 1:
            ga.sgaplot(pop_trace, '种群最优个体目标函数值', True)
            break

    # 停止计时
    End_time = time.time()
    print(End_time - start_time)

    # 整理要输出的数据的内容
    temporarylist = []
    for i, ii in enumerate(Best_trace):
        temporary2 = list(diction[int(ii - 1)][:-1])
        temporarylist.append(temporary2)
    count = 0
    for i in range(AMN):
        Length = int(temporarylist[count][4])
        if i + 1 <= big_row:
            for j in range(int(row_number)):
                Length = Length + int(temporarylist[count + j][4])
            stacker_position = Length / 2

            temporary = 0
            for j in range(int(row_number) + 1):
                temporary_1 = temporary
                temporary += int(temporarylist[count + j][4]) + safe_distance
                temporary_2 = temporary - safe_distance
                if temporary_2 <= stacker_position:
                    temporarylist[count + j].append(i + 1)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_2)
                else:
                    temporarylist[count + j].append(i + 1)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_1)
            count += row_number + 1
        else:
            for j in range(int(row_number) - 1):
                Length = Length + int(temporarylist[count + j][4])
            stacker_position = Length / 2

            temporary = 0
            for j in range(int(row_number)):
                temporary_1 = temporary
                temporary += int(temporarylist[count + j][4]) + safe_distance
                temporary_2 = temporary - safe_distance
                if temporary_2 <= stacker_position:
                    temporarylist[count + j].append(i + 1)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_2)
                else:
                    temporarylist[count + j].append(i + 1)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_1)
            count += row_number

    # 将要输出的数据导入，名为'finalresult' 的 excel
    columns = ['堆号', '编号', '名称', '堆存量', '长度', '固定搭配组合', '位置偏好', '产品相容性', '用户1', '取用频次', '消耗速度', '用户2', '取用频次', '消耗速度',
         '用户3', '取用频次', '消耗速度', '用户4', '取用频次', '消耗速度', '用户5', '取用频次', '消耗速度', '料条', '料堆', '起点', '终点', '偏好取料位置']
    # 将要输出的数据导入，名为'finalresult' 的 excel
    temporary = pd.DataFrame(temporarylist,columns=columns)
    insert_data_all(temporary,'intern','123456','10.75.4.20','financial_dept','finalresult')



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

    # 计算料堆容量约束
    a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)

    # 检查是否存在某个料堆过长导致现有的平均堆数放到每个料条的情况难以实现的情况，若出现则将其看作两个料堆
    while max(a) == -100:
        temporary = diction[:, 4]
        temporary = np.argmax(temporary)
        diction[temporary][3] = int(diction[temporary][3] / 2)
        diction[temporary][4] = int(diction[temporary][4] / 2)
        diction[temporary][9] = int(diction[temporary][9] / 2)
        diction[temporary][10] = int(diction[temporary][10] / 2)
        diction[temporary][12] = int(diction[temporary][12] / 2)
        diction[temporary][13] = int(diction[temporary][13] / 2)
        diction[temporary][15] = int(diction[temporary][15] / 2)
        diction[temporary][16] = int(diction[temporary][16] / 2)
        diction[temporary][18] = int(diction[temporary][18] / 2)
        diction[temporary][19] = int(diction[temporary][19] / 2)
        diction[temporary][21] = int(diction[temporary][21] / 2)
        diction[temporary][22] = int(diction[temporary][22] / 2)
        diction[temporary][-1] = len(diction)
        diction = np.row_stack((diction, diction[temporary]))
        lind += 1
        big_row += 1
        varlen = len(diction)
        if big_row == AMN:
            row_number += 1
            big_row = 0
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

    # 计算固定搭配组合的适应度
    d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    # 计算是否又因为固定搭配组合放到一个料堆导致的放不下的情况， 若有将固定搭配组合中最长的料堆拆成两个料堆再进行计算
    while max(d) == -100:
        temporary = diction[:, 4]
        for i,ii in enumerate(temporary):
            if diction[i][5] == 0:
                temporary[i] = 0
        temporary = np.argmax(temporary)
        diction[temporary][3] = int(diction[temporary][3] / 2)
        diction[temporary][4] = int(diction[temporary][4] / 2)
        diction[temporary][9] = int(diction[temporary][9] / 2)
        diction[temporary][10] = int(diction[temporary][10] / 2)
        diction[temporary][12] = int(diction[temporary][12] / 2)
        diction[temporary][13] = int(diction[temporary][13] / 2)
        diction[temporary][15] = int(diction[temporary][15] / 2)
        diction[temporary][16] = int(diction[temporary][16] / 2)
        diction[temporary][18] = int(diction[temporary][18] / 2)
        diction[temporary][19] = int(diction[temporary][19] / 2)
        diction[temporary][21] = int(diction[temporary][21] / 2)
        diction[temporary][22] = int(diction[temporary][22] / 2)
        diction[temporary][-1] = len(diction)
        diction = np.row_stack((diction, diction[temporary]))
        lind += 1
        big_row += 1
        varlen = len(diction)
        if big_row == AMN:
            row_number += 1
            big_row = 0
        Intialchrom = ga.crtpp(nind, lind, varlen)
        chrom_property = np.array(diction[int(Intialchrom[0][0] - 1), :])
        for i in range(len(Intialchrom) - 1):
            chrom_property = np.row_stack((chrom_property, list(diction[int(Intialchrom[i + 1][0] - 1)])))
        for j in range(len(Intialchrom[1]) - 1):
            temporary = list(diction[int(Intialchrom[0][j + 1] - 1)])
            for i in range(len(Intialchrom) - 1):
                temporary = np.row_stack((temporary, list(diction[int(Intialchrom[i + 1][j + 1] - 1)])))
            chrom_property = np.column_stack((chrom_property, temporary))
        d = dfunc(chrom_property, diction, AMN, row_number, nind, big_row, num_property)

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
    while len(chrom) < nind:
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
            if ii == 0:
                chrom.append(list(Intialchrom[i]))
                ObjV.append(list(ii))
    nind = len(chrom)
    chrom_property = list(diction[int(chrom[0][0] - 1), :])
    for i in range(len(chrom) - 1):
        chrom_property = np.row_stack((chrom_property, list(diction[int(chrom[i + 1][0] - 1)])))
    for j in range(len(chrom[1]) - 1):
        temporary = list(diction[int(chrom[0][j + 1] - 1)])
        for i in range(len(chrom) - 1):
            temporary = np.row_stack((temporary, list(diction[int(chrom[i + 1][j + 1] - 1)])))
        chrom_property = np.column_stack((chrom_property, temporary))

    # 计算相近原则的标准
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
                temporary += row_number
    b_standard = max(b)

    # 计算料条利用水平的衡量标准
    # 初始化该约束的适应度值
    g = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])
            # 计算第一料条的长度作为被减数
            if i == 0:
                # 判断分堆时剩下的料堆是否以依次插入前几个料堆
                if i + 1 <= big_row:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number)):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number + 1
                else:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number) - 1):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number
                pre_row = Length
            # 计算其余料条长度
            else:
                # 判断分堆时剩下的料堆是否以依次插入前几个料堆
                if i + 1 <= big_row:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number)):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number + 1
                else:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number) - 1):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number
                # 令其余料条长度与第一料条长度做差
                g[m] += abs(Length-pre_row)
    g_standard = min(g)

    # 计算堆取料机运行距离约束的标准值
    # 初始化该约束的适应度值
    h = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if i + 1 <= big_row:
                # 遍历每一行的每个料堆
                # 计算堆取料机假定位置
                for j in range(int(row_number)):
                    Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                stacker_position = Length/2

                # 遍历每一行的每个料堆
                # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
                temporary = 0
                for j in range(int(row_number)+1):
                    temporary_1 = temporary
                    temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
                    temporary_2 = temporary - safe_distance
                    if temporary_2 <= stacker_position:
                        h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                    else:
                        h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                # 计算堆取料机假定位置
                for j in range(int(row_number) - 1):
                    Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                stacker_position = Length / 2

                # 遍历每一行的每个料堆
                # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
                temporary = 0
                for j in range(int(row_number)):
                    temporary_1 = temporary
                    temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
                    temporary_2 = temporary - safe_distance
                    if temporary_2 <= stacker_position:
                        h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                    else:
                        h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                count += row_number
    h_standard = min(h)

    # 计算原料对公里数运行距离约束
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
                for j in range(int(row_number)):
                    # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
                    if chrom_property[m, (count + j) * num_property + 8] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
                            chrom_property[m, (count + j) * num_property + 10])
                    if chrom_property[m, (count + j) * num_property + 11] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
                            chrom_property[m, (count + j) * num_property + 13])
                    if chrom_property[m, (count + j) * num_property + 14] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
                            chrom_property[m, (count + j) * num_property + 16])
                    if chrom_property[m, (count + j) * num_property + 17] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
                            chrom_property[m, (count + j) * num_property + 19])
                    if chrom_property[m, (count + j) * num_property + 20] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
                            chrom_property[m, (count + j) * num_property + 22])
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                for j in range(int(row_number) - 1):
                    # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
                    if chrom_property[m, (count + j) * num_property + 8] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
                            chrom_property[m, (count + j) * num_property + 10])
                    if chrom_property[m, (count + j) * num_property + 11] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
                            chrom_property[m, (count + j) * num_property + 13])
                    if chrom_property[m, (count + j) * num_property + 14] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
                            chrom_property[m, (count + j) * num_property + 16])
                    if chrom_property[m, (count + j) * num_property + 17] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
                            chrom_property[m, (count + j) * num_property + 19])
                    if chrom_property[m, (count + j) * num_property + 20] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
                            chrom_property[m, (count + j) * num_property + 22])
                count += row_number
    i_standard = min(i)

    # 输出结果
    return chrom, b_standard, g_standard, i_standard, h_standard, nind, varlen, lind, diction


def aimfuc(chrom, diction, AMN, weight, row_number,nind, maxleng, standard_speed_difference,big_row, safe_distance, b_standard , g_standard, h_standard,  i_standard, digraph, diction2, num_property):
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
    if weight[1] != 0:
        a = afunc(chrom_property, AMN, row_number, nind, maxleng, big_row, safe_distance, num_property)

    if weight[2] != 0:
        b = bfunc(chrom_property, AMN, row_number, nind, standard_speed_difference, big_row, b_standard, num_property)

    if weight[3] != 0:
        c = cfunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    if weight[4] != 0:
        d = dfunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    if weight[5] != 0:
        e = efunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    if weight[6] != 0:
        f = ffunc(chrom_property, AMN, row_number, nind, big_row, num_property)

    if weight[7] != 0:
        g = gfunc(chrom_property, AMN, row_number,nind, big_row, g_standard, num_property)

    if weight[8] != 0:
        h = hfunc(chrom_property, AMN, row_number, nind, big_row, h_standard, safe_distance, num_property)

    if weight[9] != 0:
        i = ifunc(chrom_property, AMN, row_number,nind,  big_row, i_standard, digraph, diction2, num_property)

    # 根据权重系数进行计算 + weight[8] * h
    ObjV += weight[1] * a + weight[6] * b + weight[2] * c + weight[3] * d + weight[5] * e + weight[4] * f + weight[7] * g + weight[8] * h + weight[9] * i

    # 将适应度为负的个体的适应度调为0方便之后计算
    for m, i in enumerate(ObjV):
        if i < 0:
            ObjV[m] = 0
    return ObjV

def afunc(chrom_property,AMN,row_number,nind,maxleng,big_row, safe_distance, num_property):
    # 初始化该约束的适应度值
    a = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if i + 1 <= big_row:
                # 遍历每一行的每个料堆
                for j in range(int(row_number)):
                    # 计算每个料条的长度
                    Length = Length + int(chrom_property[m, (count + j + 1) * num_property + 4])
                count += row_number + 1
                # 判断每个料条的长度是否超出该料条所能容纳的最大长度
                if Length >= (maxleng[i] - row_number * safe_distance):
                        a[m] = -100  # 若该料条长度和超过 maxlength 则对应个体的 a 值为 -100
            else:
                # 遍历每一行的每个料堆
                for j in range(int(row_number) - 1):
                    Length = Length + int(chrom_property[m, (count + j + 1) * num_property + 4])
                count += row_number
                # 判断每个料条的长度是否超出该料条所能容纳的最大长度
                if Length >= (maxleng[i] - (row_number-1) * safe_distance):
                        a[m] = -100  # 若该料条长度和超过 maxlength 则对应个体的 a 值为 -100
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
                temporary += row_number
    if b_standard != 0:
        b = b/b_standard*100
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
                        c[i] = -100
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
                        c[i] = -100
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
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if m+1 <= big_row:
                temp_list = []
                # 遍历每一行的每个料堆
                for n in range(int(row_number)+1):
                    # 对每个料堆验证其是否存在固定搭配组合位于该料条
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
                    if temporaryvalue != '0' and temporaryvalue in temp_list:
                        temp_list.remove(temporaryvalue)
                    elif temporaryvalue != '0' and temporaryvalue not in temp_list:
                        temp_list.append(temporaryvalue)
                if len(temp_list) != 0:
                    d[i] = -100
                count += row_number+1
            else:
                for n in range(int(row_number)):
                    # 对每个料堆验证其是否存在固定搭配组合位于该料条
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 5]
                    if temporaryvalue != '0' and temporaryvalue in temp_list:
                        temp_list.remove(temporaryvalue)
                    elif temporaryvalue != '0' and temporaryvalue not in temp_list:
                        temp_list.append(temporaryvalue)
                if len(temp_list) != 0:
                    d[i] = -100
                count += row_number

    # 遍历所有个体
    for i in range(nind):
        count = 0
        # 对每个料条进行循环
        for m in range(AMN):
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if m + 1 <= big_row:
                # 遍历每一行的每个料堆
                for n in range(int(row_number) - 1):
                    # 对每个料堆判断其是否被分堆，若被分堆观察其被分堆的组合是否位于相邻位置
                    if chrom_property[i, (count + n + 1) *num_property + 23] != 0 and chrom_property[i, (count + n + 1) *num_property + 23] != 0.0:
                        if chrom_property[i, (count + n) * num_property + 23] != chrom_property[i, (count + n + 1) *
                                                                                                   num_property + 23] \
                                and chrom_property[i, (count + n + 1) * num_property + 23] != chrom_property[i,
                                                                                                             (
                                                                                                                     count + n + 2) * num_property + 23]:
                            d[i] = -100
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                for n in range(int(row_number) - 2):
                    # 对每个料堆判断其是否被分堆，若被分堆观察其被分堆的组合是否位于相邻位置
                    if chrom_property[i, (count + n + 1) *num_property + 23] != 0 and chrom_property[i, (count + n + 1) *num_property + 23] != 0.0:
                        if chrom_property[i, (count + n) * num_property + 23] != chrom_property[i, (count + n + 1) *
                                                                                                   num_property + 23] \
                                and chrom_property[i, (count + n + 1) * num_property + 23] != chrom_property[i,
                                                                                                             (
                                                                                                                     count + n + 2) * num_property + 23]:
                            d[i] = -100
                count += row_number
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
                        e[i] = -100
                count += row_number+1
            else:
                # 遍历每一行的每个料堆
                for n in range(int(row_number)):
                    # 判断该料堆是否存在偏好位置并注明偏好位置所在料条
                    temporaryvalue = chrom_property[i, (count + n) * num_property + 6]
                    if int(temporaryvalue) != 0 and int(temporaryvalue) != m+1:
                        e[i] = -100
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
                            f[i] = -100
                        break
                    count += row_number+1
                else:
                    # 遍历每一行的每个料堆
                    for n in range(int(row_number) - 1):
                        temporaryvalue = chrom_property[i, (count + n) * num_property + 7]
                        # 判断该料堆是否存在不相容的组合且位于同一料条
                        if int(temporaryvalue) != 0 and temporaryvalue == chrom_property[i, (count + n + 1) * num_property + 7]:
                            f[i] = -100
                        break
                    count += row_number
    return f

def gfunc(chrom_property, AMN, row_number,nind, big_row, g_standard, num_property):
    # 初始化该约束的适应度值
    g = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])
            # 计算第一料条的长度作为被减数
            if i == 0:
                # 判断分堆时剩下的料堆是否以依次插入前几个料堆
                if i + 1 <= big_row:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number)):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number + 1
                else:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number) - 1):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number
                pre_row = Length
            # 计算其余料条长度
            else:
                # 判断分堆时剩下的料堆是否以依次插入前几个料堆
                if i + 1 <= big_row:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number)):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number + 1
                else:
                    # 遍历每一行的每个料堆
                    for j in range(int(row_number) - 1):
                        Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                    count += row_number
                # 令其余料条长度与第一料条长度做差
                g[m] += abs(Length-pre_row)
    g = 2*g_standard/(g+g_standard)*100
    return g

def hfunc(chrom_property, AMN, row_number, nind, big_row, h_standard, safe_distance, num_property):
    # 初始化该约束的适应度值
    h = np.zeros((nind, 1))

    # 遍历所有个体
    for m in range(nind):
        count = 0
        # 对每个料条进行循环
        for i in range(AMN):
            Length = int(chrom_property[m, count * num_property + 4])
            # 判断分堆时剩下的料堆是否以依次插入前几个料堆
            if i + 1 <= big_row:
                # 遍历每一行的每个料堆
                # 计算堆取料机假定位置
                for j in range(int(row_number)):
                    Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                stacker_position = Length/2

                # 遍历每一行的每个料堆
                # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
                temporary = 0
                for j in range(int(row_number)+1):
                    temporary_1 = temporary
                    temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
                    temporary_2 = temporary - safe_distance
                    if temporary_2 <= stacker_position:
                        h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                    else:
                        h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                # 计算堆取料机假定位置
                for j in range(int(row_number) - 1):
                    Length = Length + int(chrom_property[m, (count + j) * num_property + 4])
                stacker_position = Length / 2

                # 遍历每一行的每个料堆
                # 计算每个料堆的两端中离堆取料机较近的一端到堆取料机的距离
                temporary = 0
                for j in range(int(row_number)):
                    temporary_1 = temporary
                    temporary += int(chrom_property[m, (count + j) * num_property + 4]) + safe_distance
                    temporary_2 = temporary - safe_distance
                    if temporary_2 <= stacker_position:
                        h[m] += abs(stacker_position - temporary_2)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                    else:
                        h[m] += abs(stacker_position - temporary_1)*(int(chrom_property[m, (count + j) * num_property + 9])
                                                                     + int(chrom_property[m, (count + j) * num_property + 12])
                                                                     + int(chrom_property[m, (count + j) * num_property + 15])
                                                                     + int(chrom_property[m, (count + j) * num_property + 18])
                                                                     + int(chrom_property[m, (count + j) * num_property + 21]))
                count += row_number
    h = 2*h_standard/(h+h_standard)*100
    return h

def ifunc ( chrom_property, AMN, row_number,nind,  big_row, i_standard, digraph, diction2, num_property):
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
                for j in range(int(row_number)):
                    # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
                    if chrom_property[m, (count + j) * num_property + 8] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
                            chrom_property[m, (count + j) * num_property + 10])
                    if chrom_property[m, (count + j) * num_property + 11] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
                            chrom_property[m, (count + j) * num_property + 13])
                    if chrom_property[m, (count + j) * num_property + 14] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
                            chrom_property[m, (count + j) * num_property + 16])
                    if chrom_property[m, (count + j) * num_property + 17] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
                            chrom_property[m, (count + j) * num_property + 19])
                    if chrom_property[m, (count + j) * num_property + 20] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
                            chrom_property[m, (count + j) * num_property + 22])
                count += row_number + 1
            else:
                # 遍历每一行的每个料堆
                for j in range(int(row_number) - 1):
                    # 对每个料堆计算其五个可能的用户对应的吨公里运输距离的和
                    if chrom_property[m, (count + j) * num_property + 8] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 8]]] * float(
                            chrom_property[m, (count + j) * num_property + 10])
                    if chrom_property[m, (count + j) * num_property + 11] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 11]]] * float(
                            chrom_property[m, (count + j) * num_property + 13])
                    if chrom_property[m, (count + j) * num_property + 14] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 14]]] * float(
                            chrom_property[m, (count + j) * num_property + 16])
                    if chrom_property[m, (count + j) * num_property + 17] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 17]]] * float(
                            chrom_property[m, (count + j) * num_property + 19])
                    if chrom_property[m, (count + j) * num_property + 20] != '0':
                        i[m] += digraph[diction2[str(n + 1)]][
                                    diction2[chrom_property[m, (count + j) * num_property + 20]]] * float(
                            chrom_property[m, (count + j) * num_property + 22])
                count += row_number
    i = 2*i_standard/(i+i_standard)*100
    return i


#mysql
# 存入所有所得数据
def insert_data_all(data,user,password,host,database,table_name):
    """
    作用：
         把全部数据存储进mysql数据库中
    输入：
        data:待存入的数据，为pandas.dataframe结构
        user:用户名
        password:密码
        host:地址
        database:数据库名
        table_name:待存入数据的表名
    """
    connect = create_engine("mysql+pymysql://"+user+":"+password+"@"+host+"/"+database)
    pd.io.sql.to_sql(data, table_name, connect, schema=database, if_exists="append",index=False)
    return


def database_connect(my_host,my_user,my_password,database_name):
    '''
    作用：
        连接数据库
    输入：
        my_host:IP地址
        my_user:用户名
        my_password:密码
        database_name:数据库名
    输出：
        conn:连接路径
    '''
    conn = pymysql.connect(host=my_host,user= my_user,password=my_password,database=database_name)
    return conn


def get_data_all(conn,table_name):
    '''
    作用：
        查询数据库中某表里的全部数据
    输入：
        conn:数据库路径
        table_name:希望查询的表名
    输出：
        data_out:查询得到的数据，格式为pandas.dataframe
    '''
    conn.ping(reconnect=True)
    cursor = conn.cursor()
    data_out = pd.DataFrame()
    sql = "select * from " + table_name
    cursor.execute(sql)
    data = cursor.fetchall()
    data_out = data_out.append(list(data))
    cursor.close()
    conn.close()
    return data_out

def get_data_month(conn,table_name,month):
    conn.ping(reconnect=True)
    cursor = conn.cursor()
    data_out = pd.DataFrame()
    sql = "select * from " + table_name + " where month=" + month
    print(sql)
    cursor.execute(sql)
    data = cursor.fetchall()
    data_out = data_out.append(list(data))
    #data_out = data_out.reset_index(drop=True)
    #data_out.columns = features
    cursor.close()
    conn.close()
    return data_out



GA()