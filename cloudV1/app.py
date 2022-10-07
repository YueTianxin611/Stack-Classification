"""
app.py 云版本-第一版
"""
from flask import Flask, jsonify
app = Flask(__name__)

import geatpy as ga
import pymysql
import numpy as np
import pandas as pd
import time
from sqlalchemy import create_engine
from optim_field import afunc,bfunc,cfunc,dfunc,efunc,ffunc,gfunc,hfunc,ifunc
from optim_field import aimfuc,filter_function,cal_i,cal_g,cal_h,users_velocity
from optim_field import MINVALUE
import os
import json

#mysql

def database_connect(my_host,my_port,my_user,my_password,database_name):
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
    print(my_host,my_port,my_user,my_password,database_name)
    conn = pymysql.connect(host=my_host,port=my_port,user= my_user,password=my_password,database=database_name,charset='utf8')
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

def insert_data_all(data):
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
    #connect = database_connect("10.75.4.20", "intern", "123456", "financial_dept")
    outputsetting = settings["OutputSettings"]
    port = outputsetting["Port"]
    user = outputsetting["User"]
    password = outputsetting["Password"]
    host = outputsetting["Host"]
    database = outputsetting["Database"]
    table_name = outputsetting["TableName"]
    connect = create_engine("mysql+pymysql://"+user+":"+password+"@"+host+":"+port+"/"+database+"?charset=utf8")
    pd.io.sql.to_sql(data,table_name, connect, schema=database, if_exists="append",index=False)
    return

# 获取连接
def get_connect():
    """
    获取链接
    """
    inputsetting = settings["InputSettings"]
    my_host = inputsetting["MyHost"]
    my_port = int(inputsetting["Port"])
    user = inputsetting["User"]
    password = inputsetting["Password"]
    database_name = inputsetting["DatabaseName"]
    conn = database_connect(my_host,my_port,user,password,database_name)
    return conn

# 输出结果到数据库
def output2db(temporarylist):  
    """
    该部分代码根据实际运行进行修改
    """
    columns = ['pile_number', 'material_id', 'pile_name', 'pile_stock', 'pile_length', 'fixed_pair', 'position_perferred', 'compatibility', 'user_1', 'frequency_1', 'speed_1', 'user_2', 'frequency_2', 'speed_2',
         'user_3', 'frequency_3', 'speed_3', 'user_4', 'frequency_4', 'speed_4', 'user_5', 'frequency_5', 'speed_5', 'pile_yard', 'pile_order', 'pile_start', 'pile_end', 'pile_picked']
    # 将要输出的数据导入，名为'finalresult' 的 excel
    temporary = pd.DataFrame(temporarylist,columns=columns)
    IDlist = [int(time.time())]*len(temporary)
    temporary.insert(0,'cal_id',IDlist)
    insert_data_all(temporary)
    return True

def get_optim(conn,nind):

    """
    conn：mysql数据库连接对象
    MAXGEN： 最大遗传代数
    nind：    # 染色体数
    """
    # 设置数据库链接路径
    
    df = get_data_all(conn, "pile_attri")
    print(118)
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
    print(141)
    # AMN 储存料堆堆数
    AMN = tempvalue.ix[0, 0]
    AMN = int(AMN)
    # print(AMN,AMN)

    # safe_distance 储存料堆间安全距离
    safe_distance = tempvalue.ix[0, 1]

    # maxleng 储存每个料条的最大长度
    maxleng = tempvalue.ix[0, 3:2+AMN]
    maxleng = list(maxleng)
    # print(maxleng)

    # standard_speed_difference 储存相近约束原则的标准
    standard_speed_difference = tempvalue.ix[0, 2]
    
    # print("OK")
    # 用料堆数处以行数并向下取整，为一般的料条的料堆个数，并计算余下的料堆数
    if len(diction) % AMN != 0:
        row_number = row // AMN
    else:
        row_number = row // AMN
    big_row = len(diction) % AMN

    # 读取料条名称
    sname = get_data_all(conn,'slice_name')
    slice_name = []  # [实际的料条名称]  # 默认对应编号0，1，2...(AMN-1)
    sname = list(sname.iloc[:AMN,-1])
    for ind in range(len(sname)):
        slice_name.append(sname[ind])        

    """=======================有向图初始化配置========================="""
    MAXVALUE = 20000000  
    # 读取 'digraph.csv' excel 文件， 生成 ndarray 格式的 path 储存所有路径
    edges = get_data_all(conn, "digraph")
    print(168)
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
    # # 判断如果没有到所有料条的路径，则报错
    # for item in range(1,AMN+1):
    #     if str(item) not in diction2:
    #         print("Error:",item,"is not approachable.")
    #         raise Exception ("Graph is incomplete！")

    temporary2 = np.zeros((len(temporary), len(temporary)))
    for i in range(len(temporary2)):
        for j in range(len(temporary2[0])):
            temporary2[i][j] = MAXVALUE   # 足够大的初始值
        temporary2[i][i] = 0
    for i in path:
        temporary2[diction2[i[0]]][diction2[i[1]]] = i[2]   # 单向图
        # temporary2[diction2[i[1]]][diction2[i[0]]] = i[2]
    # # 创建空表
    # digraph = list(np.zeros((len(temporary), len(temporary))))
    # for i, ii in enumerate(digraph):
    #     digraph[i] = list(digraph[i])

    # 得到所有的料条和用户的名称
    slicelist = slice_name
    userlist,max_ves,min_ves = users_velocity(df,AMN)

    def dijkstra(graph, s, nslice):    # 判断图是否为空，如果为空直接退出；默认为按物料流动方向建立单向图
        if graph is None:
            raise Exception ("Digraph is empty!")   
        S = [] # 存储进入过队列的元素 
        Q = [i for i in range(len(graph)) if 0<i<MAXVALUE]    
        dist = [i for i in graph[s]]
        dist[s] = 0
        count = 0  # 记录是否有环
        while Q:
            step = Q[0]
            Q.pop(0)
            S.append(step)
            for next_step in range(len(graph[step])):
                if 0<graph[step][next_step]<MAXVALUE:
                    if dist[next_step] > dist[step]+graph[step][next_step]:
                        dist[next_step] = dist[step]+graph[step][next_step]
                    if next_step not in S:
                        Q.append(next_step)
                    else:
                        count += 1 
        if count > 0:
            print("For Slice",nslice,", path graph has",count,"circle(s).")
    
        return dist
        
    # 计算料条到各节点的路径
    digraph = []
    print(220)
    for item in slicelist:
        first = diction2[item]
        distance = dijkstra(temporary2, first, item)
        digraph.append(distance)
        print(226)
    # 判断如果料条到用户的路径不全，则报错
    for item in digraph:
        if min(item)==MAXVALUE:
            print("Error:",item,"is not approachable.")
            raise Exception ("Graph is incomplete！")

    # def Floyd(G, D):
    #     D = G
    #     t = 0
    #     m = 1
    #     while m == 1: 
    #         m = 0
    #         # t += 1
    #         # print(t)
    #         # G = D.copy()
    #         for k in range(0, len(D)):
    #             for v in range(0, len(D)):
    #                 # if D[k][v] < 0:
    #                 for w in range(0, len(G)):
    #                     temp = min(D[k][w], D[w][v])
    #                     if temp >= 0 and (D[k][v] < 0 or D[k][v]>(D[k][w] + D[w][v])):
    #                         # t = t + 1
    #                         m = 1
    #                         D[k][v] = D[k][w] + D[w][v]
                        
    #     return D
    # dig_max = digraph
    
    # digraph,dig_max = Floyd(temporary2, digraph, dig_max)
    # digraph = Floyd(temporary2, digraph)
    
    # for i in range(len(digraph)):
    #     for j in range(len(digraph[0])):
    #         if digraph[i][j] == -1:
    #             digraph[i][j] = 0
    
    dig_min = digraph

    # for item in dig_min:
    #     print(item)

    # 计算i_standard
    imax, imin = cal_i(dig_min,diction2,AMN,slicelist,userlist,max_ves,min_ves)

    # 计算h_standard
    h_standard = cal_h(maxleng,df)

    """==========================初始化配置==========================="""

    
    # 遗传算法的选择方式设为"sus"——轮盘赌算法
    selectStyle = 'sus'
    # 遗传算法的重组方式，设为两点交叉
    recombinStyle = 'xovpm'
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
    set_nind = nind
    [chrom, nind, varlen, lind, diction, big_row] = filter_function(Intialchrom,\
        diction, AMN, row_number, nind, maxleng, standard_speed_difference,\
            big_row, lind, varlen, safe_distance, digraph, diction2,\
                num_property)
    
    # 如果超过循环次数且无解，则认为是无解的
    if nind == 0 :
        print("Suggest reset parameters!")
        raise Exception ("No answer!")
    
    # 初始化部分约束适应度
    b_standard = 2
    g_standard = cal_g(df,maxleng,safe_distance,AMN)
    i_min,i_max = imin, imax
    
    print("307")
    # 如果经过循环后有解，但是解小于染色体个数
    if nind < set_nind:
        print(310)
        ObjV,judge = aimfuc(chrom, diction, AMN, weight, row_number, nind, maxleng, standard_speed_difference, big_row,
                    safe_distance, b_standard, g_standard, h_standard, i_max, i_min, digraph, diction2, num_property)           
        
        # 当所有约束的最大值均为MINVALUE时认为几乎不存在合理解
        rejudge = [[] for _ in range(len(judge[0]))]
        for item in judge:
            njudge = np.max(item)
            for it in range(len(judge[0])):
                rejudge[it].append(float(item[it]))
        # if np.max(njudge) == MINVALUE:
        #     print("Probably no solution! If there were an solution, only come from mutation.\
        #         Suggest reset parameters.")
        #     raise Exception ("Probably no solution!")
        
        temp_ObjV = ObjV
        while max(temp_ObjV)>0:
            bestscore = np.max(list(temp_ObjV))  # 记录最佳适应度
            bestind = np.argmax(list(temp_ObjV))   # 记录最佳适应度对应的序号
            temporary = bestind
            if np.min(rejudge[temporary]) >= 0:   # 最好个体须满足适应度中没有MINVALUE的值
                break
            temp_ObjV[temporary] = MINVALUE
            # temp_ObjV = np.delete(temp_ObjV, temporary, 0)
        
        if max(temp_ObjV)>0:
            Best_trace = chrom[temporary]
            ind = bestind
            par = []
            for item in judge:
                item = list(item)
                # # 打印最后一代的各参数
                # for name in item:
                #     print(name,end=' ')
                # print('')
                par.append(float(item[ind]))
            print("Have at least a solution, but no enough.")
            print("Scores of best answer:")
            print(par)  # judge内依次存入初始权重不为零的a,b,c...参数的所有值
            print("Order of best answer:")
            print(Best_trace)
        else:
            print("Probably no solution! If there were an solution, only come from mutation.\
                Suggest reset parameters.")
            raise Exception ("Probably no solution!")
            
    else:  # 正常情况进行遗传算法
        print(357)
        chrom = np.array(chrom)

        # 记录每一代的最优个体
        pop_trace = np.zeros((MAXGEN, 1))

        # # 交叉、变异操作前后父代大小于子代大小关系
        # SUBPOP = 1

        # # 选择操作前后父代大小于子代大小关系
        # SUBPOP2 = 0.5

        # # 变异可能
        # PM = 0.4

        # # 交叉可能
        # RecOpt = 0.4

        # # 循环退出调节，最大稳定迭代数
        # MAX_STABLE = 40

        # 开始计时
        start_time = time.time()

        # 开始遗传算法
        have_solution = 0   # 默认无解
        for gen in range(MAXGEN):
            print(gen,384)

            # recombine 为交叉操作
            # 'xovpm' 可实现随即位置的互换，并自动为对应的位置增补实现去重
            # chrom 为用于重组的旧种群，每行对应着一个个体
            # RecOpt(可选参数) 指明了种群发生重组的概率，若为缺省值或 None，则默认概率 为 0.7。
            # SUBPOP(可选参数) 为子种群的数量，若为缺省值或 None，则默认为 1
            SelCh = ga.recombin(recombinStyle, chrom, RecOpt, SUBPOP)

            # 突变率设为0.1，进行变异
            SelCh = ga.mutpp(SelCh, varlen, PM)

            # 将子代种群与父代种群合并
            chrom = np.row_stack((SelCh, chrom))
            nind = len(chrom)

            # 计算适应度
            ObjV,judge = aimfuc(chrom, diction, AMN, weight, row_number, nind, maxleng, standard_speed_difference, big_row,
                        safe_distance, b_standard, g_standard, h_standard, i_max, i_min, digraph, diction2, num_property)

            # 整理适应度，当所有约束的最大值均为MINVALUE时认为几乎不存在合理解
            rejudge = [[] for _ in range(len(judge[0]))]
            for item in judge:
                njudge = np.max(item)
                for it in range(len(judge[0])):
                    rejudge[it].append(float(item[it]))
            # if np.max(njudge) == MINVALUE:
            #     print("Probably no solution! If there were an solution, only come from mutation.\
            #         Suggest reset parameters.")
            #     raise Exception ("Probably no solution!")

            # 记录每一代表现最好的个体，同时移除每一代表现最差的一个个体
            # temporary = np.argmax(list(ObjV))
            # 得到最好个体的序号
            temp_ObjV = ObjV
            while max(temp_ObjV)>0:
                bestscore = np.max(list(temp_ObjV))  # 记录最佳适应度
                bestind = np.argmax(list(temp_ObjV))   # 记录最佳适应度对应的序号
                temporary = bestind
                if np.min(rejudge[temporary]) >= 0:   # 最好个体须满足适应度中没有MINVALUE的值
                    break
                temp_ObjV[temporary] = MINVALUE
                # temp_ObjV = np.delete(temp_ObjV, temporary, 0)
            Best_trace = list(chrom[temporary, :])

            if MINVALUE in temp_ObjV:    # 如果适应度存在无解情况，则删除该值
                temporary = temp_ObjV.index(MINVALUE)
            else:       # 否则删除适应度最小的值
                temporary = np.argmin(list(ObjV))
            # temporary = np.argmin(list(ObjV))
            chrom = np.delete(chrom, temporary, 0)
            ObjV = np.delete(ObjV, temporary, 0)
            if MINVALUE in temp_ObjV:
                temporary = temp_ObjV.index(MINVALUE)
            else:       
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

            # # 记录当代目标函数的最优值
            # pop_trace[gen] = np.max(ObjV)
            # 记录上一代目标函数的最优值
            pop_trace[gen] = bestscore

            # ax = ga.sgaplot(pop_trace, '种群最优个体目标函数值', False, ax, gen)

            # 将本代最好的个体插入种群
            chrom = np.row_stack((chrom, Best_trace))

            # 设置提前退出循环的条件
            if gen > MAX_STABLE and abs(pop_trace[gen]-pop_trace[gen - MAX_STABLE])<10**(-6):
                have_solution = 1 # 有解
                # 判断该解是否为合理解
                ind = bestind
                par = []
                for item in judge:
                    item = list(item)
                    # # 打印最后一代的各参数
                    # for name in item:
                    #     print(name,end=' ')
                    # print('')
                    par.append(float(item[ind]))
                print("Have enough solutions.")
                print("Scores of best answer:")
                print(par)  # judge内依次存入初始权重不为零的a,b,c...参数的所有值
                print("Order of best answer:")
                print(Best_trace)
                # if len(temp_ObjV)==0:
                #     print("Suggest reset parameters.")
                #     print("No answer")
                # if float(MINVALUE) in par:
                #     print("exist MINVALUE")

                    # raise Exception ("Don't have a solution!")
                #ga.sgaplot(pop_trace, '种群最优个体目标函数值', True)
                break
            
            # 测试代码    
            # if gen == 2:
            #     ga.sgaplot(pop_trace, '种群最优个体目标函数值', True)
            #     break

    # 设置遗传代数过少，退出循环
    if have_solution == 0:
        print("Generation number is too small. Please set a larger MAXGEN.")
        raise Exception("No solution!")

    # 停止计时
    End_time = time.time()
    print('cal_id:',int(time.time()),'  Time used:', End_time - start_time)

    # 整理要输出的数据的内容
    temporarylist = []
    for i, ii in enumerate(Best_trace):
        temporary2 = list(diction[int(ii - 1)][:-1])
        temporarylist.append(temporary2)
    count = 0
    for i in range(AMN):
        uname = slice_name[i] # 转换料条名称
        Length = int(temporarylist[count][4])+2*safe_distance
        if i + 1 <= big_row:
            for j in range(int(row_number)):
                Length = Length + int(temporarylist[count + j + 1][4])+safe_distance
            stacker_position = Length / 2

            temporary = safe_distance
            for j in range(int(row_number) + 1):
                temporary_1 = temporary
                temporary += int(temporarylist[count + j][4]) + safe_distance
                temporary_2 = temporary - safe_distance
                if temporary_2 <= stacker_position:
                    temporarylist[count + j].append(uname)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_2)
                else:
                    temporarylist[count + j].append(uname)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_1)
            count += row_number + 1
        else:
            for j in range(int(row_number) - 1):
                Length = Length + int(temporarylist[count + j+1][4])+safe_distance
            stacker_position = Length / 2

            temporary = safe_distance
            for j in range(int(row_number)):
                temporary_1 = temporary
                temporary += int(temporarylist[count + j][4]) + safe_distance
                temporary_2 = temporary - safe_distance
                if temporary_2 <= stacker_position:
                    temporarylist[count + j].append(uname)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_2)
                else:
                    temporarylist[count + j].append(uname)
                    temporarylist[count + j].append(j + 1)
                    temporarylist[count + j].append(temporary_1)
                    temporarylist[count + j].append(temporary_2)
                    temporarylist[count + j].append(temporary_1)
            count += row_number

    return temporarylist
    

@app.route("/todo/api/v1.0/GA", methods=['GET'])
def get_tasks():
    print("234")
    conn = get_connect()
    print("2345")
    #print('conn is successful')
    res_array = get_optim(conn,nind)
    output2db(res_array)
    return 'great'

def read_settings(filename):
    # 读取本地json文件
    with open(filename,'r') as load_f:
        data = json.load(load_f)
        # print(list(data))
    # print(list(data))
    load_f.close()
    settings = data
    # print(settings)
    return settings

if __name__ == '__main__':
    # 初始设置
    settings = read_settings("appsettings.json")
   
    # 设置初始参数
    global nind
    para = settings["Parameters"]
    SUBPOP = para["SUBPOP"]      # 设置交叉变异操作的子种群大小
    SUBPOP2 = para["SUBPOP2"]   # 设置选择操作的子种群大小
    PM = para["PM"]        # 变异概率
    RecOpt = para["RecOpt"]    # 交叉概率
    MAX_STABLE = para["MAX_STABLE"] # 循环退出调节，最大稳定迭代数 
    MAXGEN = para["MAXGEN"]  # 遗传代数上限
    nind = para["nind"]   # 染色体数量
    
    apisetting = settings["APISettings"]
    # 运行app
    if apisetting["use"] == 1:
        api_host = apisetting["host"]
        api_port = apisetting["port"]
        api_debug = apisetting["debug"]
        if api_debug==1:
            api_debug = True
        if api_debug==0:
            api_debug = False
        app.run(host=api_host,port=api_port, debug=api_debug)
    else:
        print("2348")
        conn = get_connect()
        print("23459")
        #print('conn is successful')
        res_array = get_optim(conn,nind)
        output2db(res_array)
