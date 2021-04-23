# 该代码用来删除xyz文件每帧的前两行，对各种xyz文件都适用


import re,datetime

a = datetime.datetime.now()
with open("final.xyz","r") as f:
    with open('final2.xyz','w') as f2:
        f2.write(re.sub('\d+\nAtoms.*?\n','',f.read()))
print('完成',datetime.datetime.now() - a)

    