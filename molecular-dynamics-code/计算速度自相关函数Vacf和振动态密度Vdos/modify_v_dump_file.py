import os
import re

def file_handling(imports,prints=''):
    # 如果没有指定的话，输出路径默认桌面
    prints = os.path.join(os.path.expanduser("~"), 'Desktop') if prints == '' else prints
    # 判断路径是否正确
    if not os.path.isfile(imports):
        return '输入的路径不是一个文件，请检查路径是否正确。'

    with open(imports,'r',encoding='utf-8') as f:
        new_f = re.sub('(?s)ITEM: TIMESTEP.*?ITEM: ATOMS vx vy vz\W+','',f.read())
    filepath, fullflname = os.path.split(imports)
    fname, ext = os.path.splitext(fullflname)
    with open(os.path.join(prints,f'{fname}2{ext}'),'w',encoding='utf-8') as file:
        file.write(new_f)


if __name__ == '__main__':
    imports = input('请输入文件路径：')
    prints = input('请输入输出路径：[直接回车默认桌面]')

    file_handling(imports,prints)
