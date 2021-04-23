import re
import os

def file_handling(imports,offset,prints=''):
    # 如果没有指定的话，输出路径默认桌面
    prints = os.path.join(os.path.expanduser("~"), 'Desktop') if prints == '' else prints
    # 判断路径是否正确
    if not os.path.isfile(imports):
        return '输入的路径不是一个文件，请检查路径是否正确。'


    with open(f'{imports}','r') as f:
        list1 = f.readlines()
        # 文件很大，获取内容后就关闭 释放内存
        f.close()

    # 保存处理后的内容
    content_list = ""
    # 统计处理了几行，可以删除
    num = 0
    for lie in list1:
        # 按空格分割成小列表
        line = lie.split()
        # 判断开头，不是4000或者Atoms开头. Timestep: 2912000
        if re.match('^(?!\d|[Atoms])',line[0]):
           
            lie = f"{line[0]}"
            for lin in line[1:]:
                # 判断是否小于100
                if float(lin) < 100:
                    # 将结果加 3 后添加进lin里暂存
                    lie += ' {:.5f}'.format(float(lin) + float(offset))
                else:
                    # 大于 100 不用处理，直接暂存
                    lie += f' {lin}'

            # 处理完添加进上面保存
            content_list += lie+'\n'
        else:
            content_list += lie

        # 打印出来已处理几行
        print(num)
        num += 1

    # 文件保存
    with open(f'{prints}\坐标调整后的xyz文件.txt','w') as ff:
        ff.write(content_list)

    print('处理完成')


if __name__ == '__main__':
    imports = input('请输入文件路径：')
    offset = input('请输入偏移量：')
    prints = input('请输入输出路径：[直接回车默认桌面]')

    file_handling(imports, offset, prints)
