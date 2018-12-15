# -*- coding: utf-8 -*-
"""
Created on Sun May 20 23:15:26 2018

"""
import re
# 父类 存放 染色体号 起始和结束位置
class Genome_info:
    def __init__(self):
        self.chr = ""
        self.start = 0
        self.end = 0
        
# Gene的 id 和  方向(每个gene特有的属性)
class Gene(Genome_info):
    def __init__(self):
        
        # 加上这句是为了优先使用父类的构造函数
        # 有两种方法可以使用
        #Genome_info.__init__(self)
        super().__init__()
        self.orientation = ""
        self.id = ""
        
# 转录本的id和parent(所属gene的id)
class Transcript(Genome_info):
    def __init__(self):
        Genome_info.__init__(self)
        self.id = ""
        self.parent = ""
        
# 外显子 只有parent(所属的转录本的id)
class Exon(Genome_info):
    def __init__(self):
        Genome_info.__init__(self)
        self.parent = ""
        
# 因为gene和transcript都是有id的，所以用使用dic  
# exon的信息统计到transcript
dic_gene = {}
dic_transcript = {}
lis_exon = []
lis_chr = []
with open(r'hg38.txt','r') as f:
    for line in f:      
        if line.startswith('#'):
            continue
        lines = line.strip().split('\t')
            
        # 染色体号 类型 开始结束坐标  方向  属性
        chr_id = lines[0]
        type_ = lines[2]
        start = int(lines[3])
        end = int(lines[4])
        orientation = lines[6]
        # attr这用分号分隔了多个信息
        attr = lines[8]
            
        # 若基因不存在蛋白编码，直接跳过
        # 有蛋白编码的基因特点是都带cds
        if not re.search(r'protein_coding',attr):
            continue
                
        # 更新chr的list 
        if not chr_id in lis_chr:
            lis_chr.append(chr_id)
                
        # 若类型为基因 则初始化一个Gene对象 最后将gene信息存到字典
        if type_ == 'gene':
            # 实例化一个 Gene 对象
            gene = Gene()
            
            # attr字段  gene_id "ENSG00000223972";  只提取双引号里的值
            # gourp(1)表示返回的第一个分组匹配成功的子串　即 ENSG00000223972
            # group(0)表示返回匹配成功的整个字符串      即 gene_id "ENSG00000223972";
            gene_id = re.search(r'gene_id "(.*?)";',attr).group(1)  
            gene.chr_id = chr_id
            gene.start = start
            gene.end = end
            gene.id = gene_id
            gene.orientation = orientation
            dic_gene[gene_id] = gene   
            #print(gene_id)
            
        # 转录本
        elif type_ == 'transcript':
            transcript = Transcript()
            transcript_id = re.search(r'transcript_id "(.*?)";',attr).group(1)
            # parent就是所属的gene id
            parent = re.search(r'gene_id "(.*?)";',attr).group(1)
            # 排除没有所属基因的转录本或者文件存在缺失(向上找是找不到所属信息)
            if not parent in dic_gene:
                continue
            transcript.chr_id = chr_id
            transcript.start = start
            transcript.end = end
            transcript.id = transcript_id
            transcript.parent = parent
            # 以transcript的id作为key 其信息作为value
            dic_transcript[transcript_id] = transcript
            #print(parent)
              
            # 外显子
        elif type_ == 'exon':
            exon = Exon()
            # exon自身没有id 所以只有parent parent是所属transcript的id
            parent = re.search(r'transcript_id "([^;]+)";?', attr).group(1)
            # 判断parent在不在transcript的dic中 道理同上
            if not parent in dic_transcript:
                continue
            exon.chr_id = chr_id
            exon.start = start
            exon.end = end
            exon.parent = parent
            # 统计完后将信息添加到转录本的list中
            lis_exon.append(exon)
            #print(lis_exon)
# 以上基本信息统计完成，开始计算所需的数据
# 基因在染色体上的分布,即每条染色体上有多少个基因
def chr_gene():
    print('染色体上的分布数目')
    count_gene = {}
    '''
    因为gene的信息都存在了 dic_gene 中 
    所以只需遍历其dic中的value即可
    每条染色体上的相同id的基因
    dic_gene中存了  chr_id  start  end  gene_id  gene_orientation
    只需统计 chr_id 相同的基因即可
    '''
    for info in dic_gene.values():
        chr_id = info.chr_id
        if chr_id in count_gene:
            count_gene[chr_id] += 1
        else:
            count_gene[chr_id] = 1
    
    with open("chr_gene.txt",'w') as f:
        for chr_id,j in count_gene.items():
            s = "{0}---{1}".format(str(chr_gene),str(j))
            f.write(''.join(s)+'\n')

# 基因长度的分布情况
def gene_len():
    print('基因长度的分布情况')
    with open("gene_len.txt",'w') as f:
        for gene_id,info in dic_gene.items():
            len = info.end - info.start + 1
            s = "{0}---{1}".format(gene_id,str(len))    
            f.write(''.join(s)+'\n')
#if __name__ == '__main__':
#    gene_len()
# 每个基因转录本数量的分布
        
def gene_transcript():
    print('转录本的分布情况')
    count_transcript = {}
    for info in dic_transcript.values():
        gene_id = info.parent
        if gene_id in count_transcript:
            count_transcript[gene_id] += 1
        else:
            count_transcript[gene_id] = 1
    with open("gene_transcript.txt",'w') as f:
        for geng_id,count_transcript in count_transcript.items():
            s = "{0}---{1}".format(gene_id,str(count_transcript))
            f.write(''.join(s) + '\n')
#if __name__ == '__main__':
#    gene_transcript()
            
# 转录本下外显子数量的统计
def transcript_exon():
    print('转录本的外显子数量统计')
    count_exon = {}
    for exon in lis_exon:
        transcript_id = exon.parent
        if transcript_id in count_exon:
            count_exon[transcript_id] += 1
        else:
            count_exon[transcript_id] = 1
    with open("transcript_exon.txt",'w') as f:
        for id,num in count_exon.items():
            s = "{0}---{1}".format(id,str(num))
            f.write(''.join(s) + '\n')
#if __name__ == '__main__':
#    transcript_exon()
        
# 统计gene下exon的数量分布
'''
exon.parent = transcript  transcript.parent = gene →→  exon--gene
'''
def exon_gene():
    print('gene下exon的数量统计')
    exon_gene = {}
    for exon in lis_exon:
        for transcript_id,transcript in dic_transcript.items():
            for gene in dic_gene.keys():
                if (exon.parent == transcript and transcript.parent == gene) in exon_gene:
                    exon_gene[gene] += 1
                else:
                    exon_gene[gene] = 1
    with open("exon_gene.txt",'w') as f:
        for i,j in exon_gene.items():
            s = "{0}---{1}".format(i,str(j))
            f.write(''.join(s) + '\n')
        
# 外显子坐标统计
def exon_pos():
    print('外显子的坐标统计')
    count_exon = {}
    for exon in lis_exon:
        transcript_id = exon.parent
        if transcript_id in count_exon:
            count_exon[transcript_id] += ',%s - %s'%(str(exon.start),str(exon.end))
        else:
            count_exon[transcript_id] = '%s - %s'%(str(exon.start),str(exon.end))
    with open("exon_pos.txt",'w') as f:
        for i,j in count_exon.items():
            s = "{0}---{1}".format(i,str(j))
            f.write(''.join(s) + '\n')
#if __name__ == '__main__':
#    exon_pos()
        
if __name__ == '__main__':
    chr_gene()
    gene_transcript()
    transcript_exon()
    exon_gene()
    exon_pos()
    
    
