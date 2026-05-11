# -*- coding: utf-8 -*-
import os
import gseapy as gp
from gseapy.plot import barplot, dotplot
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np



def custom_dotplot(data, title, outfile):
    top_terms = data.nlargest(10, 'Count')
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(top_terms['Count'], 
                          top_terms['Term'], 
                          s=top_terms['Count']*10,           
                          c=-np.log10(top_terms['Adjusted P-value']), 
                          cmap='coolwarm', alpha=0.7)
    plt.colorbar(scatter, label='-log10(Adj P)')
    plt.xlabel('Count')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

def custom_barplot(data, title, outfile):
    # data 是包含 'Term', 'Adjusted P-value', 'Count' 的 DataFrame
    top_terms = data.nlargest(10, 'Count')  # 取 Count 最大的 10 个条目
    plt.figure(figsize=(8, 6))
    sns.barplot(data=top_terms, x='Count', y='Term', palette='viridis')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

# 用于识别gene列的关键词
GENE_COLUMN_KEYWORDS = ['gene']

def detect_gene_column(df):
    """自动检测DataFrame中的基因名列"""
    for col in df.columns:
        col_lower = str(col).lower().strip()
        for keyword in GENE_COLUMN_KEYWORDS:
            if keyword in col_lower:
                return col
    return None

def parse_csv_files(file_paths):
    """解析多个CSV文件，提取基因列表"""
    all_genes = set()
    parsed_info = []
    
    for file_path in file_paths:
        try:
            df = pd.read_csv(file_path)
            gene_col = detect_gene_column(df)
            
            if gene_col is None:
                return {
                    'success': False,
                    'error': f"无法在文件 {os.path.basename(file_path)} 中识别基因列。请确保基因列名包含以下关键词之一: {', '.join(GENE_COLUMN_KEYWORDS)}"
                }
            
            genes = df[gene_col].dropna().astype(str).tolist()
            genes = [g.strip() for g in genes if g.strip()]
            all_genes.update(genes)
            
            parsed_info.append({
                'filename': os.path.basename(file_path),
                'gene_column': gene_col,
                'gene_count': len(genes)
            })
            
        except Exception as e:
            return {
                'success': False,
                'error': f"解析文件 {os.path.basename(file_path)} 时出错: {str(e)}"
            }
    
    return {
        'success': True,
        'genes': list(all_genes),
        'parsed_info': parsed_info
    }

def run_kegg_go_analysis(gene_list, gene_sets=None, organism='human', outdir=None):
    """
    执行KEGG和GO富集分析
    
    参数:
    gene_list: 基因列表，可以是文件路径或基因名称列表
    gene_sets: 基因集库列表，默认为['KEGG_2021', 'GO_Biological_Process_2021', 'GO_Cellular_Component_2021', 'GO_Molecular_Function_2021']
    organism: 生物类型，默认为'human'
    outdir: 输出目录
    
    返回:
    dict: 包含分析结果和图表路径的字典
    """

    if isinstance(gene_list, str):
        if os.path.exists(gene_list):
            pass
        else:
            gene_list = [g.strip() for g in gene_list.strip().split('\n') if g.strip()]
    elif isinstance(gene_list,(list,tuple,set)):
        gene_list = list(gene_list)
    else:
        raise ValueError(f"gene_list must be str or list, tuple, set, got {type(gene_list)}")

    # 设置默认基因集库（使用最新版本）
    if gene_sets is None:
        gene_sets = [
            'KEGG_2021_Human',
            'GO_Biological_Process_2021',
            'GO_Cellular_Component_2021',
            'GO_Molecular_Function_2021'
        ]
    
    # 确保输出目录存在
    output_dir='kegg_go_results'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # 执行富集分析
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets,
            organism=organism,
            outdir=outdir,
            cutoff=0.05
        )
        
        enr.results['Count'] = enr.results['Overlap'].str.split('/').str[0].astype(int)
        # 生成图表
        results = {}
        plot_paths = []
        
        for gene_set in gene_sets:
            # 过滤当前基因集的结果
            subset = enr.results.loc[enr.results["Gene_set"] == gene_set]
            
            if not subset.empty:
                results[gene_set] = subset[['Gene_set', 'Term', 'Count', 'P-value', 'Adjusted P-value', 'Odds Ratio', 'Genes']].to_dict('records')

                if 'GO' in gene_set:
                    # 生成柱状图
                    barplot_path = os.path.join(output_dir, f'{gene_set}_barplot.png')
                    custom_barplot(subset, title=gene_set, outfile=barplot_path)
                    plot_paths.append({'name': f'{gene_set} 柱状图', 'path': barplot_path})
                elif 'KEGG' in gene_set:
                    # 生成气泡图
                    dotplot_path = os.path.join(output_dir, f'{gene_set}_dotplot.png')
                    custom_dotplot(subset, title=gene_set, outfile=dotplot_path)
                    plot_paths.append({'name': f'{gene_set} 气泡图', 'path': dotplot_path})
                
                # 保存结果数据
                results[gene_set] = subset.to_dict('records')
        
        # 保存完整结果到CSV
        csv_path = os.path.join(output_dir, 'enrichment_results.csv')
        enr.results.to_csv(csv_path, index=False)
        
        return {
            'success': True,
            'message': '富集分析完成',
            'results': results,
            'plot_paths': plot_paths,
            'csv_path': csv_path,
            'full_results': enr.results.to_dict('records')
        }
    
    except Exception as e:
        return {
            'success': False,
            'error': str(e)
        }

def get_tissue_cell_line_types(file_paths=None):
    """获取CSV文件中所有的tissue_cell_line类型"""
    import csv
    tissue_types = set()
    
    if not file_paths:
        return []
    
    for file_path in file_paths:
        if not file_path.endswith('.csv'):
            continue
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    tissue = row.get('tissue_cell_line', '').strip()
                    if tissue:
                        tissue_types.add(tissue)
        except Exception:
            continue
    
    return sorted(list(tissue_types))


def filter_csv_by_tissue(file_paths=None, tissue_cell_line=None):
    """
    根据tissue_cell_line筛选CSV文件数据
    :param file_paths: CSV文件路径列表
    :param tissue_cell_line: tissue cell line类型（None表示不筛选）
    :return: 筛选后的文件路径列表（保持不变）和筛选后的基因列表
    """
    import csv
    
    if not file_paths:
        return file_paths, []
    
    if not tissue_cell_line:
        return file_paths, []
    
    filtered_genes = []
    
    for file_path in file_paths:
        if not file_path.endswith('.csv'):
            continue
        
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                gene_col = detect_gene_column_from_header(reader.fieldnames)
                
                if gene_col is None:
                    continue
                
                for row in reader:
                    row_tissue = row.get('tissue_cell_line', '').strip()
                    
                    if tissue_cell_line and row_tissue != tissue_cell_line:
                        continue
                    
                    gene = row.get(gene_col, '').strip()
                    if gene:
                        filtered_genes.append(gene)
        except Exception:
            continue
    
    return file_paths, filtered_genes


def detect_gene_column_from_header(fieldnames):
    """从表头中检测基因列"""
    if not fieldnames:
        return None
    
    for col in fieldnames:
        col_lower = str(col).lower().strip()
        for keyword in GENE_COLUMN_KEYWORDS:
            if keyword in col_lower:
                return col
    return None


def get_available_gene_sets():
    """获取可用的基因集库列表"""
    try:
        # 获取Enrichr可用的基因集库
        libraries = gp.get_library_name()
        # 筛选常用的KEGG和GO相关库
        kegg_go_libs = [lib for lib in libraries if 'KEGG' in lib or 'GO' in lib]
        return sorted(kegg_go_libs)
    except Exception as e:
        # 如果获取失败，返回默认列表
        return [
            'KEGG_2021',
            'KEGG_2019',
            'GO_Biological_Process_2021',
            'GO_Cellular_Component_2021',
            'GO_Molecular_Function_2021',
            'GO_Biological_Process_2018',
            'GO_Cellular_Component_2018',
            'GO_Molecular_Function_2018'
        ]