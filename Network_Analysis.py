import os
import csv
from openai import OpenAI
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
import requests
import json
import re
from matplotlib.patches import Patch

output_dir='Network Analysis Results'
os.makedirs(output_dir, exist_ok=True)

def get_LCMS_result(lcms_file_path=None):

    if lcms_file_path is None:
        return {"错误": "请上传文件"}

    if not lcms_file_path.endswith('.csv'):
        return {"错误": "文件格式错误，请上传csv后缀文件"}
   
    lcms = {}
    try:
        with open(lcms_file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                lcms_compound = row.get('Tentative Identification', '').strip()
                abundance = row.get('Abundance', '').strip()
                if lcms_compound:
                    try:
                        lcms[lcms_compound] = float(abundance)
                    except ValueError:
                        continue
    except FileNotFoundError:
        return {"错误": "文件不存在"}
    except Exception as e:
        return {"错误": str(e)} 
    return lcms

def get_tissue_cell_line_types(file_paths=None):
    tissue_types = set()
    
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


def filter_data_by_tissue(file_paths=None, tissue_cell_line=None):
    """
    :param file_paths: list of CSV file paths
    :param tissue_cell_line: tissue cell line type（None for all types）
    :return: filtered data list
    """
    if not file_paths:
        return []
    
    all_data = []
    for file_path in file_paths:
        if not file_path.endswith('.csv'):
            continue
        
        filename = os.path.basename(file_path)
        
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                row_tissue = row.get('tissue_cell_line', '')
                
                if tissue_cell_line and row_tissue != tissue_cell_line:
                    continue
                
                all_data.append({
                    'file': filename,
                    'CID': row.get('CID', ''),
                    'Name': row.get('Name', ''),
                    'gene_name': row.get('gene_name', ''),
                    'pathway': row.get('pathway', ''),
                    'tissue_cell_line': row_tissue,
                })
    
    return all_data


def fetch_ppi_from_stringdb(gene_list, MRI_score):
    '''
    :param gene_list: list of gene names
    :return: list of PPI edges, each element is a dict with keys 'gene1': str, 'gene2': str, 'score': float
    '''
    
    if not gene_list or len(gene_list) < 2:
        return []

    gene_list = list(set(gene_list))
    
    string_api_url = "https://string-db.org/api/json/network"
    
    payload = {
        "identifiers": "\n".join(gene_list[:2000]), 
        "species": 9606,
        "required_score": MRI_score,
        "network_type":"physical",
        "add_nodes":0,
        "show_query_node_labels":1
    }
    
    try:
        response = requests.post(
            string_api_url, 
            data=payload, 
            timeout=30,
            headers={"Accept":"application/json"}
            )
        if response.status_code == 200:
            data = response.json()
            edges = []

            for interaction in data:
                edges.append({
                    'gene1': interaction.get('preferredName_A', {}),
                    'gene2': interaction.get('preferredName_B', {}),
                    'score': interaction.get('score', 0)
                })
            return edges
    
    except requests.exceptions.Timeout:
        print("STRING API请求超时，请稍后重试")
    except requests.exceptions.RequestException as e:
        print(f"STRING API请求异常: {e}")
    except Exception as e:
        print(f"处理STRING响应时出错: {e}")
    
    
    return []


def create_ppi_network(data, MRI_score, figsize_w, figsize_h, dpi, use_string_db=True):

    output_path='Network Analysis Results/PPI_network.png'

    gene_list = []

    for item in data:
        gene = item.get('gene_name', '').strip()
        if gene:
            gene_list.append(gene)
    
    gene_list = list(set(gene_list))
    
    if len(gene_list) < 2:
        return {"success": False, "error": "基因数量不足，无法创建网络。请检查数据中的有效基因名。"}
    
    G = nx.Graph()
    G.add_nodes_from(gene_list)
    
    edges = []
    if use_string_db:
        edges = fetch_ppi_from_stringdb(gene_list, MRI_score)
        for edge in edges:
            if edge['score'] > MRI_score:
                G.add_edge(edge['gene1'], edge['gene2'], weight=edge['score'])
    
    if G.number_of_edges() == 0:
        return {"success": False, "error": f"不能创建有效蛋白互作网络。请检查数据。"}
    
    plt.figure(figsize=(figsize_w, figsize_h))
    
    degrees = dict(G.degree())
    
    node_sizes = [np.log1p(degrees[node]) * 100 + 80 for node in G.nodes()]

    pos = nx.spring_layout(G, k=3, iterations=80, seed=607)
    
    edges_list = list(G.edges())
    weights = [G[u][v].get('weight', 1) for u, v in edges_list]
    max_weight = max(weights)
    nx.draw_networkx_edges(G, pos, 
                           alpha=0.6, 
                           width=[w/max_weight*1 + 0.35 for w in weights],
                           edge_color='black')
    
    node_colors = [degrees[node] for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos,
                          node_size=node_sizes,
                          node_color=node_colors,
                          cmap=plt.cm.YlOrRd,
                          alpha=0.8)
    
    nx.draw_networkx_labels(G, pos, 
                          font_size=8,
                          font_weight='bold')
    
    plt.title('Protein-Protein Interaction Network', fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    network_stats = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'avg_degree': sum(degrees.values()) / len(degrees),
        'density': nx.density(G),
        'connected_components': nx.number_connected_components(G)
    }
    
    return {
        "success": True,
        "output_file": output_path,
        "network_stats": network_stats,
        "graph": G
    }

def keygene_analysis(G, score_threshold, degree_cutoff, max_depth):
    """
    keygene algorithm implementation - search key genes（core modules）
    :param G: NetworkX graph object
    :param score_threshold: score threshold
    :param degree_cutoff: degree degree cutoff
    :param max_depth: max depth
    :return: key genes list and scores list
    """
    
    if G.number_of_nodes() == 0:
        return []
    
    total_nodes = G.number_of_nodes()
    node_weights = {}

    try:
        core_numbers = nx.core_number(G)
    except Exception:
        core_numbers = dict(G.degree())


    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        k = len(neighbors)
        
        if k >= degree_cutoff:
            neighbor_edges = 0
            for i, n1 in enumerate(neighbors):
                for n2 in neighbors[i+1:]:
                    if G.has_edge(n1, n2):
                        neighbor_edges += 1

            if k > 1:
                density = neighbor_edges / (k * (k-1))
            else:
                density = 0
            
            core_weight = core_numbers.get(node, 0)
            node_weights[node] = density * (core_weight + 1)
        else:
            node_weights[node] = 0
    
    sorted_nodes = sorted(node_weights.items(), key=lambda x: x[1], reverse=True)

    clusters = []
    visited = set()

    avg_core_weight = sum(core_numbers.values()) / len(core_numbers)
    
    for node, weight in sorted_nodes:
        if weight < score_threshold and node in visited:
            continue
        
        cluster = {node}
        queue = [node]
        depth_map = {node: 0}

        while queue:
            current = queue.pop(0)
            
            if current in visited:
                continue
                
            visited.add(current)
            current_depth = depth_map.get(current,0)

            if current_depth >= max_depth:
                continue
                
            for neighbor in G.neighbors(current):
                if neighbor not in visited and neighbor not in cluster:
                    neighbor_weight = node_weights.get(neighbor, 0)

                    if neighbor_weight >= score_threshold:
                        cluster.add(neighbor)
                        queue.append(neighbor)
                        depth_map[neighbor] = current_depth + 1
        
        if len(cluster) > 2:
            clusters.append(cluster)

    keygene_results = []
    for cluster in clusters:
        cluster_subgraph = G.subgraph(cluster)
        k = len(cluster)
        
        if k > 1:
            edges_in_cluster = cluster_subgraph.number_of_edges()
            max_possible_edges = k * (k - 1) / 2
            density = edges_in_cluster / max_possible_edges
        else:
            density = 0

        avg_cluster_weight = sum(node_weights.get(node, 0) for node in cluster) / k
        
        if avg_core_weight > 0:
            keygene_score = density * (avg_cluster_weight / avg_core_weight)
        else:
            keygene_score = density
        
        size_bonus = min(1.5, 1 + (k / total_nodes))
        keygene_score = keygene_score * size_bonus

        keygene_results.append({
            'genes': list(cluster),
            'score': keygene_score,
            'size': k,
            'density': density
        })
    
    keygene_results.sort(key=lambda x: x['score'], reverse=True)
    
    return keygene_results

def analyze_ppi_and_keygene(file_paths=None, tissue_cell_line=None, output_dir='Network Analysis Results',
                          MRI_score=0.4, figsize_w=16, figsize_h=12, score_threshold=0.2,
                          degree_cutoff=2, max_depth=100, dpi=300):
    """
    main function: create PPI network and perform keygene analysis
    :param file_paths: list of CSV file paths
    :param tissue_cell_line: tissue type or cell line
    :param MRI_score: MRI score threshold
    :param figsize_w: figure width
    :param figsize_h: figure height
    :param dpi: DPI resolution
    :return: analysis results
    """

    data = filter_data_by_tissue(file_paths, tissue_cell_line)

    if not data:
        return {"success": False, "error": "未找到符合要求的数据"}

    ppi_output = os.path.join(output_dir, 'PPI_network.png')
    ppi_result = create_ppi_network(data, MRI_score, figsize_w, figsize_h, dpi)
    
    if not ppi_result.get('success'):
        return ppi_result
    
    G = ppi_result.get('graph')
    keygene_results = keygene_analysis(G, score_threshold, degree_cutoff, max_depth)
    
    key_genes = []
    for cluster in keygene_results:
        for gene in cluster['genes']:
            key_genes.append({
                'gene': gene,
                'keygene_score': cluster['score'],
                'cluster_size': cluster['size'],
                'cluster_rank': keygene_results.index(cluster) + 1
            })
    
    key_genes.sort(key=lambda x: x['keygene_score'], reverse=True)
    
    keygene_csv = os.path.join(output_dir, 'key_genes_keygene.csv')
    with open(keygene_csv, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['gene', 'keygene_score', 'cluster_size', 'cluster_rank'])
        writer.writeheader()
        writer.writerows(key_genes)
    
    return {
        "success": True,
        "ppi_network": ppi_output,
        "key_genes_csv": keygene_csv,
        "key_genes": key_genes,
        "keygene_clusters": keygene_results,
        "network_stats": ppi_result.get('network_stats', {})
    }


def create_gene_compound_network(data, key_genes, output_path='Network Analysis Results/gene_compound_network.png'):
    """
    create gene-compound interaction interaction network
    :param data: filtered data list
    :param key_genes: key genes list
    :param output_path: output path
    :return: operation result
    """
    
    if not key_genes:
        return {"success": False, "error": "未找到符合要求的基因数据"}
    
    G = nx.Graph()
    
    gene_set = set()
    compound_set = set()
    edges = []
    compound_gene_count = {}
    
    for item in data:
        gene = item.get('gene_name', '').strip()
        compound = item.get('Name', '').strip()
        
        if gene and compound:
            key_gene_names = [kg['gene'] for kg in key_genes]
            if gene in key_gene_names:
                gene_set.add(gene)
                compound_set.add(compound)
                edges.append((gene, compound))
                compound_gene_count = defaultdict(int)
                compound_gene_count[compound] += 1
    
    if not edges:
        return {"success": False, "error": "未找到符合要求的基因交互数据"}
    
    G.add_nodes_from(gene_set, node_type='gene')
    G.add_nodes_from(compound_set, node_type='compound')
    G.add_edges_from(edges)
    
    plt.figure(figsize=(18, 14))
    
    pos = {}
    
    gene_list = list(gene_set)
    compound_list = list(compound_set)
    
    for i, gene in enumerate(gene_list):
        pos[gene] = (0, i * 2 - len(gene_list))
    
    for i, compound in enumerate(compound_list):
        pos[compound] = (2, i * 2 - len(compound_list))
    
    nx.draw_networkx_edges(G, pos, alpha=0.5, edge_color='steelblue', width=1.5)
    
    gene_nodes = set(gene_list) & set(G.nodes())
    nx.draw_networkx_nodes(G, pos, 
                          nodelist=gene_nodes,
                          node_color='coral',
                          node_size=800,
                          alpha=0.9,
                          edgecolors='darkred',
                          linewidths=2)
    
    compound_nodes = set(compound_list) & set(G.nodes())
    nx.draw_networkx_nodes(G, pos, 
                          nodelist=compound_nodes,
                          node_color='lightgreen',
                          node_size=600,
                          alpha=0.9,
                          edgecolors='darkgreen',
                          linewidths=2)
    
    labels = {}
    for node in G.nodes():
        if len(node) > 20:
            labels[node] = node[:17] + '...'
        else:
            labels[node] = node
    
    nx.draw_networkx_labels(G, pos, labels, font_size=7)
    
    plt.title('Key Gene - Compound Interaction Network', fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return {
        "success": True,
        "output_file": output_path,
        "genes": list(gene_set),
        "compounds": list(compound_set),
        "interactions": len(edges),
        "compound_gene_count": dict(sorted(compound_gene_count.items(), key=lambda x: x[1], reverse=True))
    }


def full_network_analysis(file_paths=None, tissue_cell_line=None, output_dir='Network Analysis Results',
                         MRI_score=0.4, figsize_w=16, figsize_h=12, score_threshold=0.2,
                         degree_cutoff=2, max_depth=100, dpi=300):
    """
    full network analysis process
    :param file_paths: list of CSV file paths
    :param tissue_cell_line: tissue/cell line type
    :param MRI_score: MRI score threshold
    :param figsize_w: figure width
    :param figsize_h: figure height
    :param score_threshold: score threshold
    :param degree_cutoff: degree cutoff
    :param max_depth: max depth
    :param dpi: DPI resolution

    :return: full analysis result
    """
    ppi_result = analyze_ppi_and_keygene(file_paths, tissue_cell_line, output_dir,
                          MRI_score, figsize_w, figsize_h, score_threshold,
                          degree_cutoff, max_depth, dpi)
    
    if not ppi_result.get('success'):
        return ppi_result
    
    key_genes = ppi_result.get('key_genes', [])
    
    if not key_genes:
        return {"success": False, "error": "未找到核心基因"}
    
    data = filter_data_by_tissue(file_paths, tissue_cell_line)
    
    gene_compound_output = os.path.join(output_dir, 'gene_compound_network.png')
    gc_result = create_gene_compound_network(data, key_genes, gene_compound_output)
    
    return {
        "success": True,
        "ppi_network": os.path.basename(ppi_result.get('ppi_network')),
        "key_genes_csv": os.path.basename(ppi_result.get('key_genes_csv')),
        "gene_compound_network": os.path.basename(gc_result.get('output_file')),
        "key_genes": key_genes,
        "keygene_clusters": ppi_result.get('keygene_clusters', []),
        "network_stats": ppi_result.get('network_stats', {}),
        "gene_compound_stats": {
            "genes": gc_result.get('genes', []),
            "compounds": gc_result.get('compounds', []),
            "interactions": gc_result.get('interactions', 0)
        }
    }



def ai_select_key_genes(data, keygene_results=None, api_key=None, api_base=None, model=None, system_prompt=None):
    """
    use AI to select key genes from gene list
    :param data: filtered data list
    :param keygene_results: keygene analysis results for reference
    :param api_key: API key
    :param api_base: API base URL
    :param model: AI model
    :param system_prompt: system prompt
    :return: AI selected key genes list
    """
    
    api_key = api_key or DEFAULT_AI_SETTINGS['api_key']
    api_base = api_base or DEFAULT_AI_SETTINGS['api_base']
    model = model or DEFAULT_AI_SETTINGS['model']
    
    if not api_key:
        return {"success": False, "error": "Please set API key first before using AI"}
    
    gene_info = {}
    for item in data:
        gene = item.get('gene_name', '').strip()
        if not gene:
            continue
        if gene not in gene_info:
            gene_info[gene] = {
                'compounds': set(),
                'pathways': set(),
                'tissues': set()
            }
        gene_info[gene]['compounds'].add(item.get('Ingredient_name', ''))
        gene_info[gene]['pathways'].add(item.get('pathway', ''))
        gene_info[gene]['tissues'].add(item.get('tissue_cell_line', ''))
    
    if not gene_info:
        return {"success": False, "error": "No gene data found"}
    
    input_text = "请分析以下基因列表，根据生物学意义选择最关键的核心基因。\n\n"
    
    if keygene_results and isinstance(keygene_results, list) and len(keygene_results) > 0:
        input_text += "【关键基因分析结果参考】\n"
        input_text += "以下是通过模拟MCODE算法分析得到的核心基因簇：\n"
        for i, cluster in enumerate(keygene_results):
            genes = cluster.get('genes', [])
            score = cluster.get('score', 0)
            size = cluster.get('size', 0)
            density = cluster.get('density', 0)
            input_text += f"簇 {i+1}: 基因={', '.join(genes)}, 得分={score:.4f}, 大小={size}, 密度={density:.4f}\n"
        input_text += "\n"
    
    input_text += "【基因详细信息】\n"
    
    for gene, info in gene_info.items():
        input_text += f"基因: {gene}\n"
        input_text += f"  关联化合物: {', '.join(info['compounds']) if info['compounds'] else '无'}\n"
        input_text += f"  关联通路: {', '.join(info['pathways']) if info['pathways'] else '无'}\n"


        input_text += "\n"
    
    input_text += """
请根据以下标准选择关键基因:
1. 与多个化合物相关联（可能是多靶点药物）
2. 参与重要信号通路
请以JSON格式返回选择的关键基因，格式如下:
{"selected_genes": [{"gene": "基因名", "reason": "选择原因", "importance_score": 0-10}], "summary": "总体分析摘要"}
"""
    
    default_prompt = system_prompt or "你是一位专业的生物信息学家，擅长基因网络分析和药物靶点识别。"
    
    client = OpenAI(api_key=api_key, base_url=api_base)
    
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": default_prompt},
                {"role": "user", "content": input_text}
            ],
            temperature=0.3
        )
        
        result = response.choices[0].message.content
        
        json_match = re.search(r'\{[\s\S]*\}', result)
        if json_match:
            parsed = json.loads(json_match.group())
            return {
                "success": True,
                "selected_genes": parsed.get('selected_genes', []),
                "summary": parsed.get('summary', ''),
                "raw_response": result
            }
        else:
            return {
                "success": True,
                "selected_genes": [],
                "summary": result,
                "raw_response": result
            }
    except Exception as e:
        return {"success": False, "error": f"AI分析失败: {str(e)}"}


def ai_predict_interactions(selected_genes, ppi_edges=None, api_key=None, api_base=None, model=None, system_prompt=None):
    """
    Use AI to predict interactions between genes based on selected genes
    :param selected_genes: AI selected key genes list
    :param ppi_edges: PPI edges from STRING DB for reference
    :param api_key: API key
    :param api_base: API base URL
    :param model: AI model
    :param system_prompt: system prompt
    :return: predicted interactions between genes
    """
    
    api_key = api_key or DEFAULT_AI_SETTINGS['api_key']
    api_base = api_base or DEFAULT_AI_SETTINGS['api_base']
    model = model or DEFAULT_AI_SETTINGS['model']
    
    if not api_key:
        return {"success": False, "error": "Please set API key first before using AI"}
    
    if not selected_genes:
        return {"success": False, "error": "No selected genes provided"}
    
    gene_names = []
    if isinstance(selected_genes, list):
        for g in selected_genes:
            if isinstance(g, dict):
                gene_names.append(g.get('gene', ''))
            else:
                gene_names.append(str(g))
    else:
        gene_names = [str(selected_genes)]
    
    gene_names = [g for g in gene_names if g]
    
    if len(gene_names) < 2:
        return {"success": False, "error": "At least 2 genes are required to predict interactions"}
    
    input_text = "请分析以下基因列表，预测它们之间可能的蛋白质相互作用关系。\n\n"
    
    if ppi_edges and isinstance(ppi_edges, list) and len(ppi_edges) > 0:
        input_text += "【STRING数据库PPI参考数据】\n"
        input_text += "以下是从STRING数据库获取的已知蛋白质相互作用（基因1-基因2，置信度得分）：\n"
        for edge in ppi_edges:
            gene1 = edge.get('gene1', '')
            gene2 = edge.get('gene2', '')
            score = edge.get('score', 0)
            if gene1 and gene2:
                input_text += f"{gene1} - {gene2} (置信度: {score:.4f})\n"
        input_text += "\n"
    
    input_text += "【待分析基因列表】\n"
    input_text += "基因列表: " + ", ".join(gene_names) + "\n\n"
    input_text += """
请根据以下原则预测潜在的相互作用:
1. 已知的蛋白质复合物成员
2. 同一信号通路中的上下游蛋白
3. 功能相似的蛋白（可能形成异源二聚体）
4. 调控关系（转录因子-靶基因）

请以JSON格式返回预测的相互作用，格式如下:
{"interactions": [{"gene1": "基因1", "gene2": "基因2", "interaction_type": "类型", "confidence": 0-10, "mechanism": "作用机制描述"}], "summary": "总体分析"}
"""
    
    default_prompt = system_prompt or "你是一位专业的系统生物学家，擅长预测蛋白质相互作用和信号通路。"
    
    client = OpenAI(api_key=api_key, base_url=api_base)
    
    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": default_prompt},
                {"role": "user", "content": input_text}
            ],
            temperature=0.3
        )
        
        result = response.choices[0].message.content
        
        
        
        json_match = re.search(r'\{[\s\S]*\}', result)
        if json_match:
            parsed = json.loads(json_match.group())
            return {
                "success": True,
                "interactions": parsed.get('interactions', []),
                "summary": parsed.get('summary', ''),
                "raw_response": result
            }
        else:
            return {
                "success": True,
                "interactions": [],
                "summary": result,
                "raw_response": result
            }
    except Exception as e:
        return {"success": False, "error": f"AI prediction failed: {str(e)}"}


def create_ai_enhanced_ppi_network(data, ai_interactions, output_path='Network Analysis Results/PPI_network_AI.png', MRI_score=0.4, figsize_w=16, figsize_h=12, dpi=300):
    """
    Create enhanced PPI network based on AI predicted interactions
    :param data: filtered data list
    :param ai_interactions: AI predicted interactions
    :param output_path: output file path
    :param MRI_score: score threshold for STRING PPI edges
    :return: network analysis result
    """
    
    gene_list = []
    for item in data:
        gene = item.get('gene_name', '').strip()
        if gene:
            gene_list.append(gene)
    gene_list = list(set(gene_list))
    
    if len(gene_list) < 2:
        return {"success": False, "error": "At least 2 genes are required to create a network"}
    
    G = nx.Graph()
    G.add_nodes_from(gene_list)
    
    string_edges_added = 0
    string_edges = fetch_ppi_from_stringdb(gene_list, MRI_score)
    for edge in string_edges:
        if edge['score'] > MRI_score:
            G.add_edge(edge['gene1'], edge['gene2'], weight=edge['score'], source='STRING')
            string_edges_added += 1
    
    ai_edges_added = 0
    if ai_interactions.get('success'):
        for interaction in ai_interactions.get('interactions', []):
            g1 = interaction.get('gene1', '')
            g2 = interaction.get('gene2', '')
            confidence = interaction.get('confidence', 5)
            interaction_type = interaction.get('interaction_type', 'unknown')
            mechanism = interaction.get('mechanism', '')

            if g1 in gene_list and g2 in gene_list:
                if G.has_edge(g1, g2):
                    current_weight = G[g1][g2].get('weight', 0)
                    current_source = G[g1][g2].get('source', '')
                    G[g1][g2]['weight'] = current_weight + confidence
                    G[g1][g2]['source'] = f"{current_source}+AI"
                    G[g1][g2]['interaction_type'] = interaction_type
                    G[g1][g2]['mechanism'] = mechanism
                else:
                    G.add_edge(g1, g2, weight=confidence, source='AI', 
                    interaction_type=interaction_type, mechanism=mechanism)
                    ai_edges_added += 1
    
    if G.number_of_edges() == 0:
        return {"success": False, "error": "No valid PPI network could be created with the given data and AI predictions"}
    
    plt.figure(figsize=(figsize_w, figsize_h))
    
    degrees = dict(G.degree())
    node_sizes = [np.log1p(degrees[node]) * 100 + 80 for node in G.nodes()]
    
    pos = nx.spring_layout(G, k=3, iterations=80, seed=607)
    
    string_edges = [(u, v) for u, v, d in G.edges(data=True) if 'STRING' in d.get('source', '')]
    ai_only_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get('source') == 'AI']
    ai_enhanced_edges = [(u, v) for u, v, d in G.edges(data=True) 
                         if 'AI' in d.get('source', '') and 'STRING' in d.get('source', '')]

    if string_edges:
        nx.draw_networkx_edges(G, pos, edgelist=string_edges, 
                               alpha=0.6, width=1.35, edge_color='gray')
    
    if ai_only_edges:
        nx.draw_networkx_edges(G, pos, edgelist=ai_only_edges, 
                               alpha=0.8, width=1.75, edge_color='red')
    
    if ai_enhanced_edges:
        nx.draw_networkx_edges(G, pos, edgelist=ai_enhanced_edges, 
                               alpha=1, width=2.25, edge_color='green')
    
    node_colors = [degrees[node] for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, 
                          node_size=node_sizes,
                          node_color=node_colors,
                          cmap=plt.cm.YlOrRd,
                          alpha=0.8)

    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')

    plt.title(f'AI-Enhanced PPI Network\n'
              f'(Gray dashed: STRING DB, Red: AI predicted, Green: STRING+AI combined)\n'
              f'Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}', 
              fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    network_stats = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'string_edges': string_edges_added,
        'ai_predicted_edges': ai_edges_added,
        'ai_enhanced_edges': len(ai_enhanced_edges), 
        'avg_degree': sum(degrees.values()) / len(degrees),
        'density': nx.density(G)
    }
    
    return {
        "success": True,
        "output_file": output_path,
        "network_stats": network_stats,
        "graph": G
    }

def create_compound_network_from_ppi(data, ppi_graph, figsize_w, figsize_h, dpi, output_dir='Network Analysis Results'):
    """
    :param data: Original data list with gene_name and Name columns
    :param ppi_graph: PPI network graph object
    :return: results of compound network analysis
    """
    
    gene_to_names = {}
    name_to_genes = {}
    
    for item in data:
        gene = item.get('gene_name', '').strip()
        name = item.get('Name', '').strip()
        if gene and name:
            if gene not in gene_to_names:
                gene_to_names[gene] = set()
            gene_to_names[gene].add(name)
            
            if name not in name_to_genes:
                name_to_genes[name] = set()
            name_to_genes[name].add(gene)
    
    compound_G = nx.Graph()
    
    all_compounds = set()
    for names in gene_to_names.values():
        all_compounds.update(names)
    compound_G.add_nodes_from(all_compounds)
    
    edges_added = set()
    
    for gene1, gene2 in ppi_graph.edges():
        compounds1 = gene_to_names.get(gene1, set())
        compounds2 = gene_to_names.get(gene2, set())

        ppi_source = 'unknown'
        if ppi_graph.has_edge(gene1, gene2):
            ppi_source = ppi_graph[gene1][gene2].get('source', 'unknown')

        for c1 in compounds1:
            for c2 in compounds2:
                if c1 != c2:
                    edge_key = tuple(sorted([c1, c2]))
                    weight = 1
                    if ppi_graph.has_edge(gene1, gene2):
                        weight = ppi_graph[gene1][gene2].get('weight', 1)

                    if edge_key not in edges_added:
                        compound_G.add_edge(c1, c2, weight=weight, ppi_source=ppi_source)
                        edges_added.add(edge_key)
                    else:
                        existing_source = compound_G[c1][c2].get('ppi_source', '')
                        if ppi_source not in existing_source:
                            compound_G[c1][c2]['ppi_source'] = f"{existing_source}+{ppi_source}"
                        compound_G[c1][c2]['weight'] += weight
    
    if compound_G.number_of_edges() == 0:
        return {"success": False, "error": "无法创建化合物网络，请检查基因-化合物映射关系"}
    
    output_path = os.path.join(output_dir, 'compound_network_from_ppi.png')
    plt.figure(figsize=(figsize_w, figsize_h))
    
    degrees = dict(compound_G.degree())
    
    node_sizes = [np.log1p(degrees[node]) * 100 + 80 for node in compound_G.nodes()]
    
    pos = nx.spring_layout(compound_G, k=3, iterations=0, seed=607)
    
    edges_list = list(compound_G.edges())
    if edges_list:
        weights = [compound_G[u][v].get('weight', 1) for u, v in edges_list]
        max_weight = max(weights) if weights else 1
        
        string_edges = [(u, v) for u, v in edges_list 
                       if 'STRING' in compound_G[u][v].get('ppi_source', '') and 'AI' not in compound_G[u][v].get('ppi_source', '')]
        ai_edges = [(u, v) for u, v in edges_list 
                   if 'AI' in compound_G[u][v].get('ppi_source', '') and 'STRING' not in compound_G[u][v].get('ppi_source', '')]
        mixed_edges = [(u, v) for u, v in edges_list 
                      if 'STRING' in compound_G[u][v].get('ppi_source', '') and 'AI' in compound_G[u][v].get('ppi_source', '')]
        
        if string_edges:
            nx.draw_networkx_edges(compound_G, pos, edgelist=string_edges,
                                   alpha=0.5, width=1.35, edge_color='gray', style='dashed')
        if ai_edges:
            nx.draw_networkx_edges(compound_G, pos, edgelist=ai_edges,
                                   alpha=0.7, width=1.75, edge_color='red')
        if mixed_edges:
            nx.draw_networkx_edges(compound_G, pos, edgelist=mixed_edges,
                                   alpha=0.9, width=2.25, edge_color='green')
    
    node_colors = [degrees[node] for node in compound_G.nodes()]
    nx.draw_networkx_nodes(compound_G, pos,
                          node_size=node_sizes,
                          node_color=node_colors,
                          cmap=plt.cm.Greens,
                          alpha=0.8)
    
    labels = {}
    for node in compound_G.nodes():
        if len(node) > 15:
            labels[node] = node[:15] + '...'
        else:
            labels[node] = node
    
    nx.draw_networkx_labels(compound_G, pos, labels, font_size=8, font_weight='bold')

    plt.title(f'Compound Interaction Network (Mapped from AI-Enhanced PPI)\n'
              f'Total Nodes: {compound_G.number_of_nodes()}, Edges: {compound_G.number_of_edges()}\n'
              f'(Gray: STRING-based, Red: AI-predicted, Green: Mixed)',
              fontsize=14, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    network_stats = {
        'nodes': compound_G.number_of_nodes(),
        'edges': compound_G.number_of_edges(),
        'avg_degree': sum(degrees.values()) / len(degrees) if degrees else 0,
        'density': nx.density(compound_G),
        'connected_components': nx.number_connected_components(compound_G)
    }
    
    return {
        "success": True,
        "output_file": output_path,
        "network_stats": network_stats,
        "graph": compound_G,
        "mapping_stats": {
            'unique_genes': len(gene_to_names),
            'unique_compounds': len(all_compounds),
            'ppi_edges_mapped': len(edges_added)
        }
    }

def analyze_compound_centrality(compound_graph, output_dir='Network Analysis Results'):
    """
    Analyze compound centrality metrics and identify core compounds.
    :param compound_graph: Compound network graph (NetworkX graph object)
    :param output_dir: Output directory
    :return: Dictionary with centrality metrics and core compounds identified
    """
    if compound_graph.number_of_nodes() == 0:
        return {"success": False, "error": "化合物网络为空"}
    
    results = []
    # Degree Centrality
    dc = nx.degree_centrality(compound_graph)
    # Betweenness Centrality
    bc = nx.betweenness_centrality(compound_graph, weight='weight')
    # Closeness Centrality
    cc = nx.closeness_centrality(compound_graph, distance='weight')
    # Eigenvector Centrality
    ec = nx.eigenvector_centrality(compound_graph, weight='weight', max_iter=1000)
    
    for node in compound_graph.nodes():
        node_data = {
            'compound': node,
            'degree': compound_graph.degree(node),
            'degree_centrality': dc.get(node, 0),
            'betweenness_centrality': bc.get(node, 0),
            'closeness_centrality': cc.get(node, 0),
            'eigenvector_centrality': ec.get(node, 0),
        }
        results.append(node_data)
    
    dc_values = [r['degree_centrality'] for r in results]
    bc_values = [r['betweenness_centrality'] for r in results]
    cc_values = [r['closeness_centrality'] for r in results]
    ec_values = [r['eigenvector_centrality'] for r in results]
    
    max_dc = max(dc_values) if dc_values else 0
    max_bc = max(bc_values) if bc_values else 0
    max_cc = max(cc_values) if cc_values else 0
    max_ec = max(ec_values) if ec_values else 0
    
    for r in results:
        # Normalize to 0-1 range
        norm_dc = r['degree_centrality'] / max_dc if max_dc > 0 else 0
        norm_bc = r['betweenness_centrality'] / max_bc if max_bc > 0 else 0
        norm_cc = r['closeness_centrality'] / max_cc if max_cc > 0 else 0
        norm_ec = r['eigenvector_centrality'] / max_ec if max_ec > 0 else 0
        
        # Composite score
        r['normalized_dc'] = norm_dc
        r['normalized_bc'] = norm_bc
        r['normalized_cc'] = norm_cc
        r['normalized_ec'] = norm_ec
        
        # Composite score = DC + BC + CC + EC
        r['composite_score'] = (norm_dc + norm_bc + norm_cc + norm_ec)/4
    
    # Sort by composite score in descending order
    results.sort(key=lambda x: x['composite_score'], reverse=True)
    
    # Save CSV results
    csv_path = os.path.join(output_dir, 'compound_centrality_analysis.csv')
    with open(csv_path, 'w', encoding='utf-8-sig', newline='') as f:
        fieldnames = ['compound', 'degree', 'degree_centrality', 'betweenness_centrality',
                     'closeness_centrality', 'eigenvector_centrality',
                     'normalized_dc', 'normalized_bc', 'normalized_cc', 
                     'normalized_ec', 'composite_score']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    # Identify core compounds
    # Method 1：Top 10% or top 10 compounds
    top_n = max(5, int(len(results) * 0.1))
    key_compounds = results[:top_n]
    
    # Method 2：Top 20% compounds in all metrics
    dc_threshold = sorted(dc_values, reverse=False)[max(1, int(len(dc_values) * 0.2))]
    bc_threshold = sorted(bc_values, reverse=False)[max(1, int(len(bc_values) * 0.2))]
    cc_threshold = sorted(cc_values, reverse=False)[max(1, int(len(cc_values) * 0.2))]
    ec_threshold = sorted(ec_values, reverse=False)[max(1, int(len(ec_values) * 0.2))]
    
    high_performance = []
    for r in results:
        if (r['degree_centrality'] >= dc_threshold and 
            r['betweenness_centrality'] >= bc_threshold and
            r['closeness_centrality'] >= cc_threshold):
            high_performance.append(r)
    
    # Save key compounds list
    key_path = os.path.join(output_dir, 'key_compounds.csv')
    with open(key_path, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['compound', 'composite_score', 'degree_centrality',
                                               'betweenness_centrality', 'closeness_centrality',
                                               'degree'])
        writer.writeheader()
        for compound in key_compounds:
            writer.writerow({
                'compound': compound['compound'],
                'composite_score': compound['composite_score'],
                'degree_centrality': compound['degree_centrality'],
                'betweenness_centrality': compound['betweenness_centrality'],
                'closeness_centrality': compound['closeness_centrality'],
                'degree': compound['degree']
            })
    
    # Network statistics summary
    network_stats = {
        'total_nodes': compound_graph.number_of_nodes(),
        'total_edges': compound_graph.number_of_edges(),
        'density': nx.density(compound_graph),
        'avg_degree': sum(dict(compound_graph.degree()).values()) / compound_graph.number_of_nodes(),
        'avg_clustering': nx.average_clustering(compound_graph),
        'num_connected_components': nx.number_connected_components(compound_graph),
        'largest_component_size': max(len(c) for c in nx.connected_components(compound_graph)) if compound_graph.number_of_nodes() > 0 else 0
    }
    
    return {
        "success": True,
        "all_compounds_analysis": results,
        "key_compounds": key_compounds,
        "high_performance_compounds": high_performance,
        "network_stats": network_stats,
        "csv_file": csv_path,
        "key_file": key_path
    }

def analyze_corestone_compounds(gc_result, lcms, output_dir='Network Analysis Results'):
    """
    Analyze corestone compounds
    :param gc_result: gene-compound network analysis result
    :param lcms: herb LC-MS list
    :param output_dir: output directory for CSV file
    :return: corestone compounds list and CSV file
    """
    # get common compounds
    common_compound = set(lcms.keys()) & set(gc_result['compound_gene_count'].keys())
    
    if not common_compound:
        return {"success": False, "错误": "未找到共同化合物"}
    
    compound_data = []
    for compound in common_compound:
        compound_dict = {
            'compound': compound,
            'gene_count': gc_result['compound_gene_count'][compound],
            'abundance': lcms[compound]
        }
        compound_data.append(compound_dict)

    value_gene = [c['gene_count'] for c in compound_data]
    value_abundance = [c['abundance'] for c in compound_data]

    max_gene = max(value_gene)if value_gene else 0
    sum_abundance = sum(value_abundance) if value_abundance else 0
    for c in compound_data:
        # normalize gene count and abundance
        norm_gene = c['gene_count']/max_gene if max_gene > 0 else 0
        norm_abundance = c['abundance']/sum_abundance if sum_abundance > 0 else 0

        c['normalized_gene'] = norm_gene
        c['normalized_abundance'] = norm_abundance
        c['corestone_score'] = (norm_gene + norm_abundance)/2

    compound_data.sort(key=lambda x: x['corestone_score'], reverse=True)

    # save CSV result
    csv_path = os.path.join(output_dir, 'corestone_compounds.csv')
    with open(csv_path, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['compound', 'gene_count', 'abundance', 'normalized_gene', 'normalized_abundance', 'corestone_score'])
        writer.writeheader()
        for compound in compound_data:
            writer.writerow(compound)

    # sort common compounds by gene count
    top_n = max(5, int(len(compound_data) * 0.1))
    corestone_compounds = compound_data[:top_n]



    return {
        "success": True,
        "corestone_compounds": corestone_compounds,
    }

def visualize_key_network(compound_graph, key_compounds, 
                                   output_path='Network Analysis Results/key_network.png',
                                   figsize_w=16, figsize_h=12, dpi=300):
    """
    Visualize key compounds network (highlight core compounds)
    :param compound_graph: Compound network graph
    :param key_compounds: Key compounds list
    :param output_path: Output path
    :param figsize_w: Image width
    :param figsize_h: Image height
    :param dpi: Image resolution
    :return: results
    """
    
    if compound_graph.number_of_nodes() == 0:
        return {"success": False, "error": "化合物网络为空"}
    
    plt.figure(figsize=(figsize_w, figsize_h))
    
    key_names = set([c['compound'] for c in key_compounds])
    
    pos = nx.spring_layout(compound_graph, k=4, iterations=80, seed=42)
    
    # 分离关键化合物和其他化合物
    key_nodes = [node for node in compound_graph.nodes() if node in key_names]
    other_nodes = [node for node in compound_graph.nodes() if node not in key_names]
    
    nx.draw_networkx_edges(compound_graph, pos, alpha=0.3, edge_color='gray', width=1)
    
    if other_nodes:
        nx.draw_networkx_nodes(compound_graph, pos, 
                              nodelist=other_nodes,
                              node_color='lightblue',
                              node_size=300,
                              alpha=0.6)
    
    if key_nodes:
        node_sizes = []
        for node in key_nodes:
            for c in key_compounds:
                if c['compound'] == node:
                    size = 500 + c['composite_score'] * 1000
                    node_sizes.append(size)
                    break
        
        nx.draw_networkx_nodes(compound_graph, pos,
                              nodelist=key_nodes,
                              node_color='red',
                              node_size=node_sizes if node_sizes else 800,
                              alpha=0.9,
                              edgecolors='darkred',
                              linewidths=2)
    
    labels = {}
    for node in compound_graph.nodes():
        if node in key_names:
            # Key compounds show full name
            labels[node] = node
        else:
            # Other compounds show first 12 characters only
            labels[node] = node if len(node) <= 15 else node[:12] + '...'
    
    nx.draw_networkx_labels(compound_graph, pos, labels, font_size=8)
    
    # Add legend
    legend_elements = [
        Patch(facecolor='red', alpha=0.9, label=f'key Compounds ({len(key_nodes)})'),
        Patch(facecolor='lightblue', alpha=0.6, label=f'Other Compounds ({len(other_nodes)})')
    ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=10)
    
    plt.title(f'Compound Network with key Compounds\n'
              f'Total: {compound_graph.number_of_nodes()} nodes, {compound_graph.number_of_edges()} edges\n'
              f'key: {len(key_nodes)} compounds', 
              fontsize=14, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    
    return {
        "success": True,
        "output_file": output_path,
        "num_key": len(key_nodes),
        "num_other": len(other_nodes)
    }

def ai_network_analysis(file_paths=None, lcms_file_path=None, tissue_cell_line=None,
                         api_key=None, api_base=None, model=None, system_prompt=None,
                         output_dir='Network Analysis Results',
                         MRI_score=0.4, figsize_w=16, figsize_h=12,
                         score_threshold=0.2, degree_cutoff=2, max_depth=100, dpi=300):
    """
    Complete AI enhanced analysis workflow for PPI networks
    :param file_paths: list of CSV file paths
    :param tissue_cell_line: Tissue or cell line type
    :param api_key: API key
    :param api_base: API base URL
    :param model: AI model
    :param system_prompt: System prompt
    :param lcms_file_path: LC-MS data file path
    :return: Complete analysis results
    """

    data = filter_data_by_tissue(file_paths, tissue_cell_line)
    if lcms_file_path is not None:
        lcms = get_LCMS_result(lcms_file_path)
    
    if not data:
        return {"success": False, "error": "No valid data found for the given tissue or cell line type"}
    
    gene_list = []
    for item in data:
        gene = item.get('gene_name', '').strip()
        if gene:
            gene_list.append(gene)
    gene_list = list(set(gene_list))
    
    ppi_edges = fetch_ppi_from_stringdb(gene_list, MRI_score)
    
    ppi_result = create_ppi_network(data, MRI_score, figsize_w, figsize_h, dpi)
    if not ppi_result.get('success'):
        return ppi_result
    
    G = ppi_result.get('graph')
    keygene_results = keygene_analysis(G, score_threshold, degree_cutoff, max_depth)
    
    ai_genes_result = ai_select_key_genes(data, keygene_results, api_key, api_base, model, system_prompt)
    
    if not ai_genes_result.get('success'):
        return ai_genes_result
    
    selected_genes = ai_genes_result.get('selected_genes', [])
    
    if not selected_genes:
        return {"success": False, "error": "AI failed to identify any key genes"}
    
    ai_interactions = ai_predict_interactions(selected_genes, ppi_edges, api_key, api_base, model, system_prompt)
    
    ppi_output = os.path.join(output_dir, 'PPI_network_AI.png')
    ppi_result = create_ai_enhanced_ppi_network(data, ai_interactions, ppi_output,
        MRI_score=MRI_score,
        figsize_w=figsize_w,
        figsize_h=figsize_h,
        dpi=dpi)
    
    if not ppi_result.get('success'):
        return ppi_result
    
    key_genes = []
    for i, g in enumerate(selected_genes):
        if isinstance(g, dict):
            key_genes.append({
                'gene': g.get('gene', ''),
                'reason': g.get('reason', ''),
                'importance_score': g.get('importance_score', 0),
                'ai_rank': i + 1
            })
        else:
            key_genes.append({
                'gene': str(g),
                'reason': '',
                'importance_score': 0,
                'ai_rank': i + 1
            })
    
    keygene_csv = os.path.join(output_dir, 'key_genes_AI.csv')
    with open(keygene_csv, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['gene', 'reason', 'importance_score', 'ai_rank'])
        writer.writeheader()
        writer.writerows(key_genes)
    
    gene_compound_output = os.path.join(output_dir, 'gene_compound_network_AI.png')
    gc_result = create_gene_compound_network(data, key_genes, gene_compound_output)
    
    compound_network_result = None
    centrality_result = None
    key_viz_result = None
    corestone_compounds = None

    ai_ppi_graph = ppi_result.get('graph')
    
    if ai_ppi_graph and ai_ppi_graph.number_of_edges() > 0:
        compound_network_result = create_compound_network_from_ppi(
            data, 
            ai_ppi_graph,
            figsize_w, 
            figsize_h, 
            dpi, 
            output_dir
        )
        if compound_network_result.get('success'):
            compound_graph = compound_network_result.get('graph')
            if compound_graph and compound_graph.number_of_nodes() > 0:
                centrality_result = analyze_compound_centrality(compound_graph, output_dir)
                
                if centrality_result.get('success') and centrality_result.get('key_compounds'):
                    key_viz_result = visualize_key_network(
                        compound_graph,
                        centrality_result['key_compounds'],
                        output_path=os.path.join(output_dir, 'key_compound_network_AI.png'),
                        figsize_w=figsize_w,
                        figsize_h=figsize_h,
                        dpi=dpi
                    )
        if lcms_file_path is not None:
            corestone_compounds = analyze_corestone_compounds(gc_result, lcms, output_dir)

    return {
        "success": True,
        "ppi_network": os.path.basename(ppi_result.get('output_file')),
        "key_genes_csv": os.path.basename(keygene_csv),
        "gene_compound_network": os.path.basename(gc_result.get('output_file', '')),
        "compound_network_from_ppi": os.path.basename(compound_network_result.get('output_file')) if compound_network_result and compound_network_result.get('success') else None,
        "compound_network_stats": compound_network_result.get('network_stats') if compound_network_result and compound_network_result.get('success') else None,
        "compound_centrality_csv": os.path.basename(centrality_result.get('csv_file')) if centrality_result and centrality_result.get('success') else None,
        "key_compounds_csv": os.path.basename(centrality_result.get('key_file')) if centrality_result and centrality_result.get('success') else None,
        "key_network_viz": os.path.basename(key_viz_result.get('output_file')) if key_viz_result and key_viz_result.get('success') else None,
        "key_genes": key_genes,
        "keygene_clusters": keygene_results[:10] if keygene_results else [],
        "ai_gene_selection": os.path.basename(ai_genes_result.get('summary', '')),
        "ai_interaction_prediction": os.path.basename(ai_interactions.get('summary', '')),
        "network_stats": ppi_result.get('network_stats', {}),
        "interactions": ai_interactions.get('interactions', []),
        "key_compounds": centrality_result.get('key_compounds', []) if centrality_result else [],
        "high_performance_compounds": centrality_result.get('high_performance_compounds', []) if centrality_result else [],
        "corestone_compounds": corestone_compounds.get('corestone_compounds', []) if corestone_compounds else []
    }
