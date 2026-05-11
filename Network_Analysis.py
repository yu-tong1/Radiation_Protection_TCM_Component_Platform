import os
import csv
from openai import OpenAI
OPENAI_AVAILABLE = True
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
NETWORKX_AVAILABLE = True
import requests
REQUESTS_AVAILABLE = True
import json
import re

output_dir='Network Analysis Results'
os.makedirs(output_dir, exist_ok=True)

def get_tissue_cell_line_types(file_paths=None):
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
                    'Ingredient_name': row.get('Ingredient_name', ''),
                    'gene_name': row.get('gene_name', ''),
                    'pathway': row.get('pathway', ''),
                    'dose': row.get('dose', ''),
                    'tissue_cell_line': row_tissue,
                    'time_after_irradiation': row.get('time_after_irradiation', ''),
                    'Drug_likeness': row.get('Drug_likeness', ''),
                    'OB_score': row.get('OB_score', '')
                })
    
    return all_data


def fetch_ppi_from_stringdb(gene_list, species=9606):
    """
    :param gene_list: gene list
    :param species: species ID（9606=human）
    :return: interaction edge list
    """
    if not REQUESTS_AVAILABLE:
        return []
    
    if not gene_list or len(gene_list) < 2:
        return []
    
    string_api_url = "https://string-db.org/api/json/network"
    
    payload = {
        "genes": gene_list[:500], 
        "species": species
    }
    
    try:
        response = requests.post(string_api_url, json=payload, timeout=30)
        if response.status_code == 200:
            data = response.json()
            edges = []
            for interaction in data.get('interactions', []):
                edges.append({
                    'gene1': interaction.get('gene1', {}).get('alias'),
                    'gene2': interaction.get('gene2', {}).get('alias'),
                    'score': interaction.get('score', 0)
                })
            return edges
    except Exception as e:
        print(f"STRING API request failed: {e}")
    
    return []


def create_ppi_network(data, output_path='Network Analysis Results/PPI_network.png', use_string_db=True):
    """
    :param data: filtered data list
    :param output_path: output file path
    :param use_string_db: whether to use STRING database for true interactions
    :return: network analysis result
    """
    if not NETWORKX_AVAILABLE:
        return {"success": False, "error": "Please install networkx: pip install networkx matplotlib numpy"}
    
    gene_list = []
    for item in data:
        gene = item.get('gene_name', '').strip()
        if gene:
            gene_list.append(gene)
    
    gene_list = list(set(gene_list))
    
    if len(gene_list) < 2:
        return {"success": False, "error": "Insufficient gene count to create network. Please check the data for valid gene names."}
    
    G = nx.Graph()
    G.add_nodes_from(gene_list)
    
    edges = []
    if use_string_db and REQUESTS_AVAILABLE:
        edges = fetch_ppi_from_stringdb(gene_list)
        for edge in edges:
            if edge['score'] > 0.4:
                G.add_edge(edge['gene1'], edge['gene2'], weight=edge['score'])
    
    if G.number_of_edges() == 0:
        gene_cooccurrence = defaultdict(lambda: defaultdict(int))
        
        for item in data:
            gene = item.get('gene_name', '').strip()
            compound = item.get('Ingredient_name', '').strip()
            pathway = item.get('pathway', '').strip()
            
            if gene and compound:
                for other_gene in gene_list:
                    if other_gene != gene:
                        gene_cooccurrence[gene][other_gene] += 1
        
        for g1, neighbors in gene_cooccurrence.items():
            for g2, count in neighbors.items():
                if count >= 1:
                    G.add_edge(g1, g2, weight=min(count, 10))
    
    if G.number_of_edges() == 0:
        return {"success": False, "error": "can't create valid PPI network. Please check the data for valid PPI network."}
    
    plt.figure(figsize=(16, 12))
    
    degrees = dict(G.degree())
    node_sizes = [degrees[node] * 100 + 200 for node in G.nodes()]
    
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    edges_list = list(G.edges())
    weights = [G[u][v].get('weight', 1) for u, v in edges_list]
    max_weight = max(weights) if weights else 1
    
    nx.draw_networkx_edges(G, pos, 
                           alpha=0.6, 
                           width=[w/max_weight*3 + 0.5 for w in weights],
                           edge_color='gray')
    
    node_colors = [degrees[node] for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, 
                          node_size=node_sizes,
                          node_color=node_colors,
                          cmap=plt.cm.YlOrRd,
                          alpha=0.8)
    
    nx.draw_networkx_labels(G, pos, 
                          font_size=8,
                          font_weight='bold')
    
    plt.title('PPI Network (Protein-Protein Interaction)', fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
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


def mcode_analysis(G, score_threshold=0.2, degree_cutoff=2, max_depth=100):
    """
    MCODE algorithm implementation - search key genes（core modules）
    :param G: NetworkX graph object
    :param score_threshold: score threshold
    :param degree_cutoff: degree degree cutoff
    :param max_depth: max depth
    :return: key genes list and scores list
    """
    if not NETWORKX_AVAILABLE:
        return []
    
    if G.number_of_nodes() == 0:
        return []
    
    node_weights = {}
    for node in G.nodes():
        neighbors = list(G.neighbors(node))
        if len(neighbors) >= degree_cutoff:
            neighbor_edges = 0
            for i, n1 in enumerate(neighbors):
                for n2 in neighbors[i+1:]:
                    if G.has_edge(n1, n2):
                        neighbor_edges += 1
            k = len(neighbors)
            if k > 1:
                density = neighbor_edges / (k * (k-1) / 2)
                node_weights[node] = density
        else:
            node_weights[node] = 0
    
    clusters = []
    visited = set()
    
    sorted_nodes = sorted(node_weights.items(), key=lambda x: x[1], reverse=True)
    
    for node, weight in sorted_nodes:
        if weight < score_threshold or node in visited:
            continue
        
        cluster = {node}
        queue = [node]
        
        while queue:
            current = queue.pop(0)
            if current in visited:
                continue
            
            visited.add(current)
            
            for neighbor in G.neighbors(current):
                if neighbor not in visited and node_weights.get(neighbor, 0) >= score_threshold:
                    cluster.add(neighbor)
                    queue.append(neighbor)
        
        if len(cluster) >= 2:
            clusters.append(cluster)
    
    mcode_results = []
    for cluster in clusters:
        cluster_subgraph = G.subgraph(cluster)
        k = len(cluster)
        
        if k > 1:
            edges_in_cluster = cluster_subgraph.number_of_edges()
            density = edges_in_cluster / (k * (k-1) / 2)
            mcode_score = density * (k / G.number_of_nodes()) ** 0.5
        else:
            mcode_score = 0
        
        mcode_results.append({
            'genes': list(cluster),
            'score': mcode_score,
            'size': k
        })
    
    mcode_results.sort(key=lambda x: x['score'], reverse=True)
    
    return mcode_results


def analyze_ppi_and_mcode(file_paths=None, tissue_cell_line=None, output_dir='Network Analysis Results'):
    """
    main function: create PPI network and perform MCODE analysis
    :param file_paths: list of CSV file paths
    :param tissue_cell_line: tissue type or cell line
    :return: analysis results
    """
    if not NETWORKX_AVAILABLE:
        return {"success": False, "error": "Please install required libraries first: pip install networkx matplotlib numpy requests"}
    
    data = filter_data_by_tissue(file_paths, tissue_cell_line)
    
    if not data:
        return {"success": False, "error": "No data found matching the criteria provided"}
    
    ppi_output = os.path.join(output_dir, 'PPI_network.png')
    ppi_result = create_ppi_network(data, ppi_output)
    
    if not ppi_result.get('success'):
        return ppi_result
    
    G = ppi_result.get('graph')
    mcode_results = mcode_analysis(G)
    
    key_genes = []
    for cluster in mcode_results:
        for gene in cluster['genes']:
            key_genes.append({
                'gene': gene,
                'mcode_score': cluster['score'],
                'cluster_size': cluster['size'],
                'cluster_rank': mcode_results.index(cluster) + 1
            })
    
    key_genes.sort(key=lambda x: x['mcode_score'], reverse=True)
    
    mcode_csv = os.path.join(output_dir, 'key_genes_mcode.csv')
    with open(mcode_csv, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['gene', 'mcode_score', 'cluster_size', 'cluster_rank'])
        writer.writeheader()
        writer.writerows(key_genes)
    
    return {
        "success": True,
        "ppi_network": ppi_output,
        "key_genes_csv": mcode_csv,
        "key_genes": key_genes,
        "mcode_clusters": mcode_results,
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
    if not NETWORKX_AVAILABLE:
        return {"success": False, "error": "Please install required libraries first: pip install networkx matplotlib"}
    
    if not key_genes:
        return {"success": False, "error": "No key genes data available"}
    
    G = nx.Graph()
    
    gene_set = set()
    compound_set = set()
    edges = []
    
    for item in data:
        gene = item.get('gene_name', '').strip()
        compound = item.get('Name', '').strip()
        
        if gene and compound:
            key_gene_names = [kg['gene'] for kg in key_genes]
            if gene in key_gene_names:
                gene_set.add(gene)
                compound_set.add(compound)
                edges.append((gene, compound))
    
    if not edges:
        return {"success": False, "error": "No key genes interaction data found for the provided criteria"}
    
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
        "interactions": len(edges)
    }


def full_network_analysis(file_paths=None, tissue_cell_line=None, output_dir='Network Analysis Results'):
    """
    full network analysis process
    :param file_paths: list of CSV file paths
    :param tissue_cell_line: tissue/cell line type

    :return: full analysis result
    """
    ppi_result = analyze_ppi_and_mcode(file_paths, tissue_cell_line)
    
    if not ppi_result.get('success'):
        return ppi_result
    
    key_genes = ppi_result.get('key_genes', [])
    
    if not key_genes:
        return {"success": False, "error": "No key genes found"}
    
    data = filter_data_by_tissue(file_paths, tissue_cell_line)
    
    gene_compound_output = os.path.join(output_dir, 'gene_compound_network.png')
    gc_result = create_gene_compound_network(data, key_genes, gene_compound_output)
    
    return {
        "success": True,
        "ppi_network": os.path.basename(ppi_result.get('ppi_network')),
        "key_genes_csv": os.path.basename(ppi_result.get('key_genes_csv')),
        "gene_compound_network": os.path.basename(gc_result.get('output_file')),
        "key_genes": key_genes,
        "mcode_clusters": ppi_result.get('mcode_clusters', []),
        "network_stats": ppi_result.get('network_stats', {}),
        "gene_compound_stats": {
            "genes": gc_result.get('genes', []),
            "compounds": gc_result.get('compounds', []),
            "interactions": gc_result.get('interactions', 0)
        }
    }



def ai_select_key_genes(data, api_key=None, api_base=None, model=None, system_prompt=None):
    """
    use AI to select key genes from gene list
    :param data: filtered data list
    :param api_key: API key
    :param api_base: API base URL
    :param model: AI model
    :param system_prompt: system prompt
    :return: AI selected key genes list
    """
    if not OPENAI_AVAILABLE:
        return {"success": False, "error": "OpenAI SDK not installed"}
    
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
    
    input_text = "请分析以下基因列表，根据生物学意义选择最关键的核心基因：\n\n"
    
    for gene, info in gene_info.items():
        input_text += f"基因: {gene}\n"
        input_text += f"  关联化合物: {', '.join(info['compounds']) if info['compounds'] else '无'}\n"
        input_text += f"  关联通路: {', '.join(info['pathways']) if info['pathways'] else '无'}\n"
        input_text += f"  关联组织: {', '.join(info['tissues']) if info['tissues'] else '无'}\n"
        input_text += "\n"
    
    input_text += """
请根据以下标准选择关键基因:
1. 与多个化合物相关联（可能是多靶点药物）
2. 参与重要信号通路
3. 在多个组织/细胞系中出现
4. 具有潜在的治疗意义

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


def ai_predict_interactions(selected_genes, api_key=None, api_base=None, model=None, system_prompt=None):
    """
    Use AI to predict interactions between genes based on selected genes
    :param selected_genes: AI selected key genes list
    :param api_key: API key
    :param api_base: API base URL
    :param model: AI model
    :param system_prompt: system prompt
    :return: predicted interactions between genes
    """
    if not OPENAI_AVAILABLE:
        return {"success": False, "error": "OpenAI SDK not installed"}
    
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
    
    input_text = "请分析以下基因列表，预测它们之间可能的蛋白质相互作用关系：\n\n"
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


def create_ai_enhanced_ppi_network(data, ai_interactions, output_path='Network Analysis Results/PPI_network_AI.png'):
    """
    Create enhanced PPI network based on AI predicted interactions
    :param data: filtered data list
    :param ai_interactions: AI predicted interactions
    :param output_path: output file path
    :return: network analysis result
    """
    if not NETWORKX_AVAILABLE:
        return {"success": False, "error": "Please install networkx: pip install networkx matplotlib"}
    
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
    
    gene_cooccurrence = defaultdict(lambda: defaultdict(int))
    for item in data:
        gene = item.get('gene_name', '').strip()
        compound = item.get('Ingredient_name', '').strip()
        if gene and compound:
            for other_gene in gene_list:
                if other_gene != gene:
                    gene_cooccurrence[gene][other_gene] += 1
    
    for g1, neighbors in gene_cooccurrence.items():
        for g2, count in neighbors.items():
            if count >= 1:
                G.add_edge(g1, g2, weight=min(count, 10), source='data')
    
    if ai_interactions.get('success'):
        for interaction in ai_interactions.get('interactions', []):
            g1 = interaction.get('gene1', '')
            g2 = interaction.get('gene2', '')
            confidence = interaction.get('confidence', 5)
            if g1 in gene_list and g2 in gene_list:
                if G.has_edge(g1, g2):
                    G[g1][g2]['weight'] += confidence
                    G[g1][g2]['source'] += '+AI'
                else:
                    G.add_edge(g1, g2, weight=confidence, source='AI')
    
    if G.number_of_edges() == 0:
        return {"success": False, "error": "No valid PPI network could be created with the given data and AI predictions"}
    
    plt.figure(figsize=(18, 14))
    
    degrees = dict(G.degree())
    node_sizes = [degrees[node] * 100 + 200 for node in G.nodes()]
    
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)
    
    data_edges = []
    ai_edges = []
    for u, v, data in G.edges(data=True):
        source = data.get('source', '')
        if 'AI' in source:
            ai_edges.append((u, v))
        else:
            data_edges.append((u, v))
    
    if data_edges:
        nx.draw_networkx_edges(G, pos, edgelist=data_edges, 
                               alpha=0.4, width=1.5, edge_color='gray')
    
    if ai_edges:
        nx.draw_networkx_edges(G, pos, edgelist=ai_edges, 
                               alpha=0.8, width=2.5, edge_color='red')
    
    node_colors = [degrees[node] for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, 
                          node_size=node_sizes,
                          node_color=node_colors,
                          cmap=plt.cm.YlOrRd,
                          alpha=0.8)
    
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold')
    
    plt.title('AI-Enhanced PPI Network\n(Red edges: AI-predicted interactions)', 
              fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    network_stats = {
        'nodes': G.number_of_nodes(),
        'edges': G.number_of_edges(),
        'data_edges': len(data_edges),
        'ai_predicted_edges': len(ai_edges),
        'avg_degree': sum(degrees.values()) / len(degrees),
        'density': nx.density(G)
    }
    
    return {
        "success": True,
        "output_file": output_path,
        "network_stats": network_stats,
        "graph": G
    }


def ai_network_analysis(file_paths=None, tissue_cell_line=None,
                         api_key=None, api_base=None, model=None, system_prompt=None,
                         output_dir='Network Analysis Results'):
    """
    Complete AI enhanced analysis workflow for PPI networks
    :param file_paths: list of CSV file paths
    :param tissue_cell_line: Tissue or cell line type
    :param api_key: API key
    :param api_base: API base URL
    :param model: AI model
    :param system_prompt: System prompt
    :return: Complete analysis results
    """

    data = filter_data_by_tissue(file_paths, tissue_cell_line)
    
    if not data:
        return {"success": False, "error": "No valid data found for the given tissue or cell line type"}
    
    ai_genes_result = ai_select_key_genes(data, api_key, api_base, model, system_prompt)
    
    if not ai_genes_result.get('success'):
        return ai_genes_result
    
    selected_genes = ai_genes_result.get('selected_genes', [])
    
    if not selected_genes:
        return {"success": False, "error": "AI failed to identify any key genes"}
    
    ai_interactions = ai_predict_interactions(selected_genes, api_key, api_base, model, system_prompt)
    
    ppi_output = os.path.join(output_dir, 'PPI_network_AI.png')
    ppi_result = create_ai_enhanced_ppi_network(data, ai_interactions, ppi_output)
    
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
    
    mcode_csv = os.path.join(output_dir, 'key_genes_AI.csv')
    with open(mcode_csv, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['gene', 'reason', 'importance_score', 'ai_rank'])
        writer.writeheader()
        writer.writerows(key_genes)
    
    gene_compound_output = os.path.join(output_dir, 'gene_compound_network_AI.png')
    gc_result = create_gene_compound_network(data, key_genes, gene_compound_output)
    
    return {
        "success": True,
        "ppi_network": os.path.basename(ppi_result.get('output_file')),
        "key_genes_csv": os.path.basename(mcode_csv),
        "gene_compound_network": os.path.basename(gc_result.get('output_file', '')),
        "key_genes": key_genes,
        "ai_gene_selection": os.path.basename(ai_genes_result.get('summary', '')),
        "ai_interaction_prediction": os.path.basename(ai_interactions.get('summary', '')),
        "network_stats": ppi_result.get('network_stats', {}),
        "interactions": ai_interactions.get('interactions', [])
    }
