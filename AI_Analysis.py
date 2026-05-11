import os
import csv
import io

try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False

DEFAULT_AI_SETTINGS = {
    'model': '',
    'api_key': '',  
    'api_base': 'https://api.deepseek.com'
}


def analyze_compound(file_contents, api_key=None, api_base=None, system_prompt=None, model=None):
    """
    分析化合物筛选结果
    :param file_contents: 文件内容列表，每个元素为 {'filename': 'xxx.csv', 'content': 'csv_content'}
    :param api_key: API密钥
    :param api_base: API基础URL
    :param system_prompt: 系统提示词
    :param model: AI模型
    :return: AI分析结果
    """
    if not OPENAI_AVAILABLE:
        return {"success": False, "error": "OpenAI SDK未安装，请运行: pip install openai"}
    

    api_key = api_key or DEFAULT_AI_SETTINGS['api_key']
    api_base = api_base or DEFAULT_AI_SETTINGS['api_base']
    model = model or DEFAULT_AI_SETTINGS['model']
    
    if not api_key:
        return {"success": False, "error": "请先在前端页面设置API密钥"}
    
    if not file_contents or len(file_contents) == 0:
        return {"success": False, "error": "未上传任何文件"}
    
    all_data = []
    for file_item in file_contents:
        filename = file_item['filename']
        content = file_item['content']
        
        f = io.StringIO(content)
        reader = csv.DictReader(f)
        for row in reader:
            all_data.append({
                'file': filename,
                'CID': row.get('CID', ''),
                'Name': row.get('Name', ''),
                'Ingredient_name': row.get('Ingredient_name', ''),
                'gene_name': row.get('gene', ''),
                'pathway': row.get('pathway', ''),
                'dose': row.get('dose', ''),
                'tissue_cell_line': row.get('tissue_cell_line', ''),
                'time_after_irradiation': row.get('time_after_irradiation', ''),
                'Drug_likeness': row.get('Drug_likeness', ''),
                'OB_score': row.get('OB_score', '')
            })

    if not all_data:
        return {"success": False, "error": "文件中没有数据"}
    
    input_text = "化合物筛选结果分析：\n\n"
    input_text += f"共找到 {len(all_data)} 条相关记录\n\n"
    
    herb_data = {}
    for item in all_data:
        try:
            herb_name = item['file'].split(' ')[4]
        except IndexError:
            herb_name = '未知药材'
        if herb_name not in herb_data:
            herb_data[herb_name] = []
        herb_data[herb_name].append(item)
    
    for herb_name, items in herb_data.items():
        input_text += f"【{herb_name}】\n"
        input_text += f"  相关化合物数量: {len(items)}\n"
        
        gene_count = {}
        pathway_count = {}
        drug_likeness_values = []
        ob_score_values = []
        compound_types = {}
        
        for item in items:
            if item['gene_name']:
                gene_count[item['gene_name']] = gene_count.get(item['gene_name'], 0) + 1
            if item['pathway']:
                pathway_count[item['pathway']] = pathway_count.get(item['pathway'], 0) + 1
            if item['Drug_likeness']:
                try:
                    drug_likeness_values.append(float(item['Drug_likeness']))
                except ValueError:
                    pass
            if item['OB_score']:
                try:
                    ob_score_values.append(float(item['OB_score']))
                except ValueError:
                    pass
            if item['Ingredient_name']:
                compound_types[item['Ingredient_name']] = compound_types.get(item['Ingredient_name'], 0) + 1
        
        top_genes = sorted(gene_count.items(), key=lambda x: x[1], reverse=True)[:]
        if top_genes:
            input_text += "  所有基因: " + ", ".join([f"{g} ({c}次)" for g, c in top_genes]) + "\n"
        
        top_pathways = sorted(pathway_count.items(), key=lambda x: x[1], reverse=True)[:]
        if top_pathways:
            input_text += "  所有通路: " + ", ".join([f"{p} ({c}次)" for p, c in top_pathways]) + "\n"
        
        if drug_likeness_values:
            avg_drug_likeness = sum(drug_likeness_values) / len(drug_likeness_values)
            input_text += f"  平均药物相似性: {avg_drug_likeness:.2f}\n"
        if ob_score_values:
            avg_ob_score = sum(ob_score_values) / len(ob_score_values)
            input_text += f"  平均口服生物利用度: {avg_ob_score:.2f}\n"
        
        top_compound_types = sorted(compound_types.items(), key=lambda x: x[1], reverse=True)[:]
        if top_compound_types:
            input_text += "  所有化合物类型: " + ", ".join([f"{t} ({c}次)" for t, c in top_compound_types]) + "\n"
        
        input_text += "\n"
    
    client = OpenAI(
        api_key=api_key,
        base_url=api_base
    )
    
    
    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": input_text}
    ]
    
    response = client.chat.completions.create(
        model=model,
        messages=messages,
        stream=False
    )
    
    result = response.choices[0].message.content
    
    return {"success": True, "result": result}