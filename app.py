# -*- coding: utf-8 -*-
from flask import Flask, request, jsonify, send_from_directory, render_template
import os
import sys
if getattr(sys, 'frozen', False):
    os.chdir(os.path.dirname(sys.executable))
import csv
import uuid
import requests
import webbrowser
import threading
import time
from Compound_screening_normal import handle_file_upload, handle_query_target
from Network_Analysis import (
    get_tissue_cell_line_types,
    full_network_analysis,
    ai_network_analysis
)
from AI_Analysis import analyze_compound
from KEGG_GO_Analysis import run_kegg_go_analysis, parse_csv_files, get_tissue_cell_line_types, filter_csv_by_tissue



app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs('Compound screening results', exist_ok=True)

BATMAN2_API_URL = 'http://batman2api.cloudna.cn/queryTarget'


@app.route('/web_database_files/<path:filename>')
def serve_static(filename):
    return send_from_directory('templates/web_database_files', filename)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/Compound_screening.html')
def compound_screening_normal():
    return render_template('Compound_screening.html')

@app.route('/Herb_miRNA_data.html')
def herb_miRNA_data():
    return render_template('Herb_miRNA_data.html')

@app.route('/Network_Analysis.html')
def network_analysis():
    return render_template('Network_Analysis.html')

@app.route('/AI_Analysis.html')
def AI_analysis():
    return render_template('AI_Analysis.html')

@app.route('/KEGG_GO.html')
def kegg_go_analysis():
    return render_template('KEGG_GO.html')

@app.route('/run', methods=['POST'])
def run_screening():
    if 'file' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'Empty filename'}), 400

    ext = os.path.splitext(file.filename)[1]
    temp_filename = str(uuid.uuid4()) + ext
    temp_path = os.path.join(UPLOAD_FOLDER, temp_filename)
    file.save(temp_path)

    original_filename = os.path.basename(file.filename)

    try:
        result = handle_file_upload(temp_path, original_filename)
        if result['success']:
            return jsonify({
                'matched': result['matched'],
                'unmatched': result['unmatched'],
                'matched_file': result['matched_file'],
                'unmatched_file': result['unmatched_file']
            })
        else:
            return jsonify({'error': result['error']}), 500
    finally:
        os.remove(temp_path)

@app.route('/query_target', methods=['POST'])
def query_target():
    try:
        data = request.get_json()
        if not data or 'content' not in data:
            return jsonify({'error': 'Invalid request: missing content'}), 400

        result = handle_query_target(data, UPLOAD_FOLDER)

        if result['success']:
            return jsonify({
                'results': result['results'],
                'matched': result['matched'],
                'unmatched': result['unmatched'],
                'matched_file': result['matched_file'],
                'unmatched_file': result['unmatched_file']
            })
        else:
            if 'timeout' in result['error'].lower():
                return jsonify({'error': result['error']}), 504
            else:
                return jsonify({'error': result['error']}), 500

    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/analyze-compound', methods=['POST'])
def analyze_compound_api():
    try:
        data = request.get_json()
        if not data or 'files' not in data:
            return jsonify({'success': False, 'error': 'Missing files parameter'}), 400
        
        files = data['files']
        api_key = data.get('api_key')
        api_base = data.get('api_base')
        system_prompt = data.get('system_prompt')
        model = data.get('model')
        result = analyze_compound(files, api_key, api_base, system_prompt, model)
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/get-tissue-types', methods=['POST'])
def get_tissue_types():
    try:
        if 'files' not in request.files:
            return jsonify({'success': False, 'error': 'No files uploaded'}), 400
        
        files = request.files.getlist('files')
        if len(files) == 0:
            return jsonify({'success': False, 'error': 'No files selected'}), 400
        
        file_paths = []
        for file in files:
            if file.filename == '':
                continue
            ext = os.path.splitext(file.filename)[1]
            if ext.lower() not in ['.csv']:
                return jsonify({'success': False, 'error': f"不支持的文件格式: {file.filename}。请上传CSV文件。"}), 400
            
            temp_filename = str(uuid.uuid4()) + ext
            temp_path = os.path.join(UPLOAD_FOLDER, temp_filename)
            file.save(temp_path)
            file_paths.append(temp_path)
        
        tissue_types = get_tissue_cell_line_types(file_paths)
        
        # 清理临时文件
        for fp in file_paths:
            os.remove(fp)
        
        return jsonify({'success': True, 'tissue_types': tissue_types})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/network-analysis', methods=['POST'])
def network_analysis_api():
    try:
        if 'files' not in request.files:
            return jsonify({'success': False, 'error': 'No files uploaded'}), 400
        
        files = request.files.getlist('files')
        if len(files) == 0:
            return jsonify({'success': False, 'error': 'No files selected'}), 400
        
        # 保存上传的文件
        file_paths = []
        for file in files:
            if file.filename == '':
                continue
            ext = os.path.splitext(file.filename)[1]
            if ext.lower() not in ['.csv']:
                # 清理已保存的文件
                for fp in file_paths:
                    os.remove(fp)
                return jsonify({'success': False, 'error': f"不支持的文件格式: {file.filename}。请上传CSV文件。"}), 400
            
            temp_filename = str(uuid.uuid4()) + ext
            temp_path = os.path.join(UPLOAD_FOLDER, temp_filename)
            file.save(temp_path)
            file_paths.append(temp_path)
        
        # 获取筛选参数
        tissue_cell_line = request.form.get('tissue_cell_line')
        use_ai = request.form.get('use_ai', 'false').lower() == 'true'
        
        try:
            if use_ai:
                api_key = request.form.get('api_key')
                api_base = request.form.get('api_base')
                model = request.form.get('model')
                system_prompt = request.form.get('system_prompt')
                result = ai_network_analysis(
                    file_paths=file_paths,
                    tissue_cell_line=tissue_cell_line,
                    output_dir='Network Analysis Results',
                    api_key=api_key,
                    api_base=api_base,
                    model=model,
                    system_prompt=system_prompt
                )
            else:
                result = full_network_analysis(
                    file_paths=file_paths,
                    tissue_cell_line=tissue_cell_line,
                    output_dir='Network Analysis Results'
                )
        finally:
            # 清理临时文件
            for fp in file_paths:
                os.remove(fp)
        
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/network-image/<path:filename>')
def serve_network_image(filename):
    return send_from_directory('Network Analysis Results', filename)

@app.route('/download/<path:filename>')
def download_file(filename):
    return send_from_directory('Compound screening results', filename, as_attachment=True)

@app.route('/api/kegg-go-analysis', methods=['POST'])
def kegg_go_analysis_api():
    try:
        data = request.get_json()
        if not data or 'gene_list' not in data:
            return jsonify({'success': False, 'error': 'Missing gene_list parameter'}), 400
        
        gene_list = data['gene_list']
        gene_sets = data.get('gene_sets', ['KEGG_2021_Human', 'GO_Biological_Process_2021'])
        organism = data.get('organism', 'human')
        
        result = run_kegg_go_analysis(gene_list, gene_sets, organism)
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/kegg-go-get-tissue-types', methods=['POST'])
def kegg_go_get_tissue_types_api():
    try:
        if 'files' not in request.files:
            return jsonify({'success': False, 'error': 'No files uploaded'}), 400
        
        files = request.files.getlist('files')
        if len(files) == 0:
            return jsonify({'success': False, 'error': 'No files selected'}), 400
        
        file_paths = []
        for file in files:
            if file.filename == '':
                continue
            ext = os.path.splitext(file.filename)[1]
            if ext.lower() not in ['.csv']:
                return jsonify({'success': False, 'error': f"不支持的文件格式: {file.filename}。请上传CSV文件。"}), 400
            
            temp_filename = str(uuid.uuid4()) + ext
            temp_path = os.path.join(UPLOAD_FOLDER, temp_filename)
            file.save(temp_path)
            file_paths.append(temp_path)
        
        tissue_types = get_tissue_cell_line_types(file_paths)
        
        # 清理临时文件
        for fp in file_paths:
            os.remove(fp)
        
        return jsonify({'success': True, 'tissue_types': tissue_types})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/api/kegg-go-analysis-files', methods=['POST'])
def kegg_go_analysis_files_api():
    try:
        if 'files' not in request.files:
            return jsonify({'success': False, 'error': 'No files uploaded'}), 400
        
        files = request.files.getlist('files')
        if len(files) == 0:
            return jsonify({'success': False, 'error': 'No files selected'}), 400
        
        gene_sets = request.form.getlist('gene_sets')
        if len(gene_sets) == 0:
            gene_sets = ['KEGG_2021_Human', 'GO_Biological_Process_2021']
        
        organism = request.form.get('organism', 'human')
        tissue_cell_line = request.form.get('tissue_cell_line', '')
        
        # 保存上传的文件
        file_paths = []
        for file in files:
            if file.filename == '':
                continue
            ext = os.path.splitext(file.filename)[1]
            if ext.lower() not in ['.csv']:
                return jsonify({'success': False, 'error': f"不支持的文件格式: {file.filename}。请上传CSV文件。"}), 400
            
            temp_filename = str(uuid.uuid4()) + ext
            temp_path = os.path.join(UPLOAD_FOLDER, temp_filename)
            file.save(temp_path)
            file_paths.append(temp_path)
        
        # 根据tissue_cell_line筛选数据
        if tissue_cell_line:
            _, filtered_genes = filter_csv_by_tissue(file_paths, tissue_cell_line)
            if not filtered_genes:
                # 清理临时文件
                for fp in file_paths:
                    os.remove(fp)
                return jsonify({'success': False, 'error': f"未找到匹配 tissue_cell_line '{tissue_cell_line}' 的数据"}), 400
            
            gene_list = filtered_genes
            parse_result = {
                'success': True,
                'genes': gene_list,
                'parsed_info': [{'filename': f'Filtered by {tissue_cell_line}', 'gene_column': 'gene', 'gene_count': len(gene_list)}]
            }
        else:
            # 解析CSV文件
            parse_result = parse_csv_files(file_paths)
            if not parse_result['success']:
                # 清理临时文件
                for fp in file_paths:
                    os.remove(fp)
                return jsonify(parse_result), 400
            
            gene_list = parse_result['genes']
        
        # 执行富集分析
        result = run_kegg_go_analysis(gene_list, gene_sets, organism)
        
        # 添加解析信息到结果中
        result['parsed_info'] = parse_result['parsed_info']
        
        # 清理临时文件
        for fp in file_paths:
            os.remove(fp)
        
        return jsonify(result)
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

@app.route('/kegg-go-plot/<path:filename>')
def serve_kegg_go_plot(filename):
    import urllib.parse
    decoded_filename = urllib.parse.unquote(filename)
    # 处理完整路径情况，提取文件名
    decoded_filename = os.path.basename(decoded_filename)
    return send_from_directory('kegg_go_results', decoded_filename)

@app.route('/kegg-go-download/<path:filename>')
def download_kegg_go_result(filename):
    import urllib.parse
    decoded_filename = urllib.parse.unquote(filename)
    # 处理完整路径情况，提取文件名
    decoded_filename = os.path.basename(decoded_filename)
    return send_from_directory('kegg_go_results', decoded_filename, as_attachment=True)

if __name__ == '__main__':
    # 在服务器启动前打开浏览器
    
    def open_browser():
        time.sleep(1)  # 给服务器1秒启动时间
        webbrowser.open('http://localhost:6011')
    
    # 启动浏览器线程
    threading.Thread(target=open_browser).start()
    
    # 运行Flask应用
    app.run(debug=False, port=6011)