import pandas as pd
import os
import csv
import requests

output_dir = "Compound screening results"
os.makedirs(output_dir, exist_ok=True)

BATMAN2_API_URL = 'http://batman2api.cloudna.cn/queryTarget'


def safe_get_column(df, column_names):

    if isinstance(column_names, str):
        column_names = [column_names]

    clean_columns = {}
    for col in df.columns:
        if isinstance(col, str):
            cleaned = col.lstrip('\ufeff').strip()
            clean_columns[cleaned] = col
        else:
            clean_columns[str(col)] = col

    for col in column_names:
        if col in clean_columns:
            original_col = clean_columns[col]
            return df[original_col].tolist(), original_col

    lower_names = [name.lower() for name in column_names]
    for cleaned, original_col in clean_columns.items():
        if cleaned.lower() in lower_names:
            return df[original_col].tolist(), original_col

    for cleaned, original_col in clean_columns.items():
        for name in column_names:
            if name.lower() in cleaned.lower():
                return df[original_col].tolist(), original_col

    return [], None

#Match and output the results.
def match_pathways(file_name1, file_name2, file_name3, file_name4, original_filename):
    file1 = pd.read_csv(file_name1)
    file2 = pd.read_csv(file_name2, low_memory=False)
    file3 = pd.read_csv(file_name3)
    file4 = pd.read_csv(file_name4)

    cid1, used_col1 = safe_get_column(file1, ['cid', 'CID', 'Cid', 'cID', 'pubchem_id', 'PubChem_id', 'PubChem_ID'])
    target1, used_col2 = safe_get_column(file1, ['gene_name', 'Gene_name'])

    base_name = os.path.splitext(original_filename)[0]
    matched1_out = os.path.join(output_dir, "Matching normal results of {} ingredients.csv".format(base_name))
    matched2_out = os.path.join(output_dir, "Unmapped normal results of {} ingredients.csv".format(base_name))

    if used_col2 is None:
        if used_col1 is None:
            raise ValueError(
                "No usable column found in the uploaded CSV. "
                "Expected a CID column (cid/CID/PubChem_id) or a gene_name column (gene_name/Gene_name). "
                "Please check your CSV headers."
            )

        matched1 = pd.merge(file1, file2, left_on=used_col1, right_on='cid0', how='inner')

        matched_cid = matched1[used_col1].unique()
        unmap_origin = file1[~file1[used_col1].isin(matched_cid)]

        matched2 = pd.merge(unmap_origin, file3, left_on=used_col1, right_on='PubChem_id', how='left')
        matched2 = matched2.drop('PubChem_id', axis=1)

        matched1.to_csv(matched1_out, index=False, encoding='utf-8-sig')
        matched2.to_csv(matched2_out, index=False, encoding='utf-8-sig')

    else:
        matched1 = pd.merge(file1, file4, left_on=used_col2, right_on='gene', how='inner')
        matched_tar = matched1[used_col2].unique()
        unmap_origin = file1[~file1[used_col2].isin(matched_tar)]

        if used_col1 is not None and used_col1 in unmap_origin.columns:
            matched2 = pd.merge(unmap_origin, file3, left_on=used_col1, right_on='PubChem_id', how='left')
            matched2 = matched2.drop('PubChem_id', axis=1)
        else:
            matched2 = unmap_origin

        if used_col1 is not None and used_col1 in matched1.columns:
            matched1 = pd.merge(matched1, file3, left_on=used_col1, right_on='PubChem_id', how='left')

        matched1.to_csv(matched1_out, index=False, encoding='utf-8-sig')
        matched2.to_csv(matched2_out, index=False, encoding='utf-8-sig')

    return matched1, matched2


def process_tcm_file_normal(tcm_file_path, original_filename):
    matched, unmatched = match_pathways(
        tcm_file_path,
        'results of target matching normal.csv',
        'herb_ingredient.csv',
        'DEGs_normal.csv',
        original_filename
    )
    return matched, unmatched


def handle_file_upload(file_path, original_filename):
    """处理文件上传筛选流程"""
    try:
        matched, unmatched = process_tcm_file_normal(file_path, original_filename)

        base_name = os.path.splitext(original_filename)[0]
        matched_file = "Matching normal results of " + base_name + " ingredients.csv"
        unmatched_file = "Unmapped normal results of " + base_name + " ingredients.csv"

        matched_data = matched.fillna('').to_dict('records') if hasattr(matched, 'to_dict') else matched
        unmatched_data = unmatched.fillna('').to_dict('records') if hasattr(unmatched, 'to_dict') else unmatched

        return {
            'success': True,
            'matched': matched_data,
            'unmatched': unmatched_data,
            'matched_file': matched_file,
            'unmatched_file': unmatched_file
        }
    except Exception as e:
        return {'success': False, 'error': str(e)}


def handle_query_target(data, upload_folder='uploads'):
    """处理靶点查询流程"""
    try:
        content = data.get('content')
        if not content or not isinstance(content, list):
            return {'success': False, 'error': 'Invalid content format'}

        first_item = content[0]
        item_name = first_item.get('clusterName', 'unnamed')
        filename = item_name.replace(' ', '_') + '.csv'
        temp_file_path = os.path.join(upload_folder, filename)

        response = requests.post(
            BATMAN2_API_URL,
            json=data,
            headers={'Content-Type': 'application/json'},
            timeout=60
        )

        if response.status_code != 200:
            return {'success': False, 'error': 'API error: ' + str(response.status_code)}

        result = response.json()

        with open(temp_file_path, 'w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['CID', 'Name', 'gene_name', 'gene_id', 'score']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            processed_results = []
            if isinstance(result, list):
                for item in result:
                    cid = item.get('cid', '')
                    name = item.get('name', '')
                    target_list = item.get('target', [])

                    if target_list:
                        for target_item in target_list:
                            gene_name = target_item.get('gene_name', '')
                            gene_id = target_item.get('gene_id', '')
                            score = target_item.get('score', '')
                            writer.writerow({'CID': cid, 'Name': name, 'gene_name': gene_name, 'gene_id': gene_id, 'score': score})
                            processed_results.append({
                                'name': name,
                                'cid': cid,
                                'gene_name': gene_name,
                                'gene_id': gene_id,
                                'score': score,
                            })
                    else:
                        writer.writerow({'CID': cid, 'Name': name, 'gene_name': '', 'gene_id': '', 'score': ''})
                        processed_results.append({
                            'name': name,
                            'cid': cid,
                            'gene_name': '',
                            'gene_id': '',
                            'score': '',
                        })

        try:
            matched, unmatched = process_tcm_file_normal(temp_file_path, filename)
        except Exception as e:
            return {'success': False, 'error': str(e)}
        finally:
            if os.path.exists(temp_file_path):
                try:
                    os.remove(temp_file_path)
                except:
                    pass

        base_name = os.path.splitext(filename)[0]
        matched_file = "Matching normal results of " + base_name + " ingredients.csv"
        unmatched_file = "Unmapped normal results of " + base_name + " ingredients.csv"

        matched_data = matched.fillna('').to_dict('records') if hasattr(matched, 'to_dict') else matched
        unmatched_data = unmatched.fillna('').to_dict('records') if hasattr(unmatched, 'to_dict') else unmatched

        return {
            'success': True,
            'results': processed_results,
            'matched': matched_data,
            'unmatched': unmatched_data,
            'matched_file': matched_file,
            'unmatched_file': unmatched_file
        }

    except requests.exceptions.Timeout:
        return {'success': False, 'error': 'API request timeout'}
    except requests.exceptions.RequestException as e:
        return {'success': False, 'error': 'API request failed: ' + str(e)}
    except Exception as e:
        return {'success': False, 'error': str(e)}