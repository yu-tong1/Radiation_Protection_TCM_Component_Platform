# 辐射防护中药成分筛选系统

## 1. 系统概述

本系统是一个基于 Flask 框架开发的辐射防护中药成分筛选平台，旨在帮助研究人员快速筛选具有辐射防护作用的中药成分和方剂。系统集成了多种功能模块，包括化合物筛选、方剂筛选、重金属诱导排泄预测等，为辐射防护研究提供了有力的工具支持。

### 1.1 技术架构

- **后端**：Flask 框架
- **前端**：HTML、CSS、JavaScript、Bootstrap
- **API 集成**：BATMAN2 API（靶点预测）、STING API（蛋白互作网络构建）

### 1.2 系统目录结构

```
website/
├── Compound screening results/    # 筛选结果存储目录
├── templates/                      # 前端模板
│   ├── web_database_files/         # 静态资源文件
│   ├── Compound_screening.html     # 化合物筛选页面
│   ├── index.html                  # 主页页面
│   ├── Herb_miRNA_data.html        # 中药 miRNA 数据页面
│   ├── Network_analysis.html       # 网络分析页面
│   ├── KEGG_GO.html                # KEGG GO 分析页面
│   └── AI_Analysis.html            # AI 分析页面
├── KEGG_GO_analysis.py  # KEGG GO 分析模块
├── Compound_screening_normal.py    # 正常模式下的化合物筛选
├── Network_analysis.py   # 网络分析模块
├── AI_Analysis.py                  # AI 分析模块
├── app.py                          # 主应用程序
└── Test-file Zisu.csv                # 测试用例文件
```

## 2. 功能模块

### 2.1 辐射防护化合物筛选 (RadioProtect Compound Screening)

**功能说明**：
- 上传中药成分 CSV 文件，筛选具有辐射防护作用的化合物
- 支持正常模式筛选
- 生成匹配和未匹配的结果文件
- 提供结果下载功能

**技术实现**：
- 使用 `Compound_screening_normal.py` 处理筛选
- 通过 pandas 进行数据处理和匹配
- 调用 BATMAN2 API 进行靶点预测

### 2.2 网络分析 (Network Analysis)

**功能说明**：
- 上传2.1辐射防护化合物筛选运行结果CSV 文件，分析其之间的蛋白互作网络
- 进行蛋白互作网络分析，筛选核心靶点及其对应的辐射防护化合物
- AI挖掘潜在的蛋白互作关系及核心靶点，筛选对应的辐射防护化合物，生成专业的分析报告。
- 提供API配置选项，支持多种AI模型

**技术实现**：
- 使用 `Network_analysis.py` 实现分析逻辑
- 通过 STING API 构建蛋白互作网络
- 筛选核心靶点及其对应的辐射防护化合物

### 2.3 KEGG GO 分析 (KEGG GO Analysis)

**功能说明**：
- 上传2.1辐射防护化合物筛选运行结果CSV 文件，进行KEGG GO 分析

**技术实现**：
- 使用 `KEGG_GO_analysis.py` 实现

### 2.4 中药 miRNA 数据 (Herb miRNA data)

**功能说明**：
- 展示中药 miRNA 相关数据

**技术实现**：
- 通过前端页面展示 miRNA 数据

### 2.5 AI 解读筛选结果 (AI Analysis of Screening Results)

**功能说明**：
- 智能分析化合物筛选结果，生成专业的分析报告
- 支持不同细胞组织类型的分析
- 提供 API 配置选项，支持多种 AI 模型
- 可自定义系统提示词，调整分析方向

**技术实现**：
- 使用 `AI_Analysis.py` 实现后端分析逻辑
- 通过前端页面 `AI_Analysis.html` 提供用户交互界面
- 集成 OpenAI SDK，调用 AI 模型进行分析
- 从 `Compound screening results` 目录读取筛选结果文件
- 按药材分组，统计主要基因、通路、化合物类型等信息
- 生成结构化的分析报告，包括化合物种类富集、主要基因通路、药物相似性和口服生物利用度、防护机制等方面

## 3. 操作指南

### 3.1 系统启动

1. 确保安装了所需的依赖包：
   ```
   pip install flask pandas requests openai numpy matplotlib seaborn gseapy networkx
   ```

   **注意**：如果使用 AI 分析功能，需要额外安装 openai 库。

2. 进入 website 目录，运行应用：
   ```
   python app.py
   ```

### 3.2 辐射防护化合物筛选操作步骤

1. 在功能区点击"RadioProtect Compound Screening"
2. 点击 "选择文件" 按钮，上传包含中药成分的 CSV 文件，或在搜索框中输入药材名称点击Search按钮
3. 点击 "Run screening" 按钮提交
4. 等待系统处理完成
5. 查看筛选结果，包括匹配和未匹配的化合物
6. 点击 "Download" 按钮下载结果文件(或在文件夹Compound screening results下查看)

### 3.3 KEGG GO 分析操作步骤

1. 上传2.1辐射防护化合物筛选运行结果CSV 文件
2. 选择细胞或组织类型
3. 选择KEGG GO 分析类型、物种
4. 查看分析结果

### 3.4 网络分析操作步骤

1. 上传2.1辐射防护化合物筛选运行结果CSV 文件
2. 选择细胞或组织类型
3. 选择网络分析类型
4. 选择AI模型（可选）
5. 点击 "Run network analysis" 按钮提交

### 3.5 中药 miRNA 数据操作步骤

1. 在功能区点击"Plant microRNA Targets in Human"
2. 进入外部链接页面进行 miRNA 靶点预测

### 3.6 AI 解读筛选结果操作步骤

1. 上传2.1辐射防护化合物筛选运行结果CSV 文件
2. 选择细胞或组织类型
3. 配置 API 设置：
   - 输入 API Key
   - 设置 API Base URL（默认为 https://api.deepseek.com）
   - 选择 AI 模型
4. 可根据需要修改系统提示词，调整分析方向
5. 点击 "Run AI Analysis" 按钮开始分析
6. 等待系统处理完成，查看 AI 生成的分析报告



## 4. 数据格式要求

### 4.1 辐射防护化合物筛选数据格式

CSV 文件应包含以下列（任意一种）：
- CID 列：`cid`, `CID`, `Cid`, `cID`, `pubchem_id`, `PubChem_id`, `PubChem_ID`
- 基因名列：`gene_name`, `Gene_name`

**示例文件**：`Test-file Zisu.csv`


## 5. 结果解释

### 5.1 辐射防护化合物筛选结果

- **Matching results**：包含与辐射防护靶点匹配的化合物
- **Unmapped results**：未匹配的化合物

结果文件保存在 `Compound screening results` 目录中，命名格式为：`Matching normal results of [文件名] ingredients.csv`

### 5.2 KEGG GO 分析结果

- **KEGG析结果**：包含辐射防护化合物的 KEGG分析结果，包括通路名称、通路 ID、富集 p-value 等信息
- **GO 分析结果**：包含辐射防护化合物的 GO 分析结果，包括基因名称、基因 ID、GO ID、GO 分类、GO 分类 p-value 等信息

### 5.3 网络分析结果

- **网络分析结果**：包含辐射防护化合物的网络分析结果，包括蛋白互作网络、核心基因、核心基因-中药化合物关系等
- **AI 模型结果**：包含 AI 模型的分析结果，AI挖掘潜在蛋白互作关系、评分筛选核心基因、核心基因-中药化合物关系等

## 6. 系统配置与维护

### 6.1 配置文件

- **app.py**：主配置文件，包含路由和核心功能
- **BATMAN2_API_URL**：BATMAN2 API 地址（默认为 'http://batman2api.cloudna.cn/queryTarget'）

### 6.2 数据文件

- **DEGs_normal.csv**：正常模式下的差异表达基因
- **herb_ingredient.csv**：中药成分数据
- **results of target matching normal.csv**：正常模式下的靶点匹配结果

## 7. 常见问题与解决方案

### 7.1 上传文件失败

- **问题**：上传文件后系统提示错误
- **解决方案**：检查文件格式是否正确，确保文件包含所需的列（如 CID 或 SMILE）

### 7.2 筛选结果为空

- **问题**：筛选后没有匹配的化合物
- **解决方案**：检查输入文件中的 CID 或基因名是否正确，确保数据格式符合要求

### 7.3 系统运行缓慢

- **问题**：处理大型文件时系统响应缓慢
- **解决方案**：减小文件大小，或分批处理数据

### 7.4 API 调用失败

- **问题**：BATMAN2 API 调用失败
- **解决方案**：检查网络连接，确保 API 地址正确

## 8. 示例使用流程

### 8.1 使用测试文件进行辐射防护化合物筛选

1. 启动系统后，进入 "RadioProtect Compound Screening" 页面
2. 选择 "normal" 模式
3. 上传 "Test-file Zisu.csv"
4. 点击 "Submit" 按钮
5. 查看筛选结果，下载 "Matching normal results of Test-file Zisu ingredients.csv"

## 9. 技术支持

- **联系人**：马家骅教授 Email: jiahuama@swust.edu.cn
- **开发者**：
            马家骅（博士，教授，西南科技大学生命科学与农林学院）；
            柒世龙（硕士在读，西南科技大学生命科学与农林学院）；
            周阳（本科在读，西南科技大学生命科学与农林学院）。
- **机构**：西南科技大学生命科学与农林学院
- **项目**：辐射防护中药成分筛选系统