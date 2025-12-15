一、 上游原始数据处理 (Upstream Data Processing)
这是最重要、最基础的缺失部分。您当前的文档几乎都是从“表达矩阵”或“比对好的文件”开始的，但如何从测序仪下机的原始数据走到这一步，是所有分析的第一公里。
	•	1. 原始数据质控 (Raw Data QC)
	◦	知识点: 如何评估原始测序数据（FASTQ文件）的质量。
	◦	核心工具: FastQC, MultiQC
	◦	重要性: 这是决定后续分析成败的第一步。Agent如果不知道什么是Phred score、接头污染、GC含量偏离，就无法回答关于数据质量的问题。
	•	2. 数据预处理 (Data Pre-processing)
	◦	知识点: 去除低质量碱基和测序接头。
	◦	核心工具: Trimmomatic, Cutadapt, fastp
	◦	重要性: “垃圾进，垃圾出”，这是保证分析结果可靠性的关键。
	•	3. 序列比对 (Sequence Alignment/Mapping)
	◦	知识点: 将处理过的短序列（reads）比对到参考基因组上。
	◦	核心工具:
	▪	RNA-seq: STAR, HISAT2
	▪	DNA-seq / ChIP-seq: BWA, Bowtie2
	◦	重要性: 这是连接原始序列和生物学位置的桥梁，产生了BAM/SAM文件，是下游分析的基础。
	•	4. 定量 (Quantification)
	◦	知识点: 统计比对到每个基因/转录本上的reads数量，生成“表达矩阵”。
	◦	核心工具: featureCounts (基于比对), Salmon, Kallisto (免比对/轻量比对)
	◦	重要性: 直接生成您当前知识库中大多数流程的输入文件。

二、 其他主流“组学”领域
您的知识库以基因组和转录组为中心，但“多组学”是当前的主流。
	•	1. 蛋白质组学 (Proteomics)
	◦	知识点: 基于质谱（Mass Spectrometry）的数据分析，包括肽段鉴定、蛋白定量、翻译后修饰（PTM）分析。
	◦	核心工具/概念: MaxQuant, Proteome Discoverer, TMT/iTRAQ标记定量，Label-free定量。
	◦	重要性: 连接转录本与生物学功能的关键一环。
	•	2. 代谢组学 (Metabolomics)
	◦	知识点: 分析生物体内的所有小分子代谢物。包括峰检测、鉴定和差异分析。
	◦	核心工具/平台: XCMS, MetaboAnalyst
	◦	重要性: 更贴近表型，是系统生物学的重要组成部分。

三、 结构生物信息学 (Structural Bioinformatics)
这部分完全缺失，但它是一个非常经典且再次变得热门的领域。
	•	知识点: 研究蛋白质、RNA等大分子的三维结构与功能的关系。
	•	核心工具/概念:
	◦	蛋白质结构预测: AlphaFold, RoseTTAFold
	◦	分子对接 (Docking): 药物设计的基础，AutoDock Vina
	◦	分子动力学模拟 (MD Simulation): GROMACS, AMBER
	•	重要性: 在药物研发、蛋白质工程等领域不可或缺。

四、 基础IT技能与编程语言
您的知识库是R语言的“深水区”，但对于初学者和处理上游数据的用户来说，更基础的技能是必须的。
	•	1. 命令行 (Shell/Bash)
	◦	知识点: Linux/macOS终端的基本操作，使用grep, awk, sed等工具对大型文本文件（如VCF, BED, GTF）进行快速处理。
	◦	重要性: 上游分析工具几乎全部是命令行工具，这是生信分析师的“基本功”。
	•	2. Python在生物信息中的应用
	◦	知识点: 使用Python进行数据处理和分析。
	◦	核心库: Biopython (处理序列文件), Pandas (数据处理), Matplotlib/Seaborn (绘图)。
	◦	重要性: Python是生信领域与R并驾齐驱的语言，尤其在流程搭建和机器学习方面更具优势。
	•	3. 环境与包管理
	◦	知识点: 如何创建可复现的、隔离的分析环境。
	◦	核心工具: Conda / Bioconda, Docker / Singularity
	◦	重要性: 解决软件安装、版本冲突的终极方案，是“可重复性研究”的基石。
