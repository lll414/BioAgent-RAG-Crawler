import os
import re

# ================= 配置路径 =================
# 检查目录：上一步处理好的 cleaned 目录
SOURCE_DIR = "/Users/lll/Desktop/实习/黑森智动/RAG_Supplement/output_cleaned"
# 报告输出位置
REPORT_FILE = os.path.expanduser("~/Desktop/noise_report.txt")

# ================= 正则表达式定义 =================
PATTERNS = {
    "R_PROMPT": re.compile(r'^[>+\$]\s+'),          # R 语言命令行提示符 (> , + , $ )
    "R_INDEX": re.compile(r'^\[\d+\]\s+'),          # R 语言输出索引 ([1] )
    "ZSH_ERROR": re.compile(r'^zsh: command not found'), # 终端错误
    "PROGRESS_BAR": re.compile(r'\[=*>*\s*\]\s*\d+%|Download:\s+\d+/\d+'), # 进度条
    "COPYRIGHT": re.compile(r'Copyright © \d{4}', re.IGNORECASE), # 版权
    "BASE64_IMG": re.compile(r'data:image\/[^;]+;base64,[a-zA-Z0-9+/=]+'), # Base64图片
    "BOILERPLATE": re.compile(r'\[¶\]|↩|\(Link to this Section\)|Skip to content', re.IGNORECASE) # 导航残留
}

def preview_file(filepath, output_log):
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
            lines = content.split('\n')

        found_issues = False
        file_report = []

        # 1. 检查 Base64 (针对全文)
        b64_matches = PATTERNS["BASE64_IMG"].findall(content)
        if b64_matches:
            found_issues = True
            for m in b64_matches:
                file_report.append(f"  [Base64图片] 发现一段长编码 (长度: {len(m)} 字符) -> 将被移除")

        # 2. 检查特定字符 (针对全文)
        bp_matches = PATTERNS["BOILERPLATE"].findall(content)
        if bp_matches:
            found_issues = True
            # 去重显示
            unique_bp = set(bp_matches)
            file_report.append(f"  [导航符号] 发现符号: {unique_bp} -> 将被移除")

        # 3. 逐行检查 (针对行级噪音)
        for i, line in enumerate(lines):
            line_num = i + 1
            stripped = line.strip()
            
            if not stripped: continue

            # 检查 ZSH 错误
            if PATTERNS["ZSH_ERROR"].match(stripped):
                found_issues = True
                file_report.append(f"  [第{line_num}行][删除整行] 终端报错: {stripped}")
                continue

            # 检查进度条
            if PATTERNS["PROGRESS_BAR"].search(stripped):
                found_issues = True
                file_report.append(f"  [第{line_num}行][删除整行] 进度条/日志: {stripped}")
                continue

            # 检查版权
            if PATTERNS["COPYRIGHT"].search(stripped):
                found_issues = True
                file_report.append(f"  [第{line_num}行][删除整行] 版权信息: {stripped}")
                continue

            # 检查 R 提示符
            if PATTERNS["R_PROMPT"].match(line):
                found_issues = True
                new_line = PATTERNS["R_PROMPT"].sub('', line)
                file_report.append(f"  [第{line_num}行][净化代码] 原文: '{line.strip()}'  ->  修改后: '{new_line.strip()}'")
            
            # 检查 R 索引
            elif PATTERNS["R_INDEX"].match(line):
                found_issues = True
                new_line = PATTERNS["R_INDEX"].sub('', line)
                file_report.append(f"  [第{line_num}行][净化输出] 原文: '{line.strip()}'  ->  修改后: '{new_line.strip()}'")

        # 如果发现问题，写入报告
        if found_issues:
            header = f"\n📄 文件: {os.path.basename(filepath)}"
            print(header)
            for msg in file_report:
                print(msg)
            
            # 同时写入文件
            output_log.write(header + "\n")
            for msg in file_report:
                output_log.write(msg + "\n")

    except Exception as e:
        print(f"⚠️ 读取文件失败: {filepath} ({e})")

# ================= 主程序 =================
print(f"🚀 开始预检查...")
print(f"📂 扫描目录: {SOURCE_DIR}")

if not os.path.exists(SOURCE_DIR):
    print("❌ 目录不存在")
    exit(1)

with open(REPORT_FILE, 'w', encoding='utf-8') as log_file:
    log_file.write("=== RAG 数据清洗预检查报告 ===\n")
    
    for root, dirs, files in os.walk(SOURCE_DIR):
        for file in files:
            if file.endswith(".md"):
                preview_file(os.path.join(root, file), log_file)

print(f"\n✅ 检查完成！")
print(f"📜 完整报告已保存至桌面的: noise_report.txt")
