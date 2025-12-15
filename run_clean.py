import os
import re

# ================= 配置路径 =================
# 输入目录：上一步处理过链接的目录
SOURCE_DIR = "/Users/lll/Desktop/实习/黑森智动/RAG_Supplement/output_cleaned"
# 输出目录：最终产出的目录
DEST_DIR = "/Users/lll/Desktop/实习/黑森智动/RAG_Supplement/output_final"

# ================= 正则表达式 =================
PATTERNS = {
    "R_PROMPT": re.compile(r'^[>+\$]\s+'),          
    "R_INDEX": re.compile(r'^\[\d+\]\s+'),          
    "ZSH_ERROR": re.compile(r'^zsh: command not found'), 
    "PROGRESS_BAR": re.compile(r'\[=*>*\s*\]\s*\d+%|Download:\s+\d+/\d+'), 
    "COPYRIGHT": re.compile(r'Copyright © \d{4}', re.IGNORECASE), 
    "BASE64_IMG": re.compile(r'data:image\/[^;]+;base64,[a-zA-Z0-9+/=]+'), 
    "BOILERPLATE": re.compile(r'\[¶\]|↩|\(Link to this Section\)|Skip to content', re.IGNORECASE) 
}

def clean_content(content):
    # 1. 全文替换 (Base64 和 导航符)
    content = PATTERNS["BASE64_IMG"].sub('', content)
    content = PATTERNS["BOILERPLATE"].sub('', content)
    
    lines = content.split('\n')
    cleaned_lines = []
    
    # 2. 逐行清洗
    for line in lines:
        stripped = line.strip()
        
        # --- 删除整行 ---
        if PATTERNS["ZSH_ERROR"].match(stripped): continue
        if PATTERNS["PROGRESS_BAR"].search(stripped): continue
        if PATTERNS["COPYRIGHT"].search(stripped): continue
        
        # --- 修改行内容 ---
        # 去除 R 提示符 (> , +)
        if PATTERNS["R_PROMPT"].match(line):
            line = PATTERNS["R_PROMPT"].sub('', line)
            
        # 去除 R 索引 ([1])
        elif PATTERNS["R_INDEX"].match(line):
            line = PATTERNS["R_INDEX"].sub('', line)
            
        cleaned_lines.append(line)
    
    # 3. 重新组合并规范化空行 (将 >3 个换行变成 2 个)
    final_text = '\n'.join(cleaned_lines)
    final_text = re.sub(r'\n{3,}', '\n\n', final_text)
    
    return final_text

# ================= 主程序 =================
if not os.path.exists(SOURCE_DIR):
    print(f"❌ 找不到目录: {SOURCE_DIR}")
    exit(1)

print(f"🚀 开始执行最终清洗...")
print(f"📂 读取: {SOURCE_DIR}")
print(f"📂 写入: {DEST_DIR}")

count = 0
for root, dirs, files in os.walk(SOURCE_DIR):
    for file in files:
        if file.endswith(".md"):
            src_path = os.path.join(root, file)
            
            # 计算目标路径
            rel_path = os.path.relpath(src_path, SOURCE_DIR)
            dest_path = os.path.join(DEST_DIR, rel_path)
            os.makedirs(os.path.dirname(dest_path), exist_ok=True)
            
            try:
                with open(src_path, 'r', encoding='utf-8', errors='ignore') as f:
                    raw = f.read()
                
                cleaned = clean_content(raw)
                
                with open(dest_path, 'w', encoding='utf-8') as f:
                    f.write(cleaned)
                count += 1
            except Exception as e:
                print(f"⚠️ 错误: {file} - {e}")

print(f"--------------------------------------------------")
print(f"✅ 清洗完成！共处理 {count} 个文件。")
print(f"🎉 最终可用数据在: {DEST_DIR}")
