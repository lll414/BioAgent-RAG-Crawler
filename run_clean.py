import os
import re

# ================= 1. 配置区域 =================
# 输入目录：存放待处理的 .md 文件的文件夹
SOURCE_DIR = "/Users/lll/Desktop/实习/黑森智动/RAG_Supplement/output"
# 输出目录：用于存放清理后文件的全新文件夹
DEST_DIR = "/Users/lll/Desktop/实习/黑森智动/RAG_Supplement/RAG_Supplement_filtered"
# 报告输出位置：如果选择诊断模式，报告将保存在此
REPORT_FILE = os.path.join(os.path.dirname(DEST_DIR), "cleaning_report.txt")


# ================= 2. 噪音定义区域 (保持不变) =================
NOISE_KEYWORDS = [
    # 账户/操作类
    "Login", "Log in", "Sign in", "登录", "登入", "Logout", "Log out", "Sign out", "登出", "退出",
    "Sign up", "Register", "注册", "创建账户", "My Account", "我的账户", "个人中心", "Submit", "提交",
    "Search", "搜索", "Download", "下载", "Go",
    # 导航/页脚类
    "Home", "首页", "About", "About Us", "关于", "关于我们", "Contact", "Contact Us", "联系", "联系我们",
    "Products", "产品", "Services", "服务", "Blog", "博客", "News", "新闻", "FAQ", "常见问题", "Help", "帮助",
    "Support", "支持", "Terms of Service", "服务条款", "Privacy Policy", "隐私政策", "Sitemap", "网站地图",
    "Disclaimer", "免责声明",
    # 互动/社交类
    "Share", "分享", "Tweet", "转发", "Like", "点赞", "Comment", "评论", "Reply", "回复", "Edit", "编辑",
    "Share on Facebook", "Share to Twitter",
    # 网页通用模板
    "Skip to content", "跳到内容", "Back to top", "回到顶部", "Advertisement", "广告", "Sponsored", "赞助内容",
    # 交互式元素
    "Accept", "同意", "Decline", "拒绝", "Got it", "知道了", "Preferences", "设置"
]

removable_link_text = "|".join(re.escape(k) for k in NOISE_KEYWORDS)
INLINE_NOISE_PATTERN = re.compile(r'\[\s*(?:' + removable_link_text + r')\s*\]\([^)]*\)', re.IGNORECASE)
WHOLE_LINE_NOISE_PATTERNS = [re.compile(r'^\s*Copyright ©|© \d{4}.*All rights reserved', re.IGNORECASE)]


# ================= 3. 功能函数定义 =================

def clean_markdown_conservatively(content):
    """
    【清理函数】实际执行清理操作，返回干净的文本。
    """
    lines = content.split('\n')
    cleaned_lines = []
    
    # 阶段一：识别并跳过应被整行删除的行
    for line in lines:
        stripped_line = line.strip()
        is_line_noise = False
        
        for pattern in WHOLE_LINE_NOISE_PATTERNS:
            if pattern.search(stripped_line):
                is_line_noise = True
                break
        if is_line_noise: continue
            
        line_after_removal = INLINE_NOISE_PATTERN.sub('', stripped_line)
        line_after_removal = re.sub(r'[\s|·/]+', '', line_after_removal)
        if len(line_after_removal) < 5 and INLINE_NOISE_PATTERN.search(stripped_line):
            is_line_noise = True
        
        if is_line_noise: continue
        
        cleaned_lines.append(line)

    content_after_pass1 = "\n".join(cleaned_lines)
    # 阶段二：净化行内链接
    content_after_pass2 = INLINE_NOISE_PATTERN.sub('', content_after_pass1)
    # 阶段三：格式化
    final_text = re.sub(r'\n{3,}', '\n\n', content_after_pass2)
    return final_text.strip()

def diagnose_only(content):
    """
    【诊断函数】只分析不清理，返回一个包含所有待办操作的日志列表。
    """
    issues = []
    lines = content.split('\n')
    
    # 阶段一：分析整行删除
    lines_to_be_deleted = set()
    for i, line in enumerate(lines):
        line_num = i + 1
        stripped_line = line.strip()
        
        for pattern in WHOLE_LINE_NOISE_PATTERNS:
            if pattern.search(stripped_line):
                issues.append(f"  [L{line_num}][删除整行] 发现版权信息: '{stripped_line}'")
                lines_to_be_deleted.add(i)
                break
        if i in lines_to_be_deleted: continue

        line_after_removal = INLINE_NOISE_PATTERN.sub('', stripped_line)
        line_after_removal = re.sub(r'[\s|·/]+', '', line_after_removal)
        if len(line_after_removal) < 5 and INLINE_NOISE_PATTERN.search(stripped_line):
            issues.append(f"  [L{line_num}][删除整行] 可能是导航菜单: '{stripped_line}'")
            lines_to_be_deleted.add(i)

    # 阶段二：分析行内净化
    for i, line in enumerate(lines):
        if i in lines_to_be_deleted: continue
        line_num = i + 1
        matches = INLINE_NOISE_PATTERN.findall(line)
        for match in matches:
            issues.append(f"  [L{line_num}][净化行内] 将移除链接: '{match}'")
            
    return sorted(list(set(issues)))


# ================= 4. 主程序 =================
def main():
    # --- 用户选择模式 ---
    print("="*50)
    print("欢迎使用数据清洗工具包 (保守版)")
    print("="*50)
    print("请选择要执行的操作：")
    print("  1. 生成诊断报告 (只读模式，检查将要发生什么)")
    print("  2. 执行清理操作 (实际写入新文件到目标目录)")
    
    choice = ""
    while choice not in ['1', '2']:
        choice = input("请输入您的选择 (1 或 2): ")

    if not os.path.isdir(SOURCE_DIR):
        print(f"❌ 错误：找不到输入目录: {SOURCE_DIR}")
        return

    # --- 根据选择执行不同任务 ---
    if choice == '1':
        # --- 诊断模式 ---
        print("\n🚀 开始生成诊断报告 (只读模式)...")
        print(f"📂 扫描源目录: {SOURCE_DIR}")
        print(f"📜 报告将保存至: {REPORT_FILE}")
        
        total_files_with_issues = 0
        with open(REPORT_FILE, 'w', encoding='utf-8') as report_file:
            report_file.write(f"=== 数据清洗诊断报告 (保守版) ===\n扫描目录: {SOURCE_DIR}\n\n")

            for root, _, files in os.walk(SOURCE_DIR):
                for filename in files:
                    if filename.lower().endswith(".md"):
                        filepath = os.path.join(root, filename)
                        try:
                            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                                content = f.read()
                            
                            found_issues = diagnose_only(content)
                            if found_issues:
                                total_files_with_issues += 1
                                print(f"📄 发现潜在问题: {filename}")
                                report_file.write(f"📄 文件: {filepath}\n")
                                for issue in found_issues:
                                    report_file.write(f"{issue}\n")
                                report_file.write("-" * 20 + "\n")
                        except Exception as e:
                            print(f"⚠️ 读取文件时发生错误: {filename} - {e}")
                            
        print("-" * 50)
        print(f"✅ 诊断完成！共在 {total_files_with_issues} 个文件中发现待办操作。")
        print(f"🎉 请查看报告文件 '{os.path.basename(REPORT_FILE)}' 进行审核。")

    elif choice == '2':
        # --- 清理模式 ---
        print("\n🚀 即将执行清理操作...")
        print(f"📂 读取源目录: {SOURCE_DIR}")
        print(f"📂 写入目标目录: {DEST_DIR}")
        print("🚨 警告：此操作将创建新的已清理文件。")
        
        confirm = input("您确定要继续吗？(y/n): ")
        if confirm.lower() != 'y':
            print("操作已取消。")
            return

        print("正在执行清理...")
        file_count = 0
        os.makedirs(DEST_DIR, exist_ok=True)
        
        for root, _, files in os.walk(SOURCE_DIR):
            for filename in files:
                if filename.lower().endswith(".md"):
                    source_path = os.path.join(root, filename)
                    relative_path = os.path.relpath(source_path, SOURCE_DIR)
                    dest_path = os.path.join(DEST_DIR, relative_path)
                    os.makedirs(os.path.dirname(dest_path), exist_ok=True)
                    
                    try:
                        with open(source_path, 'r', encoding='utf-8', errors='ignore') as f:
                            raw_content = f.read()
                        
                        cleaned_content = clean_markdown_conservatively(raw_content)
                        
                        # 只有内容发生改变，且不为空时才写入
                        if cleaned_content and raw_content != cleaned_content:
                            with open(dest_path, 'w', encoding='utf-8') as f:
                                f.write(cleaned_content)
                            file_count += 1
                            print(f"✨ 已清理并保存: {filename}")
                    except Exception as e:
                        print(f"⚠️ 处理文件时发生错误: {filename} - {e}")
        
        print("-" * 50)
        print(f"✅ 清洗完成！共对 {file_count} 个文件进行了清理并生成了新文件。")
        print(f"🎉 最终数据已保存在: {DEST_DIR}")

if __name__ == "__main__":
    main()
