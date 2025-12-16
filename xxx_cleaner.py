import os
import sys
import datetime

# --- 配置 ---
LOG_FILE = "cleaning_log.txt"  # 日志文件名
TRIGGER_STRING = "![](data:image/png;base64"

def clean_base64_blocks(root_dir, log_file):
    """
    清理.md文件中嵌入的 base64 图片块。
    删除从包含 TRIGGER_STRING 的行开始，直到下一个空行为止的所有行。

    :param root_dir: 要处理的根目录路径。
    :param log_file: 用于写入日志的文件对象。
    """
    if not os.path.isdir(root_dir):
        print(f"错误：提供的路径 '{root_dir}' 不是一个有效的文件夹。")
        return

    print(f"开始扫描并清理 Base64 图片块：{root_dir}")
    print(f"详细日志将被记录在: {os.path.abspath(LOG_FILE)}")
    
    processed_files_count = 0
    cleaned_files_count = 0

    for dirpath, _, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.lower().endswith('.md'):
                file_path = os.path.join(dirpath, filename)
                processed_files_count += 1
                
                try:
                    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                        original_content = f.read()

                    # 如果文件中不包含触发词，直接跳过，提高效率
                    if TRIGGER_STRING not in original_content:
                        print(f"[{processed_files_count}] 正在检查: {file_path} -> [无需清理]")
                        continue
                    
                    print(f"[{processed_files_count}] 正在检查: {file_path}", end="")

                    original_lines = original_content.split('\n')
                    cleaned_lines = []
                    all_deleted_blocks_info = [] # 用于日志记录
                    
                    is_deleting = False
                    
                    for line in original_lines:
                        # 检查是否触发删除模式
                        if not is_deleting and TRIGGER_STRING in line:
                            is_deleting = True
                            # 记录被删除块的开头，用于日志
                            all_deleted_blocks_info.append(line[:80]) # 只记录前80个字符作为摘要
                            continue  # 跳过当前行，不将其加入cleaned_lines
                        
                        # 如果在删除模式中，检查是否遇到空行以停止删除
                        # 空行的定义是去除首尾空白后为空字符串
                        if is_deleting and line.strip() == "":
                            is_deleting = False
                            # 这个空行本身我们是需要保留的，所以它会和正常行一样被添加
                        
                        # 如果不在删除模式，就保留该行
                        if not is_deleting:
                            cleaned_lines.append(line)

                    # 将清理后的行重新组合成文件内容
                    cleaned_content = "\n".join(cleaned_lines)

                    # 只有内容发生变化时才执行写操作和日志记录
                    if cleaned_content != original_content:
                        with open(file_path, 'w', encoding='utf-8') as f:
                            f.write(cleaned_content)
                        
                        # 写入日志
                        timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
                        log_message = f"{timestamp} - CLEARED: {file_path}\n"
                        for i, block_start in enumerate(all_deleted_blocks_info, 1):
                            log_message += f"  └─ Removed Block {i} starting with: \"{block_start}...\"\n"
                        log_file.write(log_message + "\n")
                        
                        cleaned_files_count += 1
                        print(" -> [图片块已清理 - 详情已写入日志]")
                    else:
                        # 这种情况很少见，比如图片块在文件末尾且没有内容变化
                        print(" -> [无需清理]")

                except Exception as e:
                    print(f" -> [处理失败] 错误: {e}")
                    log_file.write(f"ERROR processing {file_path}: {e}\n")

    print("\n--------------------")
    print("      处理完成      ")
    print("--------------------")
    print(f"总共检查了 {processed_files_count} 个 .md 文件。")
    print(f"其中 {cleaned_files_count} 个文件进行了清理操作。")
    print(f"请查看日志文件 '{LOG_FILE}' 获取详细的清理记录。")


if __name__ == "__main__":
    print("=" * 50)
    print("!!! 重要：此脚本将直接修改您文件夹中的 .md 文件。")
    print("!!! 强烈建议在运行前，先备份您的整个文件夹！")
    print("=" * 50)

    target_directory = input("请输入要处理的 Markdown 文件所在的文件夹路径：\n> ")
    target_directory = target_directory.strip('\"\'')
    
    confirm = input(f"您确定要清理 '{target_directory}' 文件夹中的 Base64 图片块吗？(y/n): ")

    if confirm.lower() == 'y':
        try:
            with open(LOG_FILE, 'a', encoding='utf-8') as log_f:
                log_f.write("\n" + "="*20 + " NEW BASE64 CLEANING SESSION " + "="*20 + "\n")
                log_f.write(f"Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                log_f.write(f"Target Directory: {os.path.abspath(target_directory)}\n")
                log_f.write("-" * 69 + "\n")
                
                clean_base64_blocks(target_directory, log_f)

        except Exception as e:
            print(f"\n发生严重错误，操作中断: {e}")
            print(f"请检查您对 '{LOG_FILE}' 文件是否有写入权限。")
    else:
        print("操作已取消。")
    
    if sys.platform.startswith('win'):
        input("按 Enter 键退出...")

