# 文件名: main.py (V7.0 - 智能比对模式)
import os
import sys
import yaml
import re
import time
from datetime import datetime
from bs4 import BeautifulSoup
from markdownify import markdownify as md

# 导入 Selenium 相关库
try:
    from selenium import webdriver
    from selenium.webdriver.chrome.service import Service
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    from webdriver_manager.chrome import ChromeDriverManager
    from urllib.parse import urljoin
except ImportError as e:
    print(f"!!! 关键错误: 缺少必要的库: {e.name}。")
    print(f"!!! 请在你的虚拟环境 (rag_crawler) 中运行: pip install {e.name}")
    sys.exit(1)


def sanitize_filename(text):
    text = re.sub(r'https?://[^/]+', '', text)
    text = text.strip('/').replace('/', '_').replace('.html', '')
    text = re.sub(r'^\d+\s*', '', text) 
    filename = re.sub(r'[\\/*?:"<>|]', "", text).strip()
    return f"{filename[:150]}.md" if filename else "index.md"

def get_unique_filepath(directory, base_filename):
    """如果文件已存在，则通过添加 "(n)" 后缀来创建唯一的路径。"""
    filepath = os.path.join(directory, base_filename)
    if not os.path.exists(filepath):
        return filepath, base_filename

    base, ext = os.path.splitext(base_filename)
    counter = 1
    while True:
        new_filename = f"{base} ({counter}){ext}"
        new_filepath = os.path.join(directory, new_filename)
        if not os.path.exists(new_filepath):
            return new_filepath, new_filename
        counter += 1

def is_content_identical(filepath, new_md_content):
    """
    智能比对文件内容，忽略文件头。
    返回 True 如果核心内容相同，否则返回 False。
    """
    if not os.path.exists(filepath):
        return False
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            old_full_content = f.read()
        
        # 通过 '---' 分割线分离文件头和核心内容
        parts = old_full_content.split('\n---\n', 1)
        
        # 如果文件中有分割线，取其后的内容；否则，认为整个文件都是核心内容（兼容旧文件）
        old_md_content = parts[1].strip() if len(parts) > 1 else old_full_content.strip()
        
        return old_md_content == new_md_content.strip()
    except Exception:
        # 如果读取或解析旧文件失败，则认为内容不同，需要重新保存
        return False

def main():
    if len(sys.argv) < 2:
        print(f"用法: python {os.path.basename(__file__)} <config_file.yaml>")
        sys.exit(1)
        
    config_file = sys.argv[1]
    
    with open(config_file, 'r', encoding='utf-8') as f:
        all_configs = yaml.safe_load(f)['sites']

    print("诊断: 准备 ChromeDriver...")
    try:
        driver_path = ChromeDriverManager().install()
        print("诊断: ChromeDriver 已就绪。")
    except Exception as e:
        print(f"!!! 致命错误: ChromeDriver 安装失败: {e}")
        sys.exit(1)

    for config in all_configs:
        site_name = config['name']
        output_dir = config['output']['directory']
        parser = config['parser']
        remove_selectors = parser.get('remove_selectors', [])

        print(f"\n=======================================================")
        print(f"[{site_name}] 开始爬取...")
        print(f"=======================================================")
        
        page_urls = []
        
        if 'start_url' in config and 'discovery_selector' in config:
            # ... (这部分链接发现逻辑与之前完全相同)
            print("步骤 1: 进入链接发现模式...")
            start_url = config['start_url']
            discovery_selector = config['discovery_selector']
            temp_driver = None
            try:
                print(f"  -> 正在启动临时浏览器以发现链接...")
                options = webdriver.ChromeOptions(); options.add_argument('--headless')
                service = Service(driver_path)
                temp_driver = webdriver.Chrome(service=service, options=options)
                temp_driver.get(start_url)
                print(f"  -> 等待目录链接出现 (Selector: {discovery_selector})...")
                WebDriverWait(temp_driver, 30).until(EC.presence_of_element_located((By.CSS_SELECTOR, discovery_selector)))
                link_elements = temp_driver.find_elements(By.CSS_SELECTOR, discovery_selector)
                discovered_links = set([urljoin(start_url, el.get_attribute('href')) for el in link_elements if el.get_attribute('href')])
                page_urls = sorted(list(discovered_links))
                print(f"  -> 成功发现 {len(page_urls)} 个不重复的页面链接。")
            except Exception as e:
                print(f"!!! 错误: 在链接发现阶段失败: {e}")
                continue
            finally:
                if temp_driver: temp_driver.quit()
        
        if not page_urls:
            print("警告: 未能获取任何页面链接，跳过此网站。")
            continue

        os.makedirs(output_dir, exist_ok=True)
        saved_count = 0
        skipped_count = 0
        total_pages = len(page_urls)
        
        driver = None
        try:
            print("\n步骤 2: 启动主浏览器实例，开始逐一爬取页面...")
            options = webdriver.ChromeOptions(); options.add_argument('--headless'); options.add_argument('--log-level=3'); options.add_argument('--window-size=1920,1080'); options.add_experimental_option('excludeSwitches', ['enable-logging']); options.add_argument('--disable-dev-shm-usage'); options.add_argument('--no-sandbox')
            service = Service(driver_path)
            driver = webdriver.Chrome(service=service, options=options)
            driver.set_page_load_timeout(300)
            wait = WebDriverWait(driver, 60)
            
            for i, url in enumerate(page_urls, 1):
                print(f"\n--- 正在处理页面 ({i}/{total_pages}): {url} ---")
                
                try:
                    driver.get(url)
                    wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, parser['main_container_selector'])))
                    soup = BeautifulSoup(driver.page_source, 'html.parser')
                    main_container = soup.select_one(parser['main_container_selector'])

                    if not main_container: continue

                    if remove_selectors:
                        for selector in remove_selectors:
                            for element in main_container.select(selector): element.decompose()

                    # 优先基于标题生成文件名
                    page_title_tag = main_container.find(['h1', 'h2'])
                    base_filename = sanitize_filename(page_title_tag.get_text(strip=True)) if page_title_tag else sanitize_filename(url.split('/')[-1])
                    
                    # 生成核心Markdown内容，用于比对
                    content_md = md(str(main_container), heading_style="ATX").strip()
                    
                    # ==========================================================
                    # 【核心智能比对逻辑】
                    # ==========================================================
                    output_path = os.path.join(output_dir, base_filename)
                    if os.path.exists(output_path) and is_content_identical(output_path, content_md):
                        print(f"  -> 内容无变化，跳过: {base_filename}")
                        skipped_count += 1
                        continue
                    
                    # 只有当文件不存在或内容已更新时，才继续执行保存逻辑
                    date_scraped = datetime.now().strftime("%Y-%m-%d")
                    file_header = f"Source URL: {url}\nDate Scraped: {date_scraped}\n\n---\n"
                    final_content = f"{file_header}\n{content_md}"

                    # 使用版本控制逻辑获取最终保存路径
                    unique_filepath, final_filename = get_unique_filepath(output_dir, base_filename)

                    with open(unique_filepath, 'w', encoding='utf-8') as f: 
                        f.write(final_content)
                    
                    print(f"  -> 成功保存到: {final_filename}")
                    saved_count += 1

                except Exception as e:
                    print(f"  -> 处理页面 {url} 时发生错误: {e}")
                finally:
                    time.sleep(1)
        finally:
            if driver: driver.quit()

        print(f"\n[{site_name}] 所有任务完成！")
        print(f"[{site_name}] 共成功保存 {saved_count} 个文件，内容无变化跳过 {skipped_count} 个文件。")

if __name__ == "__main__":
    main()

