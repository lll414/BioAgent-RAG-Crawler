# 文件名: main.py (终极健壮版 V3 - 增加页面加载超时)
import os
import sys
import yaml
import re
import time
from datetime import datetime
from bs4 import BeautifulSoup
from markdownify import markdownify as md

# 导入 Selenium 相关库
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager
from urllib.parse import urljoin

def sanitize_filename(url_or_title):
    """使用 URL 或标题生成唯一且安全的文件名。"""
    text = re.sub(r'https?://[^/]+', '', url_or_title)
    text = text.strip('/').replace('/', '_').replace('.html', '')
    filename = re.sub(r'[\\*?:"<>|]', "", text)
    return f"{filename[:150]}.md" if filename else "index.md"

def main():
    if len(sys.argv) < 2:
        print(f"用法: python {os.path.basename(__file__)} <config_file.yaml>")
        sys.exit(1)
        
    config_file = sys.argv[1]
    
    with open(config_file, 'r', encoding='utf-8') as f:
        all_configs = yaml.safe_load(f)['sites']

    print("正在准备 ChromeDriver...")
    driver_path = ChromeDriverManager().install()
    print("ChromeDriver 已就绪。")

    for config in all_configs:
        site_name = config['name']
        output_dir = config['output']['directory']
        parser = config['parser']

        print(f"\n=======================================================")
        print(f"[{site_name}] 开始爬取 (模板化 URL 模式 - 终极健壮版 V3)...")
        print(f"=======================================================")
        
        page_urls = []
        if 'url_template' in config:
            template = config['url_template']
            start = config.get('template_start', 0)
            end = config.get('template_end', 14)
            page_urls = [template.format(i=f"{i:02d}") for i in range(start, end + 1)]
            print(f"步骤 1: 根据模板生成了 {len(page_urls)} 个页面链接。")
        else:
            print("错误: 配置文件中未找到 'url_template'。")
            continue

        os.makedirs(output_dir, exist_ok=True)
        saved_count = 0
        total_pages = len(page_urls)
        
        print("\n步骤 2: 开始逐一爬取、清理和保存每个页面...")
        for i, url in enumerate(page_urls, 1):
            print(f"\n--- 正在处理页面 ({i}/{total_pages}): {url} ---")
            
            driver = None
            try:
                options = webdriver.ChromeOptions()
                options.add_argument('--headless')
                options.add_argument('--log-level=3')
                options.add_argument('--window-size=1920,1080')
                options.add_experimental_option('excludeSwitches', ['enable-logging'])
                options.add_argument('--disable-dev-shm-usage') 
                options.add_argument('--no-sandbox')
                
                service = Service(driver_path)
                driver = webdriver.Chrome(service=service, options=options)
                
                # --- 【核心解决方案】 ---
                # 设置一个非常长的页面加载超时时间（300秒 = 5分钟）
                # 这直接解决了因页面内容过多、加载缓慢导致的通信超时问题。
                print("  -> 设置5分钟的页面加载超时...")
                driver.set_page_load_timeout(500)
                
                # 显式等待的时间也可以相应放长一些，以增加保险
                wait = WebDriverWait(driver, 60)

                print("  -> 正在加载页面...")
                driver.get(url)
                print("  -> 页面加载指令已发送，等待页面内容出现...")

                wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, parser['main_container_selector'])))
                print("  -> 页面主内容已出现，开始解析。")
                
                soup = BeautifulSoup(driver.page_source, 'html.parser')
                main_container = soup.select_one(parser['main_container_selector'])

                if not main_container:
                    print(f"  -> 警告: 在页面 {url} 未能找到主内容容器，已跳过。")
                    continue

                # HTML层面和Markdown层面的清理逻辑保持不变
                exercise_header = main_container.find('h2', id='exercises')
                if exercise_header:
                    for sibling in exercise_header.find_all_next(string=False): sibling.decompose()
                    exercise_header.decompose()
                
                content_md = md(str(main_container), heading_style="ATX").strip()
                lines = content_md.splitlines()
                cut_off_index = -1
                for idx, line in enumerate(lines):
                    if line.strip().startswith("##") and "exercises" in line.lower():
                        cut_off_index = idx
                        break
                if cut_off_index != -1:
                    content_md = "\n".join(lines[:cut_off_index]).strip()
                
                # 文件保存逻辑保持不变
                filename = sanitize_filename(url)
                output_path = os.path.join(output_dir, filename)
                date_scraped = datetime.now().strftime("%Y-%m-%d")
                file_header = f"Source URL: {url}\nDate Scraped: {date_scraped}\n\n---\n"
                final_content = f"{file_header}\n{content_md}"
                with open(output_path, 'w', encoding='utf-8') as f: f.write(final_content)
                
                print(f"  -> 成功保存到: {filename}")
                saved_count += 1
                
            except Exception as e:
                print(f"  -> 处理页面 {url} 时发生严重错误: {e}")
            finally:
                if driver:
                    driver.quit()
                print("  -> 暂停3秒，以确保系统资源释放...")
                time.sleep(3)

        print(f"\n[{site_name}] 所有任务完成！")
        print(f"[{site_name}] 共成功保存 {saved_count}/{total_pages} 个文件到目录: {output_dir}")

if __name__ == "__main__":
    main()
