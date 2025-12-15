import asyncio
from playwright.async_api import async_playwright

URL = "https://jlsteenwyk.com/PhyKIT/tutorials/index.html"
SELECTOR = "#bd-docs-nav a.nav-link"
SCREENSHOT_PATH = "debug_screenshot.png"

async def main():
    print("==============================================")
    print("🚀 启动PhyKIT网站的独立调试脚本...")
    print("==============================================")
    
    page = None
    browser = None
    
    try:
        async with async_playwright() as p:
            print("1. 正在启动一个可见的浏览器 (headless=False)...")
            browser = await p.chromium.launch(headless=False, slow_mo=500)
            context = await browser.new_context(
                user_agent="Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
            )
            page = await context.new_page()

            print(f"2. 正在导航到: {URL} (超时时间设置为60秒)")
            await page.goto(URL, wait_until="networkidle", timeout=60000)

            print("3. 页面已加载。现在等待导航栏出现...")
            print(f"   (选择器: '{SELECTOR}')")
            await page.wait_for_selector(SELECTOR, timeout=30000)

            print("\n✅✅✅ 成功！选择器已找到！✅✅✅")
            links_count = await page.locator(SELECTOR).count()
            print(f"   -> 找到了 {links_count} 个导航链接。")
            print("   -> 这意味着网络和选择器都是正确的。问题可能出在主脚本的集成上。")

    except Exception as e:
        print("\n❌❌❌ 调试过程中发生错误！❌❌❌")
        print(f"错误类型: {type(e).__name__}")
        print(f"错误信息: {e}")
        
        if page:
            print("\n📸 正在对当前浏览器窗口进行截图...")
            await page.screenshot(path=SCREENSHOT_PATH, full_page=True)
            print(f"   -> 截图已保存为: {SCREENSHOT_PATH}")
            print("   -> 请立即打开这个图片文件！它会显示浏览器最后看到的内容。")
            print("   -> 检查图片中是否有验证码、Cookie弹窗、广告或错误信息。")

    finally:
        if browser and browser.is_connected():
            print("\n4. 按下回车键关闭浏览器...")
            input() # 等待用户确认，以便有足够时间观察浏览器
            await browser.close()
            print("   -> 浏览器已关闭。")

if __name__ == "__main__":
    # 在Windows上运行需要此策略
    if asyncio.get_event_loop().is_running() and os.name == 'nt':
         asyncio.set_event_loop_policy(asyncio.WindowsSelectorEventLoopPolicy())
    asyncio.run(main())
