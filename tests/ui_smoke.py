"""Browser smoke test for a running mmcifXplorer (python app.py).

Usage:  python tests/ui_smoke.py load,conn,atom,crystal,comp,dict,tour,id=1D3Z [screenshot_dir]
Needs:  pip install playwright && playwright install chromium
Prints Dash callback errors and browser console errors for each step."""
import sys
import time

from playwright.sync_api import sync_playwright

URL = "http://127.0.0.1:8050"
OUT = (sys.argv[2].rstrip("/") + "/") if len(sys.argv) > 2 else "./"
steps = sys.argv[1].split(",") if len(sys.argv) > 1 else ["load"]


def wait_idle(pg, ms=1500):
    pg.wait_for_timeout(300)
    for _ in range(60):
        busy = pg.evaluate("() => document.title.includes('Updating') || !!document.querySelector('._dash-loading-callback')")
        if not busy:
            break
        pg.wait_for_timeout(250)
    pg.wait_for_timeout(ms)


def errors_in_page(pg):
    return pg.eval_on_selector_all(".dash-fe-error__title, .dash-fe-error-top", "e=>e.map(x=>x.textContent)")


with sync_playwright() as p:
    b = p.chromium.launch(args=["--use-gl=swiftshader", "--enable-unsafe-swiftshader"])
    pg = b.new_page(viewport={"width": 1600, "height": 1000})
    logs = []
    pg.on("console", lambda m: logs.append(f"{m.type}: {m.text[:300]}") if m.type == "error" else None)
    pg.on("pageerror", lambda e: logs.append(f"PAGEERR {e}"))
    pg.goto(URL)
    wait_idle(pg, 6000)
    print("status:", pg.inner_text("#status"))
    for st in steps:
        if st == "load":
            pg.screenshot(path=OUT + "ui_load.png")
        elif st == "conn":
            pg.click("#cat-search")
            pg.fill("#cat-search", "struct_conn")
            pg.keyboard.press("Enter")
            wait_idle(pg)
            pg.click("button:has-text('_struct_conn')")
            wait_idle(pg)
            pg.click("#cat-table td[data-dash-column='pdbx_dist_value'] >> nth=0")
            wait_idle(pg, 2500)
            print("banner:", pg.inner_text("#hl-banner"))
            pg.screenshot(path=OUT + "ui_conn.png")
        elif st == "atom":
            box = pg.locator("#viewer-wrap canvas").bounding_box()
            hit = False
            for dx, dy in [(0, 0), (10, 5), (-15, 10), (20, -12), (-30, -20), (40, 30), (5, -40)]:
                pg.mouse.click(box["x"] + box["width"] / 2 + dx, box["y"] + box["height"] / 2 + dy)
                pg.wait_for_timeout(900)
                txt = pg.inner_text("#atom-panel")
                if "_atom_site row" in txt:
                    hit = True
                    break
            print("atom picked:", hit, txt[:200].replace("\n", " | "))
            pg.screenshot(path=OUT + "ui_atom.png")
        elif st == "crystal":
            pg.click("a:has-text('Crystal & quality')")
            wait_idle(pg, 2500)
            pg.screenshot(path=OUT + "ui_crystal.png", full_page=True)
        elif st == "comp":
            pg.click("a:has-text('Composition')")
            wait_idle(pg, 2500)
            pg.screenshot(path=OUT + "ui_comp.png", full_page=True)
        elif st == "dict":
            pg.click("a:has-text('Dictionary')")
            pg.fill("#dict-search", "R_free")
            pg.keyboard.press("Enter")
            wait_idle(pg)
            pg.click(".dict-result >> nth=1")
            wait_idle(pg)
            pg.screenshot(path=OUT + "ui_dict.png")
        elif st == "tour":
            pg.click("a:has-text('Learn & tours')")
            wait_idle(pg)
            pg.click("button:has-text('Anatomy of an X-ray entry')")
            wait_idle(pg, 2000)
            for _ in range(2):
                pg.click("#tour-next")
                wait_idle(pg, 1500)
            print("tour:", pg.inner_text("#tour-bar")[:200].replace("\n", " "))
            print("banner:", pg.inner_text("#hl-banner")[:200])
            pg.screenshot(path=OUT + "ui_tour.png")
        elif st.startswith("id="):
            pg.fill("#pdb-input", st[3:])
            pg.click("#btn-load")
            wait_idle(pg, 5000)
            print("status:", pg.inner_text("#status"))
            pg.screenshot(path=OUT + f"ui_{st[3:]}.png")
        print(f"[{st}] dash errors:", errors_in_page(pg))
    print("console errors:")
    for l in logs:
        if "3dmol.csb.pitt.edu" in l or "dash-version.plotly.com" in l or "ERR_FAILED" in l:
            continue
        print("  ", l)
    b.close()
