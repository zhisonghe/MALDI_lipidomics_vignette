"""Render Tutorial.md as a self-contained web page.

Usage: python make_html.py [Tutorial.md] [Tutorial.html] [--fragment page_fragment.html]
Needs markdown-it-py (3.0) and mdit-py-plugins (0.4). Figures from figures/ are embedded as data URIs, so the page
is a single file. The default output is a standalone HTML document; --fragment additionally writes the same page
without the <html>/<head>/<body> skeleton (for hosts that add their own).
"""
import base64, html, os, re, sys
from markdown_it import MarkdownIt
from mdit_py_plugins.anchors import anchors_plugin

args = [a for a in sys.argv[1:] if not a.startswith("--")]
here = os.path.dirname(os.path.abspath(__file__))
src = args[0] if args else os.path.join(here, "..", "Tutorial.md")
out = args[1] if len(args) > 1 else os.path.join(os.path.dirname(src), "Tutorial.html")
frag = sys.argv[sys.argv.index("--fragment") + 1] if "--fragment" in sys.argv else None
base = os.path.dirname(os.path.abspath(src))
md_text = open(src, encoding="utf-8").read()


def gh_slug(s):
    """GitHub-style heading anchor, so the tutorial's own internal links keep working."""
    s = re.sub(r"<[^>]+>", "", s).strip().lower()
    return re.sub(r"[^\w\- ]", "", s).replace(" ", "-")


# the Markdown table of contents becomes the sidebar; the "(Back to top)" links are replaced by the sidebar too
title = re.search(r"^# (.+)$", md_text, re.M).group(1)
subtitle = re.search(r"^#### (.+)$", md_text, re.M).group(1)
body = re.sub(r"^# .+\n#### .+\n", "", md_text, count=1, flags=re.M)
body = re.sub(r"### Table of Contents\n.*?(?=\n## Introduction)", "", body, count=1, flags=re.S)
body = body.replace('<sub><a href="#top">(Back to top)</a></sub></br>', "")

md = MarkdownIt("commonmark", {"html": True, "typographer": False}).enable("table")
md.use(anchors_plugin, min_level=2, max_level=4, slug_func=gh_slug, permalink=False)
html_body = md.render(body)

# figures: embed images and turn "image + italic caption" into <figure>
def embed(m):
    path = os.path.join(base, m.group(1))
    data = base64.b64encode(open(path, "rb").read()).decode()
    return f'src="data:image/png;base64,{data}"'
html_body = re.sub(r'src="(figures/[^"]+\.png)"', embed, html_body)
html_body = re.sub(r"<p><img ([^>]*?)\s*/?>\s*<em>(Figure .*?)</em></p>",
                   r'<figure><img \1 loading="lazy"><figcaption>\2</figcaption></figure>', html_body, flags=re.S)
html_body = re.sub(r"<table>", '<div class="table-wrap"><table>', html_body).replace("</table>", "</table></div>")
html_body = re.sub(r"<blockquote>\s*<p><strong>NOTE</strong></p>", '<aside class="note"><p class="note-label">Note</p>', html_body)
html_body = re.sub(r"<blockquote>\s*<p><strong>NOTE</strong>\s*", '<aside class="note"><p class="note-label">Note</p>\n<p>', html_body)
html_body = re.sub(r"(<aside class=\"note\">.*?)</blockquote>", r"\1</aside>", html_body, flags=re.S)
html_body = html_body.replace('class="language-console"', 'class="language-bash"')

# sidebar contents from the rendered h2/h3 headings
toc = []
for level, hid, text in re.findall(r'<h([23]) id="([^"]+)">(.*?)</h\1>', html_body):
    toc.append(f'<li class="l{level}"><a href="#{hid}">{text}</a></li>')
toc_html = "\n".join(toc)
n_sections = sum(1 for t in toc if 'class="l3"' in t)

PAGE = """<title>MALDI-MSI Lipidomics Tutorial</title>
<meta name="claude:viewport-fit-sides" content="cover">
<meta name="description" content="Preprocessing and analysis of MALDI mass spectrometry imaging (spatial lipidomics) data.">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;600;700&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&family=JetBrains+Mono:wght@400;600&display=swap">
<style>
/* Layout: sticky contents rail on the left, one reading column (~70ch) on the right; wide tables and code scroll in place. */
:root {
  --bg: #fbfafc; --surface: #f3f0f6; --fg: #1f1b24; --muted: #5f5869; --rule: #e1dbe7;
  --accent: #7a2e8e; --accent-soft: #f1e6f4; --code-bg: #f4f1f7;
  --hl-keyword: #7a2e8e; --hl-string: #2d6e45; --hl-number: #a3461a; --hl-title: #22559e; --hl-comment: #77707f;
  --font-head: "Archivo", "Helvetica Neue", Arial, sans-serif;
  --font-body: "Source Serif 4", "Iowan Old Style", Georgia, serif;
  --font-mono: "JetBrains Mono", "SFMono-Regular", Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #141218; --surface: #1c1922; --fg: #ebe7f0; --muted: #aaa3b4; --rule: #2f2a37;
  --accent: #d39be2; --accent-soft: #2a2032; --code-bg: #1a1720;
  --hl-keyword: #d39be2; --hl-string: #8fd0a5; --hl-number: #f0a77a; --hl-title: #93b8f2; --hl-comment: #8d8697;
  color-scheme: dark; } }
:root[data-theme="dark"] {
  --bg: #141218; --surface: #1c1922; --fg: #ebe7f0; --muted: #aaa3b4; --rule: #2f2a37;
  --accent: #d39be2; --accent-soft: #2a2032; --code-bg: #1a1720;
  --hl-keyword: #d39be2; --hl-string: #8fd0a5; --hl-number: #f0a77a; --hl-title: #93b8f2; --hl-comment: #8d8697;
  color-scheme: dark; }
html { scroll-padding-top: 16px; }
body { background: var(--bg); color: var(--fg); font-family: var(--font-body); font-size: 17px; line-height: 1.62;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px)); }
.shell { display: grid; grid-template-columns: 17rem minmax(0, 1fr); gap: 3rem; max-width: 76rem; margin: 0 auto; }
nav.toc { position: sticky; top: env(safe-area-inset-top, 0px); align-self: start; max-height: 100vh; overflow-y: auto;
  padding-block: 2rem 2rem; font-family: var(--font-head); font-size: 0.82rem; line-height: 1.35; }
nav.toc .toc-title { font-size: 0.72rem; letter-spacing: 0.08em; text-transform: uppercase; color: var(--muted); margin: 0 0 0.6rem; }
nav.toc ul { list-style: none; margin: 0; padding: 0; display: grid; gap: 0.15rem; }
nav.toc a { display: block; padding: 0.22rem 0.5rem; border-radius: 4px; color: var(--muted); text-decoration: none; }
nav.toc li.l2 a { color: var(--fg); font-weight: 600; margin-top: 0.55rem; }
nav.toc li.l3 a { padding-left: 1.1rem; }
nav.toc a:hover { color: var(--accent); }
nav.toc a.active { background: var(--accent-soft); color: var(--accent); }
details.toc-mobile { display: none; }
main { min-width: 0; padding-block: 2.5rem 5rem; }
.article { max-width: 46rem; }
header.masthead { padding-bottom: 1.5rem; margin-bottom: 2rem; border-bottom: 1px solid var(--rule); }
.eyebrow { font-family: var(--font-head); font-size: 0.75rem; letter-spacing: 0.1em; text-transform: uppercase; color: var(--accent); margin: 0 0 0.6rem; }
h1 { font-family: var(--font-head); font-weight: 700; font-size: clamp(1.7rem, 1.2rem + 2vw, 2.45rem); line-height: 1.15; margin: 0 0 0.75rem; text-wrap: balance; }
.meta { font-family: var(--font-head); font-size: 0.85rem; color: var(--muted); margin: 0; }
h2, h3, h4 { font-family: var(--font-head); line-height: 1.25; text-wrap: balance; }
h2 { font-size: 1.6rem; font-weight: 700; margin: 3.2rem 0 1rem; padding-top: 1.2rem; border-top: 1px solid var(--rule); }
h3 { font-size: 1.22rem; font-weight: 600; margin: 2.4rem 0 0.8rem; }
h4 { font-size: 1.02rem; font-weight: 600; margin: 1.8rem 0 0.6rem; }
p, li { max-width: 70ch; }
a { color: var(--accent); text-underline-offset: 2px; }
a:focus-visible, button:focus-visible, summary:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px; }
code { font-family: var(--font-mono); font-size: 0.84em; background: var(--code-bg); padding: 0.08em 0.3em; border-radius: 3px; overflow-wrap: anywhere; }
pre { position: relative; background: var(--code-bg); border: 1px solid var(--rule); border-radius: 6px; padding: 0.9rem 1rem;
  overflow-x: auto; font-size: 0.8rem; line-height: 1.5; margin: 1rem 0 1.3rem; }
pre code { background: none; padding: 0; font-size: inherit; overflow-wrap: normal; white-space: pre; }
button.copy { position: absolute; top: 0.45rem; right: 0.45rem; font-family: var(--font-head); font-size: 0.72rem;
  color: var(--muted); background: var(--bg); border: 1px solid var(--rule); border-radius: 4px; padding: 0.15rem 0.5rem; cursor: pointer; opacity: 0.85; }
button.copy:hover { color: var(--accent); border-color: var(--accent); }
.table-wrap { overflow-x: auto; margin: 1rem 0 1.4rem; border: 1px solid var(--rule); border-radius: 6px; }
table { border-collapse: collapse; font-family: var(--font-head); font-size: 0.8rem; line-height: 1.4; min-width: 100%; }
th, td { text-align: left; vertical-align: top; padding: 0.5rem 0.65rem; border-bottom: 1px solid var(--rule); }
th { background: var(--surface); font-weight: 600; white-space: nowrap; }
tr:last-child td { border-bottom: none; }
td code { font-size: 0.92em; }
aside.note { background: var(--accent-soft); border-left: 3px solid var(--accent); border-radius: 0 6px 6px 0; padding: 0.4rem 1.1rem 0.6rem; margin: 1.2rem 0 1.5rem; font-size: 0.95rem; }
aside.note .note-label { font-family: var(--font-head); font-size: 0.72rem; letter-spacing: 0.08em; text-transform: uppercase; color: var(--accent); margin: 0.5rem 0 0.2rem; }
aside.note ul { padding-left: 1.2rem; }
aside.note pre { background: var(--bg); }
blockquote { margin: 1rem 0; padding-left: 1rem; border-left: 3px solid var(--rule); color: var(--muted); }
figure { margin: 1.6rem 0 2rem; max-width: 100%; }
figure img { display: block; max-width: 100%; height: auto; border: 1px solid var(--rule); border-radius: 6px; background: #ffffff; }
figcaption { font-size: 0.86rem; color: var(--muted); margin-top: 0.6rem; max-width: 70ch; }
hr { border: none; border-top: 1px solid var(--rule); }
.hljs-keyword, .hljs-literal, .hljs-built_in, .hljs-meta { color: var(--hl-keyword); }
.hljs-string, .hljs-regexp { color: var(--hl-string); }
.hljs-number { color: var(--hl-number); }
.hljs-title, .hljs-function .hljs-title, .hljs-params, .hljs-attr { color: var(--hl-title); }
.hljs-comment { color: var(--hl-comment); font-style: italic; }
@media (max-width: 960px) {
  .shell { grid-template-columns: minmax(0, 1fr); gap: 0; }
  nav.toc { display: none; }
  details.toc-mobile { display: block; border: 1px solid var(--rule); border-radius: 6px; margin: 1rem 0 0; font-family: var(--font-head); font-size: 0.9rem; }
  details.toc-mobile summary { cursor: pointer; padding: 0.6rem 0.9rem; font-weight: 600; }
  details.toc-mobile ul { list-style: none; margin: 0; padding: 0 0.9rem 0.8rem; display: grid; gap: 0.2rem; }
  details.toc-mobile li.l3 { padding-left: 1rem; }
  details.toc-mobile a { color: var(--fg); text-decoration: none; }
  main { padding-block: 1.2rem 4rem; }
  body { font-size: 16px; }
}
@media (prefers-reduced-motion: no-preference) { html { scroll-behavior: smooth; } }
</style>

<div class="shell" id="top">
  <nav class="toc" aria-label="Contents">
    <p class="toc-title">Contents</p>
    <ul>
__TOC__
    </ul>
  </nav>
  <main>
    <div class="article">
      <header class="masthead">
        <p class="eyebrow">MALDI-MSI &middot; spatial lipidomics &middot; m/z, pixels and lipids</p>
        <h1>__TITLE__</h1>
        <p class="meta">__SUBTITLE__ &middot; __NSEC__ sections &middot; companion files: <code>envs/</code>, <code>scripts/</code>, <code>figures/</code></p>
        <details class="toc-mobile"><summary>Contents</summary><ul>
__TOC__
        </ul></details>
      </header>
__BODY__
    </div>
  </main>
</div>

<script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>
<script>
(function () {
  if (window.hljs) { hljs.configure({ ignoreUnescapedHTML: true }); document.querySelectorAll("pre code").forEach(function (el) { hljs.highlightElement(el); }); }
  document.querySelectorAll("pre").forEach(function (pre) {
    var b = document.createElement("button"); b.className = "copy"; b.type = "button"; b.textContent = "Copy";
    b.addEventListener("click", function () {
      var text = pre.querySelector("code") ? pre.querySelector("code").innerText : pre.innerText;
      var done = function () { b.textContent = "Copied"; setTimeout(function () { b.textContent = "Copy"; }, 1500); };
      var select = function () { var r = document.createRange(); r.selectNodeContents(pre.querySelector("code") || pre);
        var s = window.getSelection(); s.removeAllRanges(); s.addRange(r); b.textContent = "Selected: press Ctrl/Cmd+C"; };
      try { navigator.clipboard.writeText(text).then(done, select); } catch (e) { select(); }
    });
    pre.appendChild(b);
  });
  var links = {}; document.querySelectorAll("nav.toc a").forEach(function (a) { links[a.getAttribute("href").slice(1)] = a; });
  var heads = Array.prototype.slice.call(document.querySelectorAll("main h2[id], main h3[id]"));
  if ("IntersectionObserver" in window) {
    var current = null;
    var obs = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) { if (e.isIntersecting) { current = e.target.id; } });
      Object.keys(links).forEach(function (k) { links[k].classList.toggle("active", k === current); });
    }, { rootMargin: "0px 0px -75% 0px" });
    heads.forEach(function (h) { obs.observe(h); });
  }
  document.querySelectorAll("details.toc-mobile a").forEach(function (a) { a.addEventListener("click", function () { a.closest("details").open = false; }); });
})();
</script>
"""
page = (PAGE.replace("__TOC__", toc_html).replace("__TITLE__", html.escape(title))
        .replace("__SUBTITLE__", html.escape(subtitle)).replace("__NSEC__", str(n_sections)).replace("__BODY__", html_body))
if frag:
    open(frag, "w", encoding="utf-8").write(page)
head, rest = page.split("<style>", 1)
style, content = rest.split("</style>", 1)
standalone = ("<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
              "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
              + head + "<style>\n:root { color-scheme: light; } body { margin: 0; } img { max-width: 100%; }\n" + style
              + "</style>\n</head>\n<body>\n" + content + "\n</body>\n</html>\n")
open(out, "w", encoding="utf-8").write(standalone)
print(f"{out}: {len(standalone) / 1e6:.2f} MB, {len(toc)} contents entries, "
      f"{html_body.count('<figure>')} figures, {html_body.count('<pre>')} code blocks, {html_body.count('<table>')} tables")
