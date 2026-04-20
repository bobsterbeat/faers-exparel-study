"""
Build a self-contained HTML report that combines:
  - The full manuscript skeleton (from manuscript.md, rendered to HTML)
  - All three figures (base64-embedded PNGs)
  - All three tables (rendered as HTML tables, styled for print)

Produces a single file `report.html` you can:
  - Open in any browser to read
  - Print to PDF from the browser (File -> Print -> Save as PDF)
  - Email / share as a single attachment (no external dependencies)

Usage: python make_report.py
"""

from __future__ import annotations

import base64
import re
from pathlib import Path

import pandas as pd


def _img_to_b64(path: Path) -> str:
    data = path.read_bytes()
    b64 = base64.b64encode(data).decode()
    return f"data:image/png;base64,{b64}"


def _md_to_html(md: str) -> str:
    """Minimal, deterministic Markdown-to-HTML converter — good enough for
    a manuscript skeleton with headings, paragraphs, lists, emphasis, and
    tables. No external deps."""
    out = []
    lines = md.split("\n")
    i = 0
    in_list = False
    in_table = False
    table_buf: list[str] = []

    def _flush_list():
        nonlocal in_list
        if in_list:
            out.append("</ul>")
            in_list = False

    def _flush_table():
        nonlocal in_table, table_buf
        if not in_table:
            return
        rows = [r for r in table_buf if r.strip()]
        if len(rows) < 2:
            in_table = False
            table_buf = []
            return
        header = [c.strip() for c in rows[0].strip("|").split("|")]
        body_rows = []
        for r in rows[2:]:  # skip header + separator
            cells = [c.strip() for c in r.strip("|").split("|")]
            body_rows.append(cells)
        out.append("<table><thead><tr>")
        for h in header:
            out.append(f"<th>{_inline(h)}</th>")
        out.append("</tr></thead><tbody>")
        for cells in body_rows:
            out.append("<tr>")
            for c in cells:
                out.append(f"<td>{_inline(c)}</td>")
            out.append("</tr>")
        out.append("</tbody></table>")
        in_table = False
        table_buf = []

    def _inline(s: str) -> str:
        # Escape HTML special chars except for our own tags
        s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        # **bold**
        s = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", s)
        # *italic* or _italic_
        s = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<em>\1</em>", s)
        # [text](url)
        s = re.sub(r"\[(.+?)\]\((.+?)\)", r'<a href="\2">\1</a>', s)
        # `code`
        s = re.sub(r"`(.+?)`", r"<code>\1</code>", s)
        return s

    while i < len(lines):
        line = lines[i]

        # Table block
        if line.strip().startswith("|") and line.strip().endswith("|"):
            if not in_table:
                _flush_list()
                in_table = True
                table_buf = []
            table_buf.append(line)
            i += 1
            continue
        else:
            _flush_table()

        # Headings
        m = re.match(r"^(#{1,6})\s+(.+)$", line)
        if m:
            _flush_list()
            lvl = len(m.group(1))
            out.append(f"<h{lvl}>{_inline(m.group(2))}</h{lvl}>")
            i += 1
            continue

        # Horizontal rule
        if re.match(r"^---+$", line.strip()):
            _flush_list()
            out.append("<hr>")
            i += 1
            continue

        # Unordered list
        if line.lstrip().startswith(("- ", "* ")):
            if not in_list:
                out.append("<ul>")
                in_list = True
            content = line.lstrip()[2:]
            out.append(f"<li>{_inline(content)}</li>")
            i += 1
            continue
        else:
            _flush_list()

        # Blank line
        if not line.strip():
            i += 1
            continue

        # Paragraph (collect until blank/structural line)
        para = [line]
        j = i + 1
        while j < len(lines) and lines[j].strip() and not lines[j].lstrip().startswith((
            "#", "- ", "* ", "|", "---")):
            para.append(lines[j])
            j += 1
        out.append(f"<p>{_inline(' '.join(para))}</p>")
        i = j

    _flush_list()
    _flush_table()
    return "\n".join(out)


def _table_to_html(df: pd.DataFrame, caption: str) -> str:
    # Escape HTML chars in cell values.
    def esc(s):
        s = str(s) if pd.notna(s) else ""
        return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

    html = [f'<figure class="tablefig"><figcaption>{esc(caption)}</figcaption>',
            '<table class="data">']
    html.append("<thead><tr>")
    for c in df.columns:
        html.append(f"<th>{esc(c)}</th>")
    html.append("</tr></thead>")
    html.append("<tbody>")
    for _, row in df.iterrows():
        html.append("<tr>")
        for val in row:
            html.append(f"<td>{esc(val)}</td>")
        html.append("</tr>")
    html.append("</tbody></table></figure>")
    return "\n".join(html)


CSS = """
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  max-width: 920px;
  margin: 2em auto;
  padding: 0 1.5em;
  color: #222;
  line-height: 1.55;
  font-size: 15px;
}
h1 { font-size: 1.7em; border-bottom: 2px solid #2c3e50; padding-bottom: 0.3em; margin-top: 1em; }
h2 { font-size: 1.35em; color: #2c3e50; margin-top: 1.6em; border-bottom: 1px solid #ddd; padding-bottom: 0.2em; }
h3 { font-size: 1.1em; color: #34495e; margin-top: 1.3em; }
h4 { font-size: 1em; color: #34495e; }
code { background: #f4f4f4; padding: 0.1em 0.4em; border-radius: 3px; font-size: 0.9em; }
a { color: #2980b9; }
hr { border: none; border-top: 1px solid #ccc; margin: 1.6em 0; }

table { border-collapse: collapse; margin: 1em 0; width: 100%; font-size: 0.92em; }
th { background: #2c3e50; color: white; text-align: left; padding: 0.5em 0.7em; }
td { padding: 0.45em 0.7em; border-bottom: 1px solid #eee; }
tr:nth-child(even) td { background: #f8f8f8; }

.tablefig figcaption { font-weight: bold; margin-bottom: 0.3em; color: #2c3e50; }
figure { margin: 1.5em 0; }
figure.figure-img { text-align: center; }
figure.figure-img img { max-width: 100%; height: auto; border: 1px solid #ddd; }
figure.figure-img figcaption { margin-top: 0.5em; font-size: 0.92em; color: #555; font-style: italic; }

@media print {
  body { max-width: none; margin: 0; padding: 0.5cm; font-size: 10pt; }
  h1 { page-break-before: avoid; }
  h2 { page-break-after: avoid; }
  figure, table { page-break-inside: avoid; }
}
"""


def main() -> None:
    root = Path(__file__).parent
    md = (root / "manuscript.md").read_text()
    manuscript_html = _md_to_html(md)

    # Tables
    t1 = pd.read_csv(root / "out/tables/table_1_cohort.csv")
    t2 = pd.read_csv(root / "out/tables/table_2_pt_signals.csv")
    t3 = pd.read_csv(root / "out/tables/table_3_composite_era.csv")

    fig1_b64 = _img_to_b64(root / "fig/fig1_composite_forest.png")
    fig2_b64 = _img_to_b64(root / "fig/fig2_pt_forest.png")
    fig3_b64 = _img_to_b64(root / "fig/fig3_era_trend.png")
    fig4_b64 = _img_to_b64(root / "fig/fig4_pathology_pharmacology.png")

    figs_html = f"""
<h2>Figures</h2>
<figure class="figure-img">
  <img src="{fig1_b64}" alt="Figure 1">
  <figcaption>
    <strong>Figure 1.</strong> Era-stratified disproportionality: case-level composite endpoints
    for Exparel vs plain bupivacaine (left) and ropivacaine (right). Top row: prolonged
    sensory/motor block composite. Bottom row: LAST-spectrum composite. Red = signal meets
    all four criteria. a/b counts shown at right of each error bar.
  </figcaption>
</figure>
<figure class="figure-img">
  <img src="{fig2_b64}" alt="Figure 2">
  <figcaption>
    <strong>Figure 2.</strong> PT-level forest plot: individual MedDRA preferred terms within
    the prolonged-block and peripheral-neuropathy spectrum, Exparel vs plain bupivacaine.
    Red = four-criterion signal. a/b counts at right margin.
  </figcaption>
</figure>
<figure class="figure-img">
  <img src="{fig3_b64}" alt="Figure 3">
  <figcaption>
    <strong>Figure 3.</strong> Exparel FAERS reporting volume (top) and era-evolution of the
    two composite endpoints (bottom), 2012-2025. Dashed lines mark the 2018 ISB indication
    expansion, the February 2021 Anesthesiology meta-analysis, and the 2025 NOPAIN Act.
    The LAST rate peaks 2013-2014 then collapses (Weber effect); the prolonged-block rate
    rises steadily post-2018.
  </figcaption>
</figure>
<figure class="figure-img">
  <img src="{fig4_b64}" alt="Figure 4">
  <figcaption>
    <strong>Figure 4.</strong> Duration and permanence stratification of the prolonged-block
    composite. Pathology-implying PTs (MedDRA terms that cannot represent expected &lt;72 h
    pharmacology) signal independently vs plain bupivacaine (ROR 2.16), ruling out pure
    drug-action explanation. Expected-pharmacology PTs signal more strongly (ROR 3.10) but
    cannot by themselves distinguish normal from prolonged. Hypoaesthesia oral after
    non-oral Exparel is examined as a systemic-absorption marker. Red = signal meets all
    four criteria. None of the strata signal vs ropivacaine, arguing against a generic
    long-acting-amide effect.
  </figcaption>
</figure>
"""

    tables_html = "\n".join([
        "<h2>Tables</h2>",
        _table_to_html(t1, "Table 1. Cohort demographics and outcome distribution"),
        _table_to_html(t2, "Table 2. PT-level signals (prolonged-block and peripheral-neuropathy spectrum)"),
        _table_to_html(t3, "Table 3. Composite-endpoint disproportionality by reporting era"),
    ])

    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Exparel FAERS Disproportionality — Draft Manuscript</title>
<style>{CSS}</style>
</head>
<body>
{manuscript_html}
{figs_html}
{tables_html}
</body>
</html>"""

    out = root / "report.html"
    out.write_text(full_html)
    size_kb = out.stat().st_size / 1024
    print(f"wrote {out} ({size_kb:.0f} KB, self-contained)")


if __name__ == "__main__":
    main()
