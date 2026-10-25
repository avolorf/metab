// Сборка пояснительной записки (DOCX): node docs/build_report.js  (из каталога r_utility)
const fs = require("fs");
const path = require("path");
const D = require("docx");
const { Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType, ShadingType, AlignmentType,
        HeadingLevel, ImageRun, PageBreak, LevelFormat, Footer, PageNumber, TableOfContents, BorderStyle } = D;

const ROOT = path.resolve(__dirname, "..");
const rd = (p) => fs.readFileSync(path.join(ROOT, p), "utf8");
const csv = (p) => {
  const lines = rd(p).trim().split(/\r?\n/);
  const split = (l) => { const out = []; let cur = "", q = false;
    for (const ch of l) { if (ch === '"') q = !q; else if (ch === "," && !q) { out.push(cur); cur = ""; } else cur += ch; }
    out.push(cur); return out; };
  const h = split(lines[0]); return lines.slice(1).map((l) => { const v = split(l); const o = {}; h.forEach((k, i) => (o[k] = v[i])); return o; });
};
const content = require("./report_content.js")({ rd, csv, ROOT });

const FONT = "Times New Roman", MONO = "Consolas";
const PAGE_W = 11906, MARG = 1134, TW = PAGE_W - 1701 - 850;   // A4, поля 3 / 1.5 / 2 / 2 см

function runs(text, opts = {}) {
  // **жирный**, `код`
  const parts = String(text).split(/(\*\*[^*]+\*\*|`[^`]+`)/);
  return parts.filter((s) => s.length).map((s) => {
    if (s.startsWith("**")) return new TextRun({ text: s.slice(2, -2), bold: true, font: FONT, ...opts });
    if (s.startsWith("`")) return new TextRun({ text: s.slice(1, -1), font: MONO, size: (opts.size || 24) - 2 });
    return new TextRun({ text: s, font: FONT, ...opts });
  });
}
const P = (t, o = {}) => new Paragraph({ children: runs(t, o.run || {}), spacing: { after: 120, line: 300 },
  alignment: o.align || AlignmentType.JUSTIFIED, indent: o.noIndent ? undefined : { firstLine: 567 }, ...o.para });
const cell = (t, w, head = false) => new TableCell({
  width: { size: w, type: WidthType.DXA }, margins: { top: 40, bottom: 40, left: 80, right: 80 },
  shading: head ? { fill: "E7E6E6", type: ShadingType.CLEAR, color: "auto" } : undefined,
  children: String(t).split("\n").map((line) => new Paragraph({ children: runs(line, { size: 20, bold: head }), spacing: { after: 0 } })) });
function table(head, rows, widths) {
  const tot = widths.reduce((a, b) => a + b, 0); const ws = widths.map((w) => Math.round((w / tot) * TW));
  ws[ws.length - 1] += TW - ws.reduce((a, b) => a + b, 0);
  return new Table({ width: { size: TW, type: WidthType.DXA }, columnWidths: ws,
    rows: [new TableRow({ tableHeader: true, children: head.map((h, i) => cell(h, ws[i], true)) }),
           ...rows.map((r) => new TableRow({ children: r.map((c, i) => cell(c, ws[i])) }))] });
}
const code = (text) => text.replace(/\t/g, "  ").split("\n").map((l, i, a) => new Paragraph({
  children: [new TextRun({ text: l.length ? l : " ", font: MONO, size: 15 })], spacing: { after: 0, line: 220 },
  shading: { fill: "F4F4F4", type: ShadingType.CLEAR, color: "auto" }, keepLines: true }));

function image(p, wcm, caption) {
  const buf = fs.readFileSync(path.join(ROOT, p));
  const w = buf.readUInt32BE(16), h = buf.readUInt32BE(20);
  const wpx = Math.round(wcm / 2.54 * 96), hpx = Math.round(wpx * h / w);
  return [new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 120, after: 60 },
           children: [new ImageRun({ type: "png", data: buf, transformation: { width: wpx, height: hpx },
                                     altText: { title: caption, description: caption, name: path.basename(p) } })] }),
          new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 }, children: runs(caption, { size: 22 }) })];
}

const children = [];
for (const b of content) {
  const [k, ...a] = b;
  if (k === "title") children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 200 }, children: runs(a[0], { bold: true, size: 30 }) }));
  else if (k === "center") children.push(new Paragraph({ alignment: AlignmentType.CENTER, spacing: { after: 120 }, children: runs(a[0]) }));
  else if (k === "h1") children.push(new Paragraph({ heading: HeadingLevel.HEADING_1, pageBreakBefore: !a[1], children: [new TextRun({ text: a[0], font: FONT })] }));
  else if (k === "h2") children.push(new Paragraph({ heading: HeadingLevel.HEADING_2, children: [new TextRun({ text: a[0], font: FONT })] }));
  else if (k === "h3") children.push(new Paragraph({ heading: HeadingLevel.HEADING_3, children: [new TextRun({ text: a[0], font: FONT })] }));
  else if (k === "p") children.push(P(a[0]));
  else if (k === "caption") children.push(new Paragraph({ spacing: { before: 160, after: 80 }, keepNext: true, children: runs(a[0], { size: 22 }) }));
  else if (k === "bullets") a[0].forEach((t) => children.push(new Paragraph({ numbering: { reference: "bul", level: 0 }, spacing: { after: 60, line: 280 }, alignment: AlignmentType.JUSTIFIED, children: runs(t) })));
  else if (k === "numbers") a[0].forEach((t) => children.push(new Paragraph({ numbering: { reference: a[1] || "num", level: 0 }, spacing: { after: 60, line: 280 }, alignment: AlignmentType.JUSTIFIED, children: runs(t) })));
  else if (k === "table") { children.push(table(a[0], a[1], a[2])); children.push(new Paragraph({ spacing: { after: 120 }, children: [] })); }
  else if (k === "code") { children.push(...code(a[0])); children.push(new Paragraph({ spacing: { after: 120 }, children: [] })); }
  else if (k === "img") children.push(...image(a[0], a[1], a[2]));
  else if (k === "toc") children.push(new TableOfContents("Содержание", { hyperlink: true, headingStyleRange: "1-2" }));
  else if (k === "pagebreak") children.push(new Paragraph({ children: [new PageBreak()] }));
}

const numCfg = (ref, fmt, text) => ({ reference: ref, levels: [{ level: 0, format: fmt, text, alignment: AlignmentType.LEFT,
  style: { paragraph: { indent: { left: 720, hanging: 360 } }, run: { font: FONT } } }] });
const doc = new Document({
  creator: "UroMetabR", title: "UroMetabR: код и пояснения (2 этап НИР)",
  styles: { default: { document: { run: { font: FONT, size: 24 } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 28, bold: true, font: FONT }, paragraph: { spacing: { before: 240, after: 200 }, outlineLevel: 0, keepNext: true } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 26, bold: true, font: FONT }, paragraph: { spacing: { before: 200, after: 120 }, outlineLevel: 1, keepNext: true } },
      { id: "Heading3", name: "Heading 3", basedOn: "Normal", next: "Normal", quickFormat: true,
        run: { size: 24, bold: true, italics: true, font: FONT }, paragraph: { spacing: { before: 160, after: 100 }, outlineLevel: 2, keepNext: true } }] },
  numbering: { config: [numCfg("bul", LevelFormat.BULLET, "–"), numCfg("num", LevelFormat.DECIMAL, "%1."),
                        numCfg("num2", LevelFormat.DECIMAL, "%1."), numCfg("num3", LevelFormat.DECIMAL, "%1."), numCfg("num4", LevelFormat.DECIMAL, "%1.")] },
  features: { updateFields: true },
  sections: [{ properties: { page: { size: { width: PAGE_W, height: 16838 }, margin: { top: MARG, bottom: MARG, left: 1701, right: 850 } } },
    footers: { default: new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] }) },
    children }] });
const out = path.join(ROOT, "docs", "UroMetabR_2_etap_kod_i_poyasneniya.docx");
Packer.toBuffer(doc).then((b) => { fs.writeFileSync(out, b); console.log("→", out); });
