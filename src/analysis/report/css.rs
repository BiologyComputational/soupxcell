// ============================================================================
// report/css.rs  —  §0  Professional full-width report CSS
// ============================================================================

pub fn report_css() -> String {
    r#"
/* ─── Reset ──────────────────────────────────────────────────────────────── */
*,*::before,*::after{box-sizing:border-box;margin:0;padding:0}
html{scroll-behavior:smooth;font-size:15px}
body{background:#F0F4F8;color:#0F172A;
     font-family:'DM Sans',system-ui,sans-serif;
     line-height:1.7;min-width:320px}

/* ─── Typography ─────────────────────────────────────────────────────────── */
h1{font-family:'Source Serif 4',serif;font-size:2.5rem;font-weight:700;
   line-height:1.2;color:#fff;margin-bottom:12px}
h2{font-family:'Source Serif 4',serif;font-size:1.5rem;font-weight:600;
   color:#0F172A;margin:0 0 8px}
h3{font-size:1.05rem;font-weight:600;color:#1E293B;margin:0 0 10px}
h4{font-size:.95rem;font-weight:600;color:#334155;margin:0 0 6px}
p{color:#334155;margin-bottom:12px;line-height:1.75}
p:last-child{margin-bottom:0}
code{font-family:'JetBrains Mono',monospace;font-size:.82em;
     background:#E8EFF7;color:#1E40AF;border-radius:4px;padding:1px 6px}
pre{font-family:'JetBrains Mono',monospace;font-size:.82em;
    background:#0F172A;color:#E2E8F0;border-radius:8px;
    padding:16px 20px;overflow-x:auto;line-height:1.6;margin:14px 0}
strong{color:#0F172A;font-weight:600}
a{color:#2563EB;text-decoration:none}
a:hover{text-decoration:underline}
ol,ul{padding-left:20px}
li{margin-bottom:6px;color:#334155;font-size:.9rem}

/* ─── Page wrapper ───────────────────────────────────────────────────────── */
.cover{background:linear-gradient(140deg,#0F172A 0%,#1E3A5F 55%,#0C4A6E 100%);
       color:#F1F5F9;width:100%}
.cover-inner{max-width:1240px;margin:0 auto;padding:72px 56px 60px}
.page{max-width:1240px;margin:0 auto;padding:0 56px}
.page-section{padding:40px 0 0}
@media(max-width:900px){
  .cover-inner,.page{padding-left:24px;padding-right:24px}
  .two-col,.three-col{grid-template-columns:1fr!important}
  .finding-grid{grid-template-columns:1fr!important}
  .interp-grid{grid-template-columns:1fr!important}
}

/* ─── Grid helpers ───────────────────────────────────────────────────────── */
.two-col{display:grid;grid-template-columns:1fr 1fr;gap:28px;align-items:start}
.three-col{display:grid;grid-template-columns:1fr 1fr 1fr;gap:20px;align-items:start}

/* ─── Section typography ─────────────────────────────────────────────────── */
.section-label{font-size:.72rem;font-weight:700;letter-spacing:.14em;
               text-transform:uppercase;color:#64748B;margin-bottom:4px}
.section-desc{font-size:.88rem;color:#64748B;margin-bottom:22px;
              line-height:1.7}
.section-divider{border:none;border-top:2px solid #E2E8F0;margin:36px 0}

/* ─── Cover ──────────────────────────────────────────────────────────────── */
.cover-eyebrow{font-size:.74rem;letter-spacing:.16em;text-transform:uppercase;
               color:#94A3B8;margin-bottom:16px}
.cover-meta{display:flex;flex-wrap:wrap;gap:8px;margin-top:24px}
.cover-badge{background:rgba(255,255,255,.1);border:1px solid rgba(255,255,255,.2);
             border-radius:20px;padding:5px 16px;font-size:.78rem;color:#CBD5E1;
             white-space:nowrap}
.cover-badge.pass{background:rgba(34,197,94,.18);border-color:rgba(34,197,94,.35);
                  color:#86EFAC}
.cover-badge.info{background:rgba(6,182,212,.15);border-color:rgba(6,182,212,.3);
                  color:#67E8F9}

/* ─── TOC ────────────────────────────────────────────────────────────────── */
.toc{background:#FFF;border:1px solid #CBD5E1;border-radius:12px;
     padding:24px 28px;margin:36px 0 0}
.toc h3{font-size:.78rem;text-transform:uppercase;letter-spacing:.12em;
        color:#64748B;margin-bottom:14px;font-family:'DM Sans',sans-serif}
.toc ol{list-style:decimal;padding-left:22px;
        columns:2;column-gap:40px}
@media(max-width:640px){.toc ol{columns:1}}
.toc li{margin-bottom:6px;font-size:.88rem;break-inside:avoid}
.toc a{color:#2563EB;font-weight:500}
.toc a:hover{color:#1D4ED8}

/* ─── Metric cards ───────────────────────────────────────────────────────── */
.metric-grid{display:grid;
             grid-template-columns:repeat(auto-fill,minmax(210px,1fr));
             gap:14px;margin-bottom:0}
.metric-card{background:#FFF;border:1px solid #E2E8F0;border-radius:12px;
             padding:20px 22px;position:relative;overflow:hidden;
             box-shadow:0 1px 4px rgba(0,0,0,.06)}
.metric-card::before{content:'';position:absolute;top:0;left:0;right:0;height:4px}
.metric-card.good::before{background:linear-gradient(90deg,#22C55E,#16A34A)}
.metric-card.info::before{background:linear-gradient(90deg,#06B6D4,#0891B2)}
.metric-card.warn::before{background:linear-gradient(90deg,#F59E0B,#D97706)}
.metric-card.neutral::before{background:#CBD5E1}
.mc-value{font-family:'JetBrains Mono',monospace;font-size:1.6rem;font-weight:700;
          color:#0F172A;line-height:1.1;margin-bottom:6px}
.mc-label{font-size:.74rem;color:#64748B;text-transform:uppercase;
          letter-spacing:.09em}
.mc-delta{font-size:.78rem;margin-top:7px;font-family:'JetBrains Mono',monospace}
.mc-delta.pos{color:#15803D}.mc-delta.neg{color:#DC2626}.mc-delta.zero{color:#64748B}

/* ─── Card wrapper ───────────────────────────────────────────────────────── */
.card{background:#FFF;border:1px solid #E2E8F0;border-radius:12px;
      padding:24px 28px;box-shadow:0 1px 4px rgba(0,0,0,.05)}

/* ─── Story / narrative blocks ───────────────────────────────────────────── */
.story-block{border-left:4px solid #2563EB;padding:16px 20px;
             background:#EFF6FF;border-radius:0 10px 10px 0;margin:18px 0}
.story-block.green{border-left-color:#16A34A;background:#F0FDF4}
.story-block.amber{border-left-color:#D97706;background:#FFFBEB}
.story-block.red  {border-left-color:#DC2626;background:#FEF2F2}
.story-block.purple{border-left-color:#7C3AED;background:#F5F3FF}
.story-block.teal {border-left-color:#0891B2;background:#ECFEFF}
.story-label{font-size:.72rem;font-weight:700;letter-spacing:.1em;
             text-transform:uppercase;margin-bottom:8px;display:block}
.story-label.blue  {color:#1D4ED8}
.story-label.green {color:#15803D}
.story-label.amber {color:#B45309}
.story-label.red   {color:#B91C1C}
.story-label.purple{color:#6D28D9}
.story-label.teal  {color:#0E7490}

/* ─── Finding boxes ──────────────────────────────────────────────────────── */
.finding-grid{display:grid;grid-template-columns:1fr 1fr 1fr;
              gap:16px;margin:20px 0}
@media(max-width:1000px){.finding-grid{grid-template-columns:1fr 1fr}}
@media(max-width:640px) {.finding-grid{grid-template-columns:1fr}}
.finding-box{border-radius:12px;padding:20px 22px;
             box-shadow:0 1px 3px rgba(0,0,0,.05)}
.finding-box.blue  {background:#EFF6FF;border:1px solid #BFDBFE}
.finding-box.green {background:#F0FDF4;border:1px solid #BBF7D0}
.finding-box.amber {background:#FFFBEB;border:1px solid #FDE68A}
.finding-box.purple{background:#F5F3FF;border:1px solid #DDD6FE}
.finding-box.teal  {background:#ECFEFF;border:1px solid #A5F3FC}
.finding-box.red   {background:#FEF2F2;border:1px solid #FECACA}
.finding-icon {font-size:1.6rem;margin-bottom:10px;display:block}
.finding-title{font-size:.82rem;font-weight:700;margin-bottom:6px;
               color:#0F172A;text-transform:uppercase;letter-spacing:.06em}
.finding-value{font-family:'JetBrains Mono',monospace;font-size:1.25rem;
               font-weight:700;margin-bottom:8px}
.finding-desc {font-size:.84rem;color:#334155;line-height:1.65}
.finding-plain{font-size:.78rem;color:#475569;margin-top:10px;padding-top:10px;
               border-top:1px solid rgba(0,0,0,.08);line-height:1.6;
               font-style:italic}

/* ─── Hypothesis / objective boxes ──────────────────────────────────────── */
.hypo-box{background:#F8F5FF;border:1px solid #DDD6FE;border-radius:10px;
          padding:18px 22px;margin-bottom:14px}
.hypo-box.question {background:#EFF6FF;border-color:#BFDBFE}
.hypo-box.objective{background:#F0FDF4;border-color:#BBF7D0}
.hypo-box.finding  {background:#FFFBEB;border-color:#FDE68A}
.hypo-label{font-size:.7rem;font-weight:700;letter-spacing:.12em;
            text-transform:uppercase;margin-bottom:8px;display:block}
.hypo-label.q{color:#1D4ED8}
.hypo-label.o{color:#15803D}
.hypo-label.f{color:#B45309}

/* ─── Timeline ───────────────────────────────────────────────────────────── */
.timeline{position:relative;padding-left:36px;margin-top:16px}
.timeline::before{content:'';position:absolute;left:11px;top:12px;
                  bottom:12px;width:2px;
                  background:linear-gradient(to bottom,#2563EB,#7C3AED,#16A34A)}
.tl-step{position:relative;margin-bottom:22px;padding-top:2px}
.tl-step::before{content:attr(data-n);position:absolute;left:-36px;top:0;
                 width:22px;height:22px;background:#2563EB;color:#FFF;
                 border-radius:50%;font-size:.72rem;font-weight:700;
                 display:flex;align-items:center;justify-content:center}
.tl-step:nth-child(2)::before{background:#7C3AED}
.tl-step:nth-child(3)::before{background:#0891B2}
.tl-step:nth-child(4)::before{background:#16A34A}
.tl-step:nth-child(5)::before{background:#D97706}
.tl-step h4{font-size:.9rem;font-weight:600;color:#0F172A;margin-bottom:4px}
.tl-step p{font-size:.84rem;margin-bottom:0;color:#334155}

/* ─── Tables ─────────────────────────────────────────────────────────────── */
.table-wrap{overflow-x:auto;margin:14px 0;
            border:1px solid #E2E8F0;border-radius:10px;
            box-shadow:0 1px 3px rgba(0,0,0,.04)}
table{width:100%;border-collapse:collapse;font-size:.86rem}
thead tr{background:#F8FAFC}
th{padding:11px 16px;text-align:left;font-weight:600;color:#475569;
   font-size:.76rem;text-transform:uppercase;letter-spacing:.07em;
   border-bottom:2px solid #E2E8F0}
td{padding:10px 16px;border-bottom:1px solid #F1F5F9;
   color:#334155;vertical-align:top}
tbody tr:last-child td{border-bottom:none}
tbody tr:hover td{background:#FAFBFD}
td.num,th.num{text-align:right;font-family:'JetBrains Mono',monospace}
td.mono{font-family:'JetBrains Mono',monospace;font-size:.82em}
td.good{color:#15803D;font-weight:600}
td.warn{color:#B45309;font-weight:600}
td.info{color:#0E7490;font-weight:600}

/* ─── Formula block ──────────────────────────────────────────────────────── */
.formula-block{background:#0F172A;border-radius:12px;
               padding:24px 28px;color:#E2E8F0}
.formula-line{display:grid;grid-template-columns:28px 1fr 1fr;
              gap:12px;align-items:start;margin-bottom:16px}
.formula-line:last-child{margin-bottom:0}
.formula-step{background:#2563EB;color:#FFF;width:22px;height:22px;
              border-radius:50%;font-size:.72rem;font-weight:700;
              display:flex;align-items:center;justify-content:center;
              flex-shrink:0;margin-top:2px}
.formula-expr{font-family:'JetBrains Mono',monospace;font-size:.9rem;
              color:#E2E8F0;line-height:1.5}
.formula-note{font-size:.8rem;color:#94A3B8;line-height:1.55;padding-top:3px}

/* ─── Interpretation cards ───────────────────────────────────────────────── */
.interp-grid{display:grid;grid-template-columns:1fr 1fr;gap:16px}
.interp-card{background:#FFF;border:1px solid #E2E8F0;border-radius:12px;
             padding:20px 24px;box-shadow:0 1px 3px rgba(0,0,0,.04)}
.ic-title{font-size:.76rem;font-weight:700;text-transform:uppercase;
          letter-spacing:.1em;margin-bottom:10px;padding-bottom:8px;
          border-bottom:2px solid #F1F5F9}
.ic-title.teal {color:#0E7490;border-bottom-color:#CFFAFE}
.ic-title.green{color:#15803D;border-bottom-color:#DCFCE7}
.ic-title.amber{color:#B45309;border-bottom-color:#FEF3C7}
.ic-title.red  {color:#B91C1C;border-bottom-color:#FEE2E2}
.ic-title.blue {color:#1D4ED8;border-bottom-color:#DBEAFE}
.ic-verdict{margin-top:12px;padding:10px 14px;border-radius:8px;
            font-size:.84rem;line-height:1.65}
.ic-verdict.good{background:#F0FDF4;color:#15803D;
                 border-left:3px solid #22C55E}
.ic-verdict.info{background:#ECFEFF;color:#0E7490;
                 border-left:3px solid #06B6D4}
.ic-verdict.warn{background:#FFFBEB;color:#92400E;
                 border-left:3px solid #F59E0B}

/* ─── Simulation animated bars ───────────────────────────────────────────── */
.sim-bar-row{display:grid;grid-template-columns:48px 1fr 100px;
             align-items:center;gap:10px;margin-bottom:10px}
.sim-bar-label{font-size:.82rem;color:#475569;
               font-family:'JetBrains Mono',monospace;text-align:right}
.sim-bar-track{height:20px;background:#F1F5F9;border-radius:10px;
               overflow:hidden;border:1px solid #E2E8F0}
.sim-bar-fill{height:100%;border-radius:10px;width:0%;
              transition:width .9s cubic-bezier(.4,0,.2,1);
              display:flex;align-items:center;justify-content:flex-end;
              padding-right:8px}
.sim-bar-fill span{font-size:.72rem;font-weight:700;color:#FFF;
                   font-family:'JetBrains Mono',monospace;white-space:nowrap}
.sim-bar-sil{font-size:.78rem;color:#64748B;
             font-family:'JetBrains Mono',monospace;text-align:right}

/* ─── Badge ──────────────────────────────────────────────────────────────── */
.badge{display:inline-block;border-radius:20px;padding:3px 12px;
       font-size:.74rem;font-weight:600}
.badge.pass  {background:#DCFCE7;color:#15803D}
.badge.info  {background:#CFFAFE;color:#0E7490}
.badge.warn  {background:#FEF3C7;color:#92400E}
.badge.fail  {background:#FEE2E2;color:#B91C1C}

/* ─── Plot iframes ───────────────────────────────────────────────────────── */
.plots-grid{display:grid;
            grid-template-columns:repeat(auto-fill,minmax(540px,1fr));
            gap:20px;margin-top:16px}
.plot-card{background:#FFF;border:1px solid #E2E8F0;border-radius:12px;
           overflow:hidden;box-shadow:0 1px 4px rgba(0,0,0,.06)}
.plot-title{font-size:.8rem;font-weight:600;color:#475569;padding:12px 16px;
            background:#F8FAFC;border-bottom:1px solid #E2E8F0;
            text-transform:uppercase;letter-spacing:.07em}
.plot-card iframe{display:block;width:100%;border:none}

/* ─── Scroll reveal ──────────────────────────────────────────────────────── */
/* NOTE: content is always visible (no opacity:.85 on load)
   Animation uses transform only so nothing is ever hidden */
.reveal{transform:translateY(16px);
        transition:transform .55s ease,opacity .55s ease;
        opacity:.85}
.reveal.visible{transform:none;opacity:1}

/* ─── Footer ─────────────────────────────────────────────────────────────── */
.footer{background:#1E293B;color:#94A3B8;
        padding:32px 56px;font-size:.8rem;
        display:flex;justify-content:space-between;
        flex-wrap:wrap;gap:12px;margin-top:48px}
.footer strong{color:#CBD5E1}

/* ─── Responsive iframe containers ──────────────────────────────────────── */
.iframe-wrap{border-radius:8px;overflow:hidden;border:1px solid #E2E8F0;
             background:#F8FAFC}

/* ─── Print ──────────────────────────────────────────────────────────────── */
@media print{
  .cover{-webkit-print-color-adjust:exact;print-color-adjust:exact}
  .section-divider{margin:24px 0}
  .reveal{transform:none!important;opacity:1!important}
  .toc{display:none}
}
"#.to_string()
}
