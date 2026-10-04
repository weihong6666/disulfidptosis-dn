#!/usr/bin/env python3
"""
Generate Figure 1: Complete analytical workflow diagram for
Disulfidptosis Biomarkers in Diabetic Nephropathy study.

Publication-quality scientific schematic with blue/teal color scheme.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

# ── Configuration ──────────────────────────────────────────────────
OUTPUT_PATH = r"C:\Users\weihong\Desktop\生物竞赛\Fig1_study_design.png"
DPI = 300
FIG_W, FIG_H = 14, 20  # inches

# Color palette (blue/teal)
C_STAGE1 = '#E3F2FD'   # lightest blue
C_STAGE2 = '#BBDEFB'   # light blue
C_STAGE3 = '#80CBC4'   # teal
C_STAGE4 = '#4DB6AC'   # medium teal
C_STAGE5 = '#00838F'   # deep teal
C_HEADER = '#263238'   # dark header bg
C_HEADER_TEXT = '#FFFFFF'
C_BORDER = '#37474F'
C_ARROW = '#546E7A'
C_SUBBOX = '#E0F2F1'
C_SUBBOX_BORDER = '#00695C'
C_ACCENT1 = '#1565C0'  # blue accent
C_ACCENT2 = '#00838F'  # teal accent
C_WHITE = '#FFFFFF'
C_TEXT = '#212121'
C_SUBTEXT = '#546E7A'

# ── Helper Functions ───────────────────────────────────────────────

def draw_rounded_box(ax, xy, width, height, facecolor, edgecolor,
                     linewidth=1.5, corner_radius=0.15, zorder=1):
    """Draw a rounded rectangle box."""
    box = FancyBboxPatch(
        xy, width, height,
        boxstyle=f"round,pad={corner_radius}",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=linewidth,
        zorder=zorder
    )
    ax.add_patch(box)
    return box


def draw_stage_box(ax, x, y, w, h, title, content_lines, header_color,
                   body_color, border_color, title_color=C_WHITE):
    """Draw a stage box with header and content area."""
    # Body
    draw_rounded_box(ax, (x, y), w, h, body_color, border_color,
                     linewidth=2.0, corner_radius=0.12, zorder=2)

    # Header bar (flat top section painted over)
    header_h = 0.55
    header = FancyBboxPatch(
        (x + 0.08, y + h - header_h - 0.04), w - 0.16, header_h,
        boxstyle=f"round,pad=0.06",
        facecolor=header_color,
        edgecolor='none',
        zorder=3
    )
    ax.add_patch(header)

    ax.text(x + w / 2, y + h - header_h / 2 - 0.04, title,
            ha='center', va='center', fontsize=11, fontweight='bold',
            color=title_color, zorder=4, fontfamily='sans-serif')

    # Content
    for i, line in enumerate(content_lines):
        ax.text(x + 0.35, y + h - header_h - 0.65 - i * 0.45,
                line, ha='left', va='top', fontsize=8.0,
                color=C_TEXT, zorder=4, fontfamily='sans-serif')


def draw_arrow(ax, start_xy, end_xy, color=C_ARROW, lw=2.5):
    """Draw a downward arrow between stages."""
    ax.annotate('', xy=end_xy, xytext=start_xy,
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=lw, connectionstyle='arc3,rad=0'),
                zorder=5)


def draw_stage_number(ax, x, y, number):
    """Draw a circled stage number."""
    circle = plt.Circle((x, y), 0.28, facecolor=C_ACCENT1,
                         edgecolor='white', linewidth=1.5, zorder=6)
    ax.add_patch(circle)
    ax.text(x, y, str(number), ha='center', va='center',
            fontsize=12, fontweight='bold', color='white', zorder=7)


# ── Build Figure ───────────────────────────────────────────────────

fig, ax = plt.subplots(1, 1, figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 14)
ax.set_ylim(0, 20)
ax.set_aspect('equal')
ax.axis('off')

# Background
fig.patch.set_facecolor('#FAFAFA')
ax.set_facecolor('#FAFAFA')

# ── Title ──────────────────────────────────────────────────────────

ax.text(7, 19.3, 'Study Design: Disulfidptosis-Related Biomarkers',
        ha='center', va='center', fontsize=16, fontweight='bold',
        color=C_TEXT, fontfamily='sans-serif')
ax.text(7, 18.6, 'in Diabetic Nephropathy',
        ha='center', va='center', fontsize=14, fontweight='normal',
        color=C_SUBTEXT, fontfamily='sans-serif')

# Divider line
ax.plot([2.5, 11.5], [18.35, 18.35], color=C_ARROW, lw=1, zorder=1)

# ── Stage dimensions ───────────────────────────────────────────────
BOX_X = 1.5
BOX_W = 11.0
GAP_BETWEEN = 0.58
START_Y = 17.5

# ── Stage 1: Data Acquisition ──────────────────────────────────────
S1_Y = 14.3
S1_H = 3.0

draw_stage_box(ax, BOX_X, S1_Y, BOX_W, S1_H,
    'Stage 1: Data Acquisition',
    [
        '▶ Training set: GSE96804 (41 DN + 20 Control, glomeruli)',
        '▶ Validation set 1: GSE30528 (9 DN + 13 Ctrl, glomeruli)',
        '▶ Validation set 2: GSE30529 (10 DN + 12 Ctrl, tubulointerstitium)',
        '▶ Single-cell: GSE131882 (3 DN + 3 Control, kidney)',
    ],
    C_ACCENT1, C_STAGE1, C_BORDER)

draw_stage_number(ax, BOX_X + 0.3, S1_Y + S1_H - 0.4, 1)

# Arrow 1->2
draw_arrow(ax, (7, S1_Y - 0.15), (7, S1_Y - GAP_BETWEEN + 0.15))

# ── Stage 2: Differential Expression & WGCNA ───────────────────────
S2_Y = 10.65
S2_H = 2.95

draw_stage_box(ax, BOX_X, S2_Y, BOX_W, S2_H,
    'Stage 2: Differential Expression & WGCNA',
    [
        '▶ limma differential analysis: |log₂FC| > 0.5, adj.P < 0.05',
        '▶ Disulfidptosis gene set (35 genes, Liu et al. 2023)',
        '▶ Intersection → 9 DR-DEGs (disulfidptosis-related DEGs)',
        '▶ WGCNA co-expression network: β = 10, R² = 0.884, 12 modules',
    ],
    C_ACCENT1, C_STAGE2, C_BORDER)

draw_stage_number(ax, BOX_X + 0.3, S2_Y + S2_H - 0.4, 2)

# Arrow 2->3
draw_arrow(ax, (7, S2_Y - 0.15), (7, S2_Y - GAP_BETWEEN + 0.15))

# ── Stage 3: Machine Learning Biomarker Selection ──────────────────
S3_Y = 7.0
S3_H = 2.95

draw_stage_box(ax, BOX_X, S3_Y, BOX_W, S3_H,
    'Stage 3: Machine Learning Biomarker Selection',
    [
        '▶ LASSO logistic regression (10-fold cross-validation, lambda.min)',
        '▶ SVM-RFE (linear kernel, recursive feature elimination)',
        '▶ Random Forest (500 trees, mean decrease Gini)',
        '▶ Triple intersection → 7 final biomarkers',
        '     NCKAP1L, THSD7A, MYH10, PDLIM1, IL1B, FLNB, SLC3A2',
    ],
    C_ACCENT1, C_STAGE3, C_BORDER)

draw_stage_number(ax, BOX_X + 0.3, S3_Y + S3_H - 0.4, 3)

# Arrow 3->4
draw_arrow(ax, (7, S3_Y - 0.15), (7, S3_Y - GAP_BETWEEN + 0.15))

# ── Stage 4: Diagnostic Model & Validation ─────────────────────────
S4_Y = 3.8
S4_H = 2.5

draw_stage_box(ax, BOX_X, S4_Y, BOX_W, S4_H,
    'Stage 4: Diagnostic Model & Validation',
    [
        '▶ Logistic regression-based diagnostic Nomogram',
        '▶ ROC evaluation: training set (GSE96804) + 2 external validation sets',
        '▶ Calibration (1000× bootstrap) + DCA decision curve analysis',
    ],
    C_ACCENT1, C_STAGE4, C_BORDER)

draw_stage_number(ax, BOX_X + 0.3, S4_Y + S4_H - 0.4, 4)

# Arrow 4->5
draw_arrow(ax, (7, S4_Y - 0.15), (7, S4_Y - GAP_BETWEEN + 0.15))

# ── Stage 5: Multi-Dimensional Downstream Analysis ─────────────────
S5_Y = 0.3
S5_H = 3.0

# Main box
draw_rounded_box(ax, (BOX_X, S5_Y), BOX_W, S5_H, C_STAGE5, C_BORDER,
                 linewidth=2.0, corner_radius=0.12, zorder=2)

# Header
header_h = 0.55
header = FancyBboxPatch(
    (BOX_X + 0.08, S5_Y + S5_H - header_h - 0.04), BOX_W - 0.16, header_h,
    boxstyle="round,pad=0.06",
    facecolor=C_HEADER,
    edgecolor='none',
    zorder=3
)
ax.add_patch(header)
ax.text(7, S5_Y + S5_H - header_h / 2 - 0.04,
        'Stage 5: Multi-Dimensional Downstream Analysis',
        ha='center', va='center', fontsize=11, fontweight='bold',
        color=C_WHITE, zorder=4, fontfamily='sans-serif')

draw_stage_number(ax, BOX_X + 0.3, S5_Y + S5_H - 0.4, 5)

# 6 sub-boxes in a 3x2 grid
sub_items = [
    ('Immune Infiltration',
     'CIBERSORT + ssGSEA\nImmune Checkpoint\nESTIMATE score',
     '❶'),
    ('Single-Cell Analysis',
     'Seurat → CellChat\nmonocle3 pseudotime\n20,681 cells, 10 types',
     '❷'),
    ('Enrichment Analysis',
     'GO + KEGG\nclusterProfiler\nBP, CC, MF, pathways',
     '❸'),
    ('Drug Prediction',
     'DGIdb drug-gene\ninteraction database\nAutoDock Vina docking',
     '❹'),
    ('Regulatory Network',
     'TF → miRNA → ceRNA\nmulti-layer regulatory\nnetwork construction',
     '❺'),
    ('Protein Validation',
     'HPA IHC\nimmunohistochemistry\nprotein-level verification',
     '❻'),
]

sub_w = (BOX_W - 0.8) / 3
sub_h = 1.45
sub_start_x = BOX_X + 0.25
sub_start_y = S5_Y + 0.5

for i, (title, desc, num) in enumerate(sub_items):
    col = i % 3
    row = 1 - (i // 3)  # top row first
    sx = sub_start_x + col * (sub_w + 0.15)
    sy = sub_start_y + row * (sub_h + 0.15)

    # Sub-box
    draw_rounded_box(ax, (sx, sy), sub_w, sub_h,
                     C_SUBBOX, C_SUBBOX_BORDER,
                     linewidth=1.2, corner_radius=0.08, zorder=4)

    # Sub-header
    sub_hdr = FancyBboxPatch(
        (sx + 0.05, sy + sub_h - 0.42), sub_w - 0.10, 0.38,
        boxstyle="round,pad=0.04",
        facecolor=C_ACCENT2,
        edgecolor='none',
        zorder=5
    )
    ax.add_patch(sub_hdr)

    ax.text(sx + sub_w / 2 + 0.1, sy + sub_h - 0.23,
            title, ha='center', va='center',
            fontsize=7.5, fontweight='bold', color=C_WHITE, zorder=6)

    # Number circle for sub-box
    ax.text(sx + 0.18, sy + sub_h - 0.23, num,
            ha='center', va='center', fontsize=7, color=C_WHITE,
            fontweight='bold', zorder=7)

    # Description
    ax.text(sx + sub_w / 2, sy + 0.52, desc,
            ha='center', va='center', fontsize=6.8,
            color=C_TEXT, zorder=6, fontfamily='sans-serif',
            linespacing=1.35)


# ── Left-side annotation bar ───────────────────────────────────────

# Stage type labels on the left
stage_labels = [
    ('Data\nPreprocessing', 15.8, 13.3),
    ('Feature\nSelection', 12.1, 9.65),
    ('Machine\nLearning', 8.5, 6.0),
    ('Model\nValidation', 5.05, 3.3),
    ('Downstream\nAnalysis', 1.8, 0.3),
]

for label, y1, y2 in stage_labels:
    mid_y = (y1 + y2 + 0.5) / 2
    ax.text(0.65, mid_y, label, ha='center', va='center',
            fontsize=7.5, fontweight='bold', color=C_ACCENT1,
            fontfamily='sans-serif', rotation=0, zorder=5)
    # Small bracket line
    ax.plot([1.15, 1.15], [y2 + 0.2, y1 - 0.2], color=C_ACCENT1, lw=1.5, zorder=1)
    ax.plot([1.05, 1.15], [y2 + 0.2, y2 + 0.2], color=C_ACCENT1, lw=1.5, zorder=1)
    ax.plot([1.05, 1.15], [y1 - 0.2, y1 - 0.2], color=C_ACCENT1, lw=1.5, zorder=1)


# ── Key findings callout box ──────────────────────────────────────

callout_x, callout_y = 9.2, 17.5
callout_w, callout_h = 3.5, 1.2

draw_rounded_box(ax, (callout_x, callout_y), callout_w, callout_h,
                 '#E8EAF6', '#3949AB', linewidth=1.2, corner_radius=0.1, zorder=5)
ax.text(callout_x + callout_w / 2, callout_y + callout_h - 0.2,
        'Key Result', ha='center', va='center',
        fontsize=8, fontweight='bold', color='#1A237E', zorder=6)
ax.text(callout_x + callout_w / 2, callout_y + 0.35,
        '7 Disulfidptosis-Related\nDiagnostic Biomarkers\nfor Diabetic Nephropathy',
        ha='center', va='center', fontsize=7.2,
        color=C_TEXT, zorder=6, linespacing=1.3)


# ── Footer ─────────────────────────────────────────────────────────
ax.text(7, 0.05, 'Figure 1. Overview of the study design and analytical workflow.',
        ha='center', va='bottom', fontsize=8.5, fontstyle='italic',
        color=C_SUBTEXT, fontfamily='sans-serif')


# ── Save ───────────────────────────────────────────────────────────

plt.tight_layout(pad=0.5)
fig.savefig(OUTPUT_PATH, dpi=DPI, bbox_inches='tight',
            facecolor=fig.get_facecolor(), edgecolor='none')
plt.close()

print(f"Figure saved to: {OUTPUT_PATH}")
print(f"Dimensions: {FIG_W}x{FIG_H} inches at {DPI} DPI")
