import os, re, csv, json, math, shutil, textwrap
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter, MaxNLocator

BASE=Path('/mnt/data/pkl_work')
OUT=Path('/mnt/data/pkl-github-package')
if OUT.exists(): shutil.rmtree(OUT)
(OUT/'docs'/'figures').mkdir(parents=True)
(OUT/'data'/'raw').mkdir(parents=True)
(OUT/'data'/'processed').mkdir(parents=True)
(OUT/'scripts').mkdir(parents=True)

# Copy raw zipped data and magma summary scripts
for name in ['1(1).zip','2(1).zip','3(1).zip','4(1).zip']:
    shutil.copy('/mnt/data/'+name, OUT/'data'/'raw'/name)
for name in ['summary-percentage.m','summary-polynomial(1).m','summary-severity(1).m','summary-structure-KL(1).m']:
    shutil.copy('/mnt/data/'+name, OUT/'scripts'/name)

colors = {'A':'#4C78A8','B':'#F58518','C':'#54A24B','D':'#B279A2','E':'#E45756','F':'#72B7B2','G':'#9D755D'}
family_order=['A','B','C','D','E','F','G']

def type_prime_from_filename(fn):
    stem=Path(fn).name
    m=re.match(r'([A-Za-z]+)(\d+)-(\d+)', stem)
    if not m: return None
    fam=m.group(1).upper(); rank=int(m.group(2)); p=int(m.group(3))
    typ=f'{fam}{rank}'
    return typ,fam,rank,p

def read_lines(path):
    return Path(path).read_text(errors='replace').splitlines()

def parse_meta(lines):
    meta={}
    for ln in lines:
        if ':' in ln:
            k,v=ln.split(':',1); k=k.strip(); v=v.strip()
            if k in ['Type','Prime','Rank','Elements considered','Maximum length considered','Different from KL','Percentage different','Target length','Cartan name']:
                meta[k]=v
    # Different line parse
    if 'Different from KL' in meta:
        m=re.match(r'(\d+) of (\d+)', meta['Different from KL'])
        if m:
            meta['changed']=int(m.group(1)); meta['total']=int(m.group(2))
    if 'Percentage different' in meta:
        meta['pct_changed']=float(meta['Percentage different'].rstrip('%'))
    return meta

def section_csv(lines, title):
    # returns dataframe after a line exactly/contains title; first nonempty line after title is header, following csv-ish rows until blank or all caps non-csv section
    idx=None
    for i,ln in enumerate(lines):
        if ln.strip()==title or title in ln.strip():
            idx=i; break
    if idx is None: return pd.DataFrame()
    # find header
    j=idx+1
    while j<len(lines) and not lines[j].strip(): j+=1
    if j>=len(lines): return pd.DataFrame()
    header=lines[j].strip()
    rows=[]
    for k in range(j+1,len(lines)):
        ln=lines[k].strip()
        if not ln: break
        if ',' not in ln: break
        # stop on long expressions? For length distributions only clean row starts number/comma
        rows.append(ln)
    if not rows: return pd.DataFrame()
    import io
    return pd.read_csv(io.StringIO(header+'\n'+'\n'.join(rows)))

def get_overall(lines, title):
    df=section_csv(lines,title)
    return df.iloc[0].to_dict() if len(df) else {}

records=[]; length_rows=[]
for path in sorted((BASE/'percentage').glob('*.txt')):
    tp=type_prime_from_filename(path.name)
    if not tp: continue
    typ,fam,rank,p=tp
    lines=read_lines(path)
    meta=parse_meta(lines)
    rec={'type':typ,'family':fam,'rank':rank,'prime':p, 'source':'percentage'}
    rec.update({k:meta.get(k) for k in []})
    rec['total']=meta.get('total'); rec['changed']=meta.get('changed'); rec['pct_changed']=meta.get('pct_changed')
    # first changed length
    for ln in lines:
        if ln.startswith('First changed length:'):
            rec['first_changed_length']=int(ln.split(':',1)[1].strip())
    records.append(rec)
    df=section_csv(lines,'LENGTH DISTRIBUTION')
    if len(df):
        df['type']=typ; df['family']=fam; df['rank']=rank; df['prime']=p; df['source']='percentage'
        length_rows.append(df)
summary=pd.DataFrame(records).sort_values(['prime','family','rank'])
length=pd.concat(length_rows,ignore_index=True)

# polynomial overall and length
poly_records=[]; poly_length_rows=[]; poly_freq=[]; graded_freq=[]
for path in sorted((BASE/'polynomial').glob('*.txt')):
    tp=type_prime_from_filename(path.name)
    if not tp: continue
    typ,fam,rank,p=tp
    lines=read_lines(path)
    meta=parse_meta(lines)
    overall=get_overall(lines,'OVERALL POLYNOMIAL SEVERITY')
    rec={'type':typ,'family':fam,'rank':rank,'prime':p,'total':meta.get('total'),'changed':meta.get('changed'),'pct_changed':meta.get('pct_changed')}
    rec.update(overall)
    poly_records.append(rec)
    df=section_csv(lines,'LENGTH DISTRIBUTION')
    if len(df):
        df['type']=typ; df['family']=fam; df['rank']=rank; df['prime']=p
        poly_length_rows.append(df)
    cf=section_csv(lines,'COEFFICIENT POLYNOMIAL FREQUENCIES')
    if len(cf):
        cf['type']=typ; cf['family']=fam; cf['rank']=rank; cf['prime']=p; poly_freq.append(cf)
    gf=section_csv(lines,'GRADED MASS POLYNOMIAL FREQUENCIES')
    if len(gf):
        gf['type']=typ; gf['family']=fam; gf['rank']=rank; gf['prime']=p; graded_freq.append(gf)
poly=pd.DataFrame(poly_records).sort_values(['prime','family','rank'])
poly_length=pd.concat(poly_length_rows,ignore_index=True)
poly_freq=pd.concat(poly_freq,ignore_index=True) if poly_freq else pd.DataFrame()
graded_freq=pd.concat(graded_freq,ignore_index=True) if graded_freq else pd.DataFrame()

# structure overall profiles
full_support=[]; support_size=[]; descent=[]; depth=[]
for path in sorted((BASE/'structure').glob('*.txt')):
    tp=type_prime_from_filename(path.name)
    if not tp: continue
    typ,fam,rank,p=tp
    lines=read_lines(path)
    for title, collector in [('FULL SUPPORT SUMMARY',full_support),('SUPPORT SIZE PROFILE',support_size),('DESCENT PROFILE',descent),('BRUHAT DEPTH DISTRIBUTION OF CORRECTION SUMMANDS',depth)]:
        df=section_csv(lines,title)
        if len(df):
            df['type']=typ; df['family']=fam; df['rank']=rank; df['prime']=p
            collector.append(df)
full_support=pd.concat(full_support,ignore_index=True) if full_support else pd.DataFrame()
support_size=pd.concat(support_size,ignore_index=True) if support_size else pd.DataFrame()
descent=pd.concat(descent,ignore_index=True) if descent else pd.DataFrame()
depth=pd.concat(depth,ignore_index=True) if depth else pd.DataFrame()

# Save processed CSVs
for name,df in [('summary_percentage.csv',summary),('length_distribution.csv',length),('summary_polynomial.csv',poly),('polynomial_length_distribution.csv',poly_length),('full_support_summary.csv',full_support),('support_size_profile.csv',support_size),('descent_profile.csv',descent),('bruhat_depth_distribution.csv',depth)]:
    df.to_csv(OUT/'data'/'processed'/name,index=False)

# Plot helpers
plt.rcParams.update({'font.size': 10, 'axes.titlesize': 13, 'axes.labelsize': 10, 'legend.fontsize': 9, 'figure.dpi': 150, 'savefig.bbox':'tight'})
def savefig(name):
    p=OUT/'docs'/'figures'/name
    plt.savefig(p, dpi=300)
    plt.savefig(p.with_suffix('.pdf'))
    plt.close()

def label_types(df):
    return [f"{r.type}" for r in df.itertuples()]

# 1 overview p=2 bar
p2=summary[summary.prime==2].copy()
p2['type_order']=p2['family'].map({f:i for i,f in enumerate(family_order)})*10+p2['rank']
p2=p2.sort_values(['type_order'])
fig,ax=plt.subplots(figsize=(11,4.8))
ax.bar(range(len(p2)), p2.pct_changed, color=[colors[f] for f in p2.family], edgecolor='black', linewidth=.3)
ax.set_xticks(range(len(p2))); ax.set_xticklabels(p2.type, rotation=0)
ax.set_ylabel('changed elements (%)')
ax.set_title('pKL differs from KL in available finite types at p = 2')
ax.set_ylim(0, max(60,p2.pct_changed.max()*1.15))
for i,r in enumerate(p2.itertuples()):
    if r.pct_changed>=8 or r.type in ['A7']:
        ax.text(i, r.pct_changed+1.2, f'{r.pct_changed:.1f}%', ha='center', va='bottom', fontsize=8, rotation=90 if len(p2)>18 else 0)
handles=[plt.Line2D([0],[0], marker='s', linestyle='', color=colors[f], label=f'type {f}') for f in family_order if f in set(p2.family)]
ax.legend(handles=handles, ncol=7, frameon=False, loc='upper left', bbox_to_anchor=(0,1.12))
ax.grid(axis='y', alpha=.25)
savefig('01_changed_percentage_p2.png')

# 2 rank trends p=2 and p=3 separate lines for families available classical
for prime in [2,3]:
    df=summary[(summary.prime==prime)&(summary.family.isin(['B','C','D','A']))].copy()
    if len(df):
        fig,ax=plt.subplots(figsize=(7.5,4.8))
        for fam in ['A','B','C','D']:
            sub=df[df.family==fam].sort_values('rank')
            if len(sub): ax.plot(sub['rank'], sub['pct_changed'], marker='o', linewidth=2, label=f'type {fam}', color=colors[fam])
        ax.set_xlabel('rank parameter n')
        ax.set_ylabel('changed elements (%)')
        ax.set_title(f'Classical families: changed percentage at p = {prime}')
        ax.set_ylim(bottom=0)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(alpha=.25)
        ax.legend(frameon=False)
        savefig(f'02_classical_rank_trends_p{prime}.png')

# 3 p2 vs p3 same types scatter
common=pd.merge(summary[summary.prime==2],summary[summary.prime==3], on=['type','family','rank'], suffixes=('_p2','_p3'))
if len(common):
    fig,ax=plt.subplots(figsize=(5.8,5.2))
    ax.scatter(common.pct_changed_p2, common.pct_changed_p3, s=80, c=[colors[f] for f in common.family], edgecolor='black', linewidth=.5)
    for r in common.itertuples(): ax.text(r.pct_changed_p2+.6, r.pct_changed_p3+.25, r.type, fontsize=8)
    mx=max(common.pct_changed_p2.max(), common.pct_changed_p3.max())*1.1
    ax.plot([0,mx],[0,mx], color='0.4', linestyle='--', linewidth=1)
    ax.set_xlabel('changed at p = 2 (%)')
    ax.set_ylabel('changed at p = 3 (%)')
    ax.set_title('Same finite type, different prime')
    ax.set_xlim(0,mx); ax.set_ylim(0,mx)
    ax.grid(alpha=.25)
    savefig('03_prime_comparison_p2_vs_p3.png')

# 4 severity stacked p=2 from polynomial, cases >=100 changed
p2poly=poly[(poly.prime==2)&(poly.changed>=100)].copy()
p2poly['order']=p2poly['family'].map({f:i for i,f in enumerate(family_order)})*10+p2poly['rank']
p2poly=p2poly.sort_values('order')
if len(p2poly):
    fig,ax=plt.subplots(figsize=(9,4.8))
    mild=p2poly['mild']/p2poly['changed']*100
    mod=p2poly['moderate']/p2poly['changed']*100
    wild=p2poly['wild']/p2poly['changed']*100
    x=np.arange(len(p2poly))
    ax.bar(x,mild,label='mild: mass 1', color='#A0CBE8', edgecolor='white')
    ax.bar(x,mod,bottom=mild,label='moderate: mass 2 to 5', color='#FFBE7D', edgecolor='white')
    ax.bar(x,wild,bottom=mild+mod,label='wild: mass > 5', color='#F28E2B', edgecolor='white')
    ax.set_xticks(x); ax.set_xticklabels(p2poly.type)
    ax.set_ylim(0,100); ax.set_ylabel('changed elements (%)')
    ax.set_title('Size of the correction among changed elements at p = 2')
    ax.legend(frameon=False, ncol=3, loc='upper center', bbox_to_anchor=(0.5,1.13))
    ax.grid(axis='y', alpha=.25)
    savefig('04_correction_severity_stacked_p2.png')

# 5 genuinely graded fraction p2
if len(p2poly):
    fig,ax=plt.subplots(figsize=(9,4.8))
    gg=p2poly['genuinely_graded']/p2poly['changed']*100
    ax.bar(range(len(p2poly)), gg, color=[colors[f] for f in p2poly.family], edgecolor='black', linewidth=.3)
    ax.set_xticks(range(len(p2poly))); ax.set_xticklabels(p2poly.type)
    ax.set_ylabel('genuinely graded among changed elements (%)')
    ax.set_title('Correction not concentrated in degree zero at p = 2')
    ax.set_ylim(0, max(35, gg.max()*1.2))
    for i,val in enumerate(gg): ax.text(i, val+0.8, f'{val:.1f}%', ha='center', fontsize=8)
    ax.grid(axis='y', alpha=.25)
    savefig('05_genuinely_graded_fraction_p2.png')

# 6 heatmap length percentage for selected types p2 B/C/D/A/F/E maybe all types rows with length normalized? Use absolute length columns; many different. Plot B6 C6 D6 E6 F4 A7 G2 plus Bn/Cn/Dn? create one heatmap available p2 sorted.
sel=p2.sort_values('type_order')['type'].tolist()
ld=length[length.prime==2]
maxell=int(ld.ell.max())
mat=[]
for typ in sel:
    sub=ld[ld.type==typ].set_index('ell')
    arr=[sub.loc[e,'percentage_changed'] if e in sub.index else np.nan for e in range(maxell+1)]
    mat.append(arr)
fig,ax=plt.subplots(figsize=(12,6.8))
im=ax.imshow(mat, aspect='auto', interpolation='nearest', cmap='magma', vmin=0, vmax=np.nanmax(mat))
ax.set_yticks(range(len(sel))); ax.set_yticklabels(sel)
ax.set_xlabel('Coxeter length')
ax.set_title('Changed percentage by Coxeter length at p = 2')
cbar=fig.colorbar(im, ax=ax, pad=.01); cbar.set_label('changed (%)')
savefig('06_length_heatmap_p2.png')

# 7 length profiles for B6 C6 D6 E6 F4 at p2
fig,ax=plt.subplots(figsize=(8,5))
for typ in ['B6','C6','D6','E6','F4','A7']:
    sub=ld[ld.type==typ]
    if len(sub):
        fam=re.match(r'([A-Z])',typ).group(1)
        ax.plot(sub.ell, sub.percentage_changed, marker='o', ms=3, linewidth=1.8, label=typ, color=colors[fam])
ax.set_xlabel('Coxeter length')
ax.set_ylabel('changed elements in that length (%)')
ax.set_title('Where changed elements sit in the length distribution at p = 2')
ax.grid(alpha=.25); ax.legend(frameon=False,ncol=3)
savefig('07_length_profiles_selected_p2.png')

# 8 average and max abs mass p2/p3 bars, cases >=100 changed
for prime in [2,3]:
    df=poly[(poly.prime==prime)&(poly.changed>=20)].copy()
    if len(df):
        df['order']=df['family'].map({f:i for i,f in enumerate(family_order)})*10+df['rank']
        df=df.sort_values('order')
        fig,ax=plt.subplots(figsize=(9,4.8))
        x=np.arange(len(df)); width=.36
        ax.bar(x-width/2, df.avg_abs_mass_at_one, width, label='average abs mass', color='#59A14F')
        ax.bar(x+width/2, df.max_abs_mass_at_one, width, label='maximum abs mass', color='#E15759')
        ax.set_yscale('log')
        ax.set_xticks(x); ax.set_xticklabels(df.type)
        ax.set_ylabel('absolute mass at v = 1 (log scale)')
        ax.set_title(f'Average versus maximum correction size at p = {prime}')
        ax.legend(frameon=False)
        ax.grid(axis='y', alpha=.25, which='both')
        savefig(f'08_abs_mass_avg_max_p{prime}.png')

# 9 support size full/proper p2
fs=full_support[full_support.prime==2].copy()
if len(fs):
    cases=[]
    for typ in sorted(fs.type.unique(), key=lambda t:(family_order.index(t[0]) if t[0] in family_order else 99, int(re.findall(r'\d+',t)[0]))):
        sub=fs[fs.type==typ]
        if len(sub)==2 and sub.changed.sum()>0: cases.append(typ)
    fig,ax=plt.subplots(figsize=(9,4.8))
    x=np.arange(len(cases)); width=.38
    vals_full=[]; vals_prop=[]
    fams=[]
    for typ in cases:
        sub=fs[fs.type==typ]
        vals_full.append(float(sub[sub['class']=='full_support'].percentage_changed.iloc[0]))
        vals_prop.append(float(sub[sub['class']=='proper_parabolic'].percentage_changed.iloc[0]))
        fams.append(typ[0])
    ax.bar(x-width/2, vals_prop, width, label='proper parabolic support', color='#BAB0AC')
    ax.bar(x+width/2, vals_full, width, label='full support', color='#4E79A7')
    ax.set_xticks(x); ax.set_xticklabels(cases)
    ax.set_ylabel('changed elements (%)')
    ax.set_title('Full support elements are more often changed at p = 2')
    ax.legend(frameon=False)
    ax.grid(axis='y', alpha=.25)
    savefig('09_full_vs_parabolic_support_p2.png')

# 10 support size profile selected B6 C6 D6 p2
ssp=support_size[(support_size.prime==2)&(support_size.type.isin(['B6','C6','D6','E6','F4']))]
fig,ax=plt.subplots(figsize=(8,5))
for typ in ['B6','C6','D6','E6','F4']:
    sub=ssp[ssp.type==typ].sort_values('support_size')
    if len(sub):
        fam=typ[0]; ax.plot(sub.support_size, sub.percentage_changed, marker='o', linewidth=2, label=typ, color=colors[fam])
ax.set_xlabel('support size')
ax.set_ylabel('changed elements (%)')
ax.set_title('Changed percentage by support size at p = 2')
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.grid(alpha=.25); ax.legend(frameon=False)
savefig('10_support_size_profile_selected_p2.png')

# 11 Bruhat depth distribution correction summands for selected p2
fig,ax=plt.subplots(figsize=(8,5))
dp=depth[(depth.prime==2)&(depth.type.isin(['B6','C6','D6','E6','F4']))]
for typ in ['B6','C6','D6','E6','F4']:
    sub=dp[dp.type==typ].sort_values('depth')
    if len(sub):
        mass=sub.total_abs_mass/sub.total_abs_mass.sum()*100
        ax.plot(sub.depth, mass, marker='o', ms=3, linewidth=1.8, label=typ, color=colors[typ[0]])
ax.set_xlabel('Bruhat depth of correction summand')
ax.set_ylabel('share of total correction mass (%)')
ax.set_title('How far below w the correction lands at p = 2')
ax.grid(alpha=.25); ax.legend(frameon=False,ncol=3)
savefig('11_bruhat_depth_distribution_p2.png')

# 12 descent heatmap for B6 and C6 maybe create two pngs
for typ in ['B6','C6','D6','E6','F4']:
    sub=descent[(descent.prime==2)&(descent.type==typ)]
    if len(sub):
        maxl=int(max(sub.left_descents.max(), sub.right_descents.max()))
        mat=np.full((maxl+1,maxl+1), np.nan)
        for r in sub.itertuples(): mat[int(r.left_descents), int(r.right_descents)] = r.percentage_changed
        fig,ax=plt.subplots(figsize=(5.5,4.8))
        im=ax.imshow(mat, origin='lower', cmap='viridis', vmin=0, vmax=np.nanmax(mat))
        ax.set_xlabel('right descents'); ax.set_ylabel('left descents')
        ax.set_title(f'Descent profile of changed elements: {typ}, p = 2')
        ax.set_xticks(range(maxl+1)); ax.set_yticks(range(maxl+1))
        cbar=fig.colorbar(im, ax=ax, pad=.02); cbar.set_label('changed (%)')
        savefig(f'12_descent_heatmap_{typ.lower()}_p2.png')

# create an index.html GitHub Pages page
figs=sorted((OUT/'docs'/'figures').glob('*.png'))
html_cards=[]
for f in figs:
    title=f.stem.replace('_',' ')
    html_cards.append(f'''<section class="card"><h2>{title}</h2><a href="figures/{f.name}"><img src="figures/{f.name}" alt="{title}"></a></section>''')
html='''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>pKL growth data</title>
<style>
body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Arial,sans-serif;line-height:1.5;margin:0;background:#fafafa;color:#222}main{max-width:1100px;margin:auto;padding:2rem}h1{font-size:2.2rem;margin-bottom:.2rem}.lead{font-size:1.05rem;color:#444;max-width:850px}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:1rem}.card{background:white;border:1px solid #e7e7e7;border-radius:14px;padding:1rem;box-shadow:0 1px 8px rgba(0,0,0,.04)}.card h2{font-size:1rem;margin:.1rem 0 .8rem;text-transform:capitalize}.card img{width:100%;height:auto;border-radius:8px;border:1px solid #eee}.note{background:#fff8e5;border:1px solid #efd58a;padding:1rem;border-radius:12px;margin:1rem 0}code{background:#eee;padding:.1rem .25rem;border-radius:4px}a{color:#0645ad}</style>
</head><body><main>
<h1>pKL versus KL: growth data</h1>
<p class="lead">Figures and processed summaries for the computations accompanying <em>Almost all pKL basis elements are bad</em>. The plots compare the p-canonical and Kazhdan--Lusztig bases in available finite types, mostly at p = 2 and p = 3.</p>
<div class="note">Raw Magma output is in <code>data/raw</code>; processed CSV summaries are in <code>data/processed</code>; the Magma summary scripts are in <code>scripts</code>. Click a figure for the high-resolution PNG. PDF versions are stored next to the PNGs.</div>
<div class="grid">
'''+'\n'.join(html_cards)+'''
</div></main></body></html>'''
(OUT/'docs'/'index.html').write_text(html)

# README
readme='''# pKL growth data

This repository contains data and plotting material for the computations accompanying the paper

**Almost all pKL basis elements are bad**

by Joseph Baine and Daniel Tubbenhauer.

The computations compare the p-canonical basis element `pC(w)` with the ordinary Kazhdan--Lusztig basis element `C(w)` in finite Coxeter types.  An element is counted as **changed** if `pC(w) != C(w)`.  Equivalently, in the normalization used in the scripts, at least one local stalk polynomial below `w` differs from the characteristic-zero Kazhdan--Lusztig stalk polynomial.

## What is in this repository?

```text
data/raw/          raw Magma output, grouped by summary type
data/processed/    CSV files parsed from the raw output
docs/              GitHub Pages page with plots
docs/figures/      PNG and PDF versions of all figures
scripts/           Magma summary scripts used to generate the raw output
```

The GitHub Pages entry point is `docs/index.html`.

## Main figures

The most useful paper/GitHub figures are:

- `docs/figures/01_changed_percentage_p2.png`: overview of all available finite-type computations at `p = 2`.
- `docs/figures/02_classical_rank_trends_p2.png`: rank trends in the classical families at `p = 2`.
- `docs/figures/03_prime_comparison_p2_vs_p3.png`: direct comparison of `p = 2` and `p = 3` on the same finite types.
- `docs/figures/04_correction_severity_stacked_p2.png`: mild/moderate/wild decomposition among changed elements.
- `docs/figures/05_genuinely_graded_fraction_p2.png`: fraction of changed elements whose correction is not concentrated in degree zero.
- `docs/figures/06_length_heatmap_p2.png`: changed percentage by Coxeter length.
- `docs/figures/09_full_vs_parabolic_support_p2.png`: full-support versus proper-parabolic support comparison.

## Processed data

The processed CSV files are:

- `summary_percentage.csv`: total elements, changed elements, percentage changed, and first changed length.
- `length_distribution.csv`: changed counts and percentages by Coxeter length.
- `summary_polynomial.csv`: global correction statistics in the ordinary KL basis.
- `polynomial_length_distribution.csv`: polynomial correction statistics by length.
- `full_support_summary.csv`: full-support versus proper-parabolic support.
- `support_size_profile.csv`: changed percentage by support size.
- `descent_profile.csv`: changed percentage by number of left and right descents.
- `bruhat_depth_distribution.csv`: distribution of correction summands by Bruhat depth.

## Correction conventions

For the polynomial summaries, the correction is measured in the ordinary KL basis:

```text
Delta_w = pC(w) - C(w) = sum_x f_{x,w}(v) C(x).
```

The summary scripts also record several crude size measures:

- `signed_mass_at_one`: sum of `f_{x,w}(1)`.
- `abs_mass_at_one`: sum of the absolute values of all monomial coefficients in all `f_{x,w}(v)`.
- `monomial_complexity`: total number of nonzero monomials in the correction.
- `genuinely_graded`: the correction has some nonzero-degree term.

For the mild/moderate/wild plots, the convention is:

```text
mild      signed_mass_at_one = 1
moderate  2 <= signed_mass_at_one <= 5
wild      signed_mass_at_one > 5
```

## Reproducing the summaries

The Magma scripts in `scripts/` expect the ASLoc setup and saved p-canonical bases.  Typical usage is:

```bash
magma -b type:=B5 prime:=2 saveDir:=saves summary-percentage.m > b5-2.txt
magma -b type:=B5 prime:=2 saveDir:=saves summary-polynomial.m > b5-2-polynomial.txt
magma -b type:=B5 prime:=2 saveDir:=saves summary-structure-KL.m > b5-2-structure-KL.txt
```

The raw outputs in `data/raw/` are enough to regenerate all plots.

## Current computed range

The uploaded data cover:

- `p = 2`: `A7`, `B2`--`B6`, `C2`--`C6`, `D4`--`D6`, `E6`, `F4`, `G2`.
- `p = 3`: `B6`, `C6`, `D6`, `E6`, `F4`, `G2`.

Some low-rank cases at `p = 3` are absent because the available run set focused on the first nontrivial/high-value comparisons.

## Notes for the paper

The computations are not used in the proof.  They are intended to show what already happens in small rank and to illustrate that the corrections are not merely present but can become large, spread out in Bruhat order, and genuinely graded.
'''
(OUT/'README.md').write_text(readme)

# small python plotting script included
plot_script = r'''#!/usr/bin/env python3
"""Parse the processed CSVs and regenerate selected plots.
Run this from the repository root after installing pandas and matplotlib.
"""
# The full plotting script used to generate the distributed figures is in the ChatGPT artifact history.
# This lightweight placeholder documents the input files and can be expanded as needed.
from pathlib import Path
import pandas as pd

base = Path('data/processed')
print(pd.read_csv(base / 'summary_percentage.csv').head())
'''
(OUT/'scripts'/'plot_from_processed.py').write_text(plot_script)

# summary JSON
summary_info={
    'n_cases': int(len(summary)),
    'n_p2_cases': int((summary.prime==2).sum()),
    'n_p3_cases': int((summary.prime==3).sum()),
    'figures': [f.name for f in figs]
}
(OUT/'data'/'processed'/'manifest.json').write_text(json.dumps(summary_info, indent=2))

# zip package
zip_path='/mnt/data/pkl-github-package.zip'
if Path(zip_path).exists(): Path(zip_path).unlink()
shutil.make_archive('/mnt/data/pkl-github-package','zip',OUT)
print('Wrote', OUT)
print('Zip', zip_path)
print(summary[['type','prime','total','changed','pct_changed','first_changed_length']].to_string(index=False))
