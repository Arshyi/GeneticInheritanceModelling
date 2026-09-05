"""Reproduce numerical experiments and figures; all simulation parameters are explicit."""
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.stats import binom, kstest
from genetics.core import InheritanceModel, hwe, select, mutation, migrate
from genetics.extensions import additive_pmf, mixture_cdf, hwe_exact, linked_gametes, linked_cross, linkage_equilibrium_update
from experiments.reproduce_version1 import legacy_sickle

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'; FIGURES=ROOT/'figures'
SEED=20260905
COLORS=['#136f79','#d18a25','#293c61','#b44949','#79706e']


def dump(name,obj):
    (RESULTS/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def savefig(name):
    plt.savefig(FIGURES/(name+'.png'),dpi=190,bbox_inches='tight')
    plt.savefig(FIGURES/(name+'.svg'),bbox_inches='tight')
    plt.close()


def population_validation():
    path=ROOT/'data/observed_genotypes.csv'
    raw=list(csv.DictReader(path.open(encoding='utf-8-sig')))
    assert len(raw)==26 and len({r['population'] for r in raw})==26
    out=[]
    def evaluate(label,counts,scope):
        n=sum(counts); q=(2*counts[2]+counts[1])/(2*n)
        exp=n*hwe([1-q,q])
        return dict(population=label,scope=scope,n=n,counts=counts,q=q,expected=exp.tolist(),
                    heterozygote_residual=counts[1]-float(exp[1]),
                    hwe_exact_p=hwe_exact(counts),
                    fitted_frequency_total_variation=float(np.abs(np.array(counts)/n-exp/n).sum()/2))
    for r in raw:
        counts=[int(r[k]) for k in ('hom_1','het','hom_2')]
        assert sum(counts)==int(r['n']) and r['locus']=='rs334'
        row=evaluate(r['population'],counts,'component population')
        row['superpopulation']=r['superpopulation']
        out.append(row)
    # Holm family-wise adjustment across all 26 tests, including monomorphic rows.
    order=sorted(range(len(out)),key=lambda i:out[i]['hwe_exact_p'])
    previous=0
    for rank,i in enumerate(order):
        previous=max(previous,min(1,(len(out)-rank)*out[i]['hwe_exact_p']))
        out[i]['holm_adjusted_p_26_tests']=previous
    total=[sum(r['counts'][j] for r in out) for j in range(3)]
    pooled=evaluate('ALL (pooled)',total,'pooled diagnostic only')
    afr=[r for r in out if r['superpopulation']=='AFR']
    african=evaluate('AFR (pooled)',[sum(r['counts'][j] for r in afr) for j in range(3)],'pooled diagnostic only')
    assert total==[2367,137,0] and sum(total)==2504
    summary=dict(source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                 component_populations=out,pooled=pooled,african_pool=african,
                 stratified_expected_rare_homozygotes=sum(r['expected'][2] for r in out),
                 interpretation='Fitted HWE compatibility audit of sampled genomic calls; no clinical labels, world prevalence, or predictive accuracy.')
    dump('population_validation.json',summary)
    with (RESULTS/'population_validation.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['population','superpopulation','n','q','hwe_exact_p','holm_adjusted_p_26_tests','heterozygote_residual','fitted_frequency_total_variation']
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(out)
    displayed=[r for r in out if r['q']>0]
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    positions=np.arange(len(displayed))
    axes[0].bar(positions,[r['q']*100 for r in displayed],color=COLORS[0])
    axes[0].set_xticks(positions,[r['population'] for r in displayed],rotation=45)
    axes[0].set_ylabel('rs334 A allele frequency (%)')
    axes[0].set_title('Observed reference-panel counts')
    vals=[sum(total),pooled['expected'][2],summary['stratified_expected_rare_homozygotes']]
    axes[1].bar(['Observed A/A','Pooled HWE','Within-population HWE'],[total[2],vals[1],vals[2]],color=[COLORS[3],COLORS[2],COLORS[0]])
    axes[1].set_ylabel('Rare-homozygote count')
    axes[1].set_title('Pooling changes fitted expectations')
    axes[1].tick_params(axis='x',labelrotation=15)
    fig.suptitle('Population-model audit: 2,504 calls, 26 populations',fontweight='bold')
    fig.tight_layout();savefig('population_hwe')
    return summary


def generation_scenarios():
    m=InheritanceModel((2,)); initial=np.array([.25,.5,.25])
    scenarios={}
    for label,fitness in [('Neutral random mating',[1,1,1]),('Zero SS reproductive weight',[1,1,0]),('Illustrative heterozygote advantage',[.9,1,.2])]:
        x=initial.copy(); trajectory=[x.tolist()]
        for _ in range(19):
            x=m.next_generation(select(x,fitness));trajectory.append(x.tolist())
        scenarios[label]=trajectory
    legacy=[[float(v) for v in row] for row in legacy_sickle('AS','AS',20)]
    assert np.max(np.abs(np.array(scenarios['Neutral random mating'])-initial))<1e-13
    assert abs(scenarios['Zero SS reproductive weight'][2][2]-.0625)<1e-13
    # Explicit operator-order experiment: selection -> migration -> gametes -> mutation -> random fertilization.
    x=hwe([.8,.2]); combined=[]
    for gen in range(21):
        combined.append(dict(generation=gen,AA=float(x[0]),AS=float(x[1]),SS=float(x[2]),q=float(x[2]+x[1]/2)))
        breeders=migrate(select(x,[1,1,.5]),hwe([.6,.4]),.01)
        q=breeders[2]+breeders[1]/2
        gametes=mutation([1-q,q],[[1-1e-5,1e-5],[1e-6,1-1e-6]])
        x=hwe(gametes)
    data=dict(initial_parental_cross='AS x AS',generation1=initial.tolist(),
              scenarios=scenarios,legacy_screenshot_algorithm=legacy,operator_order_synthetic=combined,
              parameters={'fitness':[1,1,.5],'migration':.01,'immigrant_q':.4,'forward_mutation':1e-5,'reverse_mutation':1e-6},
              status='Synthetic assumption comparisons; parameters not estimated from biological data')
    dump('generation_scenarios.json',data)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for k,(label,trajectory) in enumerate(scenarios.items()):
        axes[0].plot(range(1,21),np.array(trajectory)[:,2]*100,label=label,color=COLORS[k],lw=2)
    axes[0].plot(range(1,21),np.array(legacy)[:,2]*100,'--',label='Version I displayed algorithm',color=COLORS[3])
    axes[0].set(xlabel='Offspring generation',ylabel='SS birth fraction (%)',title='Dynamics depend on a mating/selection law')
    axes[0].legend(fontsize=7)
    axes[1].plot(range(1,21),np.sum(legacy,axis=1),label='Version I displayed algorithm',color=COLORS[3],lw=2)
    axes[1].axhline(1,color=COLORS[0],label='Normalized population model',lw=2)
    axes[1].set(xlabel='Offspring generation',ylabel='Sum of reported genotype fractions',ylim=(.9,1.01),title='Audit of probability conservation')
    axes[1].legend(fontsize=8)
    fig.tight_layout();savefig('generation_models')
    return data


def epistasis_and_linkage():
    m=InheritanceModel((3,2)) # A=0,B=1,O=2; H=0,h=1.
    ao=m.indices[0][(0,2)];bo=m.indices[0][(1,2)];hetero=m.indices[1][(0,1)]
    offspring=m.cross(m.encode((ao,hetero)),m.encode((bo,hetero)))
    phen=defaultdict(float)
    for code,p in offspring.items():
        ag,hg=m.decode(code);g=m.loci[0][ag]
        if m.loci[1][hg]==(1,1):label='O-like (hh masked)'
        elif g==(2,2):label='O (ABO OO)'
        elif 0 in g and 1 in g:label='AB'
        elif 0 in g:label='A'
        else:label='B'
        phen[label]+=p
    assert phen['A']==phen['B']==phen['AB']==3/16
    assert phen['O (ABO OO)']+phen['O-like (hh masked)']==7/16
    # M/N is a separate named Mendelian coding example, no population frequencies used.
    mn=InheritanceModel((2,)).cross(1,1)
    assert mn=={0:.25,1:.5,2:.25}
    coupling=((0,0),(1,1));repulsion=((0,1),(1,0))
    linkage=[]
    for r in (0,.01,.1,.25,.5):
        a,b=linked_gametes(coupling,r),linked_gametes(repulsion,r)
        tv=sum(abs(a.get(g,0)-b.get(g,0)) for g in set(a)|set(b))/2
        linkage.append(dict(r=r,phase_total_variation=tv,coupling={''.join(map(str,g)):p for g,p in a.items()}))
    data=dict(ABO_FUT1_simplified=dict(parents=['AO/Hh','BO/Hh'],phenotype_probabilities=dict(phen),G=m.G,U=m.U),
              MNS_simplified=dict(parents=['MN','MN'],offspring={'MM':mn[0],'MN':mn[1],'NN':mn[2]}),
              linkage=linkage,status='Exact predictions within simplified models; no independent trait phenotype cohort validation')
    dump('additional_systems.json',data)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    axes[0].bar(['A','B','AB','O / O-like'],[3/16,3/16,3/16,7/16],color=[COLORS[0],COLORS[1],COLORS[2],COLORS[3]])
    axes[0].set(ylabel='Offspring phenotype probability',title='Simplified ABO + FUT1 epistasis')
    axes[1].plot([r['r'] for r in linkage],[r['phase_total_variation'] for r in linkage],marker='o',color=COLORS[0])
    axes[1].set(xlabel='Recombination fraction r',ylabel='Total variation between gamete distributions',title='Same unphased genotype; different phase')
    fig.tight_layout();savefig('epistasis_linkage')
    return data


def polygenic():
    n=200;p=[.25,.5,.25];start=time.perf_counter();score=additive_pmf([p]*n);elapsed=time.perf_counter()-start
    np.testing.assert_allclose(score,binom.pmf(np.arange(401),400,.5),atol=1e-15)
    sd=10.
    def quantile(q):return brentq(lambda x:float(mixture_cdf([x],score,sd)[0])-q,-100,500)
    lo,hi=quantile(.05),quantile(.95)
    rng=np.random.default_rng(SEED);N=10000
    # Independent simulator uses Bernoulli-sum binomial sampling, not PMF resampling.
    samples=rng.binomial(2*n,.5,N)+rng.normal(0,sd,N)
    covered=int(((samples>=lo)&(samples<=hi)).sum());coverage=covered/N
    z=1.959963984540054;den=1+z*z/N
    center=(coverage+z*z/(2*N))/den
    half=z*math.sqrt(coverage*(1-coverage)/N+z*z/(4*N*N))/den
    pit=mixture_cdf(samples,score,sd);ks=kstest(pit,'uniform')
    data=dict(n_loci=n,genotype_catalog=str(3**n),score_bins=len(score),pmf_bytes=score.nbytes,
              construction_seconds_single_observation=elapsed,score_mean=float(np.dot(np.arange(401),score)),
              score_variance=float(np.dot((np.arange(401)-200)**2,score)),residual_sd=sd,
              nominal_interval_mass=.9,interval=[lo,hi],independent_synthetic_draws=N,seed=SEED,
              observed_interval_coverage=coverage,coverage_wilson_95=[center-half,center+half],
              pit_KS_statistic=float(ks.statistic),pit_KS_p=float(ks.pvalue),
              predictive_validation_status='Synthetic only: coefficients, independence, and residual distribution stipulated; no human height predictions')
    dump('polygenic_synthetic.json',data)
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    axes[0].plot(np.arange(401),score,color=COLORS[0],lw=2)
    axes[0].set(xlim=(150,250),xlabel='Synthetic additive score',ylabel='Probability',title='200 loci; 401 score bins')
    axes[1].hist(samples,bins=55,density=True,color=COLORS[0],alpha=.65)
    axes[1].axvline(lo,color=COLORS[3],ls='--');axes[1].axvline(hi,color=COLORS[3],ls='--')
    axes[1].set(xlabel='Score + stipulated Gaussian residual',ylabel='Density',title=f'90% interval: measured coverage {coverage:.2%}')
    fig.tight_layout();savefig('polygenic_distribution')
    return data


def run():
    RESULTS.mkdir(exist_ok=True);FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'white'})
    outputs={}
    for fn in [population_validation,generation_scenarios,epistasis_and_linkage,polygenic]:
        outputs[fn.__name__]=fn();print(fn.__name__+' complete',flush=True)
    dump('science_metadata.json',dict(completed_at_utc=datetime.now(timezone.utc).isoformat(),seed=SEED,
         numpy=np.__version__,input_pdf_sha256=hashlib.sha256((ROOT.parent/'Bioinformatics-Arshyia Mehran.pdf').read_bytes()).hexdigest(),
         status='Executed local deterministic computations and seeded synthetic validation'))
    print('All science experiments completed.')


if __name__=='__main__':run()
