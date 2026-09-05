from fractions import Fraction
from itertools import product, combinations_with_replacement
from unittest.mock import patch
import math
import numpy as np
import pytest
from genetics.core import InheritanceModel, ResourceLimitError, hwe, mutation, select, migrate
from genetics.extensions import linked_gametes, linked_cross, additive_pmf, hwe_exact, linkage_equilibrium_update


def oracle(model,i,j):
    """Independent enumeration of all maternal/paternal allele-copy choices."""
    out = {}
    for copies in product(range(2), repeat=2*len(model.loci)):
        genotype = []
        for l,(a,b) in enumerate(zip(model.decode(i),model.decode(j))):
            g = tuple(sorted((model.loci[l][a][copies[2*l]],model.loci[l][b][copies[2*l+1]])))
            genotype.append(model.loci[l].index(g))
        code = model.encode(genotype)
        out[code] = out.get(code,Fraction(0))+Fraction(1,4**len(model.loci))
    return {k:float(v) for k,v in out.items()}


@pytest.mark.parametrize('alleles',[(2,),(3,),(2,2),(3,2),(2,2,2),(4,)])
def test_exhaustive_mendelian_oracle(alleles):
    m = InheritanceModel(alleles)
    nnz=0
    for i,j in combinations_with_replacement(range(m.G),2):
        actual=m.cross(i,j)
        assert actual == oracle(m,i,j)
        assert actual == m.cross(j,i)
        assert sum(actual.values()) == 1
        nnz += len(actual)
    assert nnz == m.nnz


@pytest.mark.parametrize('alleles',[(2,),(3,),(3,2),(2,2,2)])
def test_representation_and_population_agreement(alleles):
    m=InheritanceModel(alleles)
    dense,csr,hashed=(m.kernel(k) for k in ('dense','csr','hash'))
    np.testing.assert_array_equal(dense,csr.toarray())
    np.testing.assert_allclose(dense.sum(axis=0),1)
    assert csr.nnz==m.nnz and not (csr.data==0).any()
    x=np.random.default_rng(234).dirichlet(np.ones(m.G))
    ref=m.next_generation(x)
    for k in (dense,csr,hashed): np.testing.assert_allclose(m.next_generation(x,k),ref,atol=1e-14)
    assert abs(ref.sum()-1)<1e-12


def test_hwe_is_fixed_point_and_alleles_conserved():
    m=InheritanceModel((3,))
    for p in ([.26,.077,.663],[.2,.3,.5],[1.,0,0]):
        np.testing.assert_allclose(m.next_generation(hwe(p)),hwe(p),atol=1e-14)
    x=np.array([.1,.2,.3,.1,.2,.1])
    allele=lambda x:sum((x[i]*np.bincount(g,minlength=3)/2 for i,g in enumerate(m.loci[0])),start=np.zeros(3))
    np.testing.assert_allclose(allele(x),allele(m.next_generation(x)))


def test_encode_roundtrip():
    m=InheritanceModel((3,2,4))
    for c in range(m.G): assert m.encode(m.decode(c))==c


def test_budgets_and_allocation_failure():
    m=InheritanceModel((2,)*12,max_bytes=100_000)
    with pytest.raises(ResourceLimitError): m.kernel('dense')
    with pytest.raises(ResourceLimitError): m.cross(m.encode((1,)*12),m.encode((1,)*12))
    assert m.probability(0,0,0)==1
    small=InheritanceModel((2,))
    with patch('genetics.core.np.zeros',side_effect=MemoryError('injected')):
        with pytest.raises(ResourceLimitError): small.kernel('dense')
    assert small.cross(1,1)=={0:.25,1:.5,2:.25}


def test_underflow_stays_representable_in_log_space():
    m=InheritanceModel((2,)*1000)
    parent=m.encode((1,)*1000)
    assert math.isfinite(m.log_probability(parent,parent,0))
    with pytest.raises(FloatingPointError): m.probability(parent,parent,0)
    with pytest.raises(FloatingPointError): next(m.iter_cross(parent,parent))
    assert math.isfinite(next(m.iter_log_cross(parent,parent))[1])
    assert m.log_probability(0,0,m.G-1)==-math.inf


def test_population_operators():
    np.testing.assert_allclose(select([.25,.5,.25],[1,1,0]),[1/3,2/3,0])
    np.testing.assert_allclose(mutation([1,0],[[.99,.01],[.1,.9]]),[.99,.01])
    np.testing.assert_allclose(migrate([1,0],[0,1],.2),[.8,.2])
    for bad in ([[1,1],[0,1]],[[1,-1],[0,1]]):
        with pytest.raises(ValueError): mutation([1,0],bad)
    with pytest.raises(ValueError): select([.5,.5],[0,0])
    with pytest.raises(ValueError): hwe([.6,.5])
    np.testing.assert_allclose(select([.5,.5],[5e-324,5e-324]),[.5,.5])


def test_linkage_phase_and_recombination():
    coupling=((0,0),(1,1)); repulsion=((0,1),(1,0))
    assert linked_gametes(coupling,0)=={(0,0):.5,(1,1):.5}
    assert linked_gametes(coupling,0)!=linked_gametes(repulsion,0)
    assert linked_gametes(coupling,.5)==linked_gametes(repulsion,.5)
    assert sum(linked_cross(coupling,repulsion,.1).values())==pytest.approx(1)


def test_polygenic_dp_against_bruteforce_and_moments():
    p=[[.25,.5,.25]]*3
    result=additive_pmf(p,[1,2,3])
    oracle=np.zeros(13)
    for g in product(range(3),repeat=3): oracle[sum(a*b for a,b in zip(g,[1,2,3]))]+=math.prod(p[l][v] for l,v in enumerate(g))
    np.testing.assert_allclose(result,oracle)
    large=additive_pmf([[.25,.5,.25]]*200)
    assert len(large)==401
    np.testing.assert_allclose(np.dot(np.arange(401),large),200)
    np.testing.assert_allclose(np.dot((np.arange(401)-200)**2,large),100)
    with pytest.raises(ResourceLimitError): additive_pmf(p,[100,200,300],max_bins=10)


def test_hwe_exact_endpoints_symmetry():
    assert hwe_exact([10,0,0])==pytest.approx(1)
    assert hwe_exact([0,10,0])<.02
    assert hwe_exact([12,8,2])==pytest.approx(hwe_exact([2,8,12]))
    assert 0<=hwe_exact([25,50,25])<=1
    # Extreme tail must not be dominated by an arbitrary absolute tolerance.
    n=100
    weights=[math.comb(n,h)*math.comb(n-h,(n-h)//2)*2**h for h in range(0,n+1,2)]
    expected=float(Fraction(sum(w for w in weights if w<=2**n),sum(weights)))
    assert hwe_exact([0,n,0])==pytest.approx(expected,rel=1e-10,abs=0)


def test_population_underflow_is_explicit():
    m=InheritanceModel((2,))
    for kind in (None,m.kernel('dense'),m.kernel('csr'),m.kernel('hash')):
        with pytest.raises(FloatingPointError):m.next_generation([1.,1e-200,0.],kind)
    with pytest.raises(FloatingPointError):hwe([1.,1e-200])
    with pytest.raises(ValueError):additive_pmf([[.25,.5,.25]],max_bins=float('nan'))


def test_fixed_mates_and_factor_population_have_explicit_semantics():
    m=InheritanceModel((2,))
    op=m.fixed_mate_operator([1,0,0])
    np.testing.assert_allclose(op[:,2],[0,1,0])
    np.testing.assert_allclose(op.sum(axis=0),1)
    local=([.1,.7,.2],[.4,.2,.4])
    joint=InheritanceModel((2,2))
    factored=linkage_equilibrium_update(local,[2,2])
    np.testing.assert_allclose(joint.next_generation(np.kron(*local)),np.kron(*factored))
    correlated=np.array([.5,0,0,0,0,0,0,0,.5])
    false_factor=np.kron(*linkage_equilibrium_update([[.5,0,.5]]*2,[2,2]))
    assert np.max(abs(joint.next_generation(correlated)-false_factor))>.01
