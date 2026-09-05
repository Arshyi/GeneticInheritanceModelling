"""Independent integer-combinatorial verification of the exact HWE test."""
from fractions import Fraction
import math
import pytest
from genetics.extensions import hwe_exact


def integer_oracle(aa,het,bb):
    n=aa+het+bb;copies=het+2*bb
    weights=[];observed=None
    for h in range(copies%2,min(copies,2*n-copies)+1,2):
        b=(copies-h)//2
        # Integer multinomial count times the two orientations of each heterozygote.
        weight=math.comb(n,h)*math.comb(n-h,b)*2**h
        weights.append(weight)
        if h==het:observed=weight
    return float(Fraction(sum(w for w in weights if w<=observed),sum(weights)))


def test_every_genotype_table_through_sample_size_30():
    tested=0;largest_error=0
    for n in range(1,31):
        for aa in range(n+1):
            for het in range(n-aa+1):
                bb=n-aa-het
                largest_error=max(largest_error,abs(hwe_exact([aa,het,bb])-integer_oracle(aa,het,bb)))
                tested+=1
    assert tested==5455
    assert largest_error<7e-15


@pytest.mark.parametrize('counts',[[0,100,0],[0,1000,0],[500,0,500]])
def test_extreme_tails_against_integer_oracle(counts):
    assert hwe_exact(counts)==pytest.approx(integer_oracle(*counts),rel=1e-10,abs=0)
