/* Rational matrices: put-backs, exact exceptional loci and an independent Q route. */
#include "nm.h"
#include <setjmp.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

extern jmp_buf nm_on_error;
extern char nm_error_msg[512];
static int checks, fail, putbacks, lowest, numeric;
static Q qi(int64_t n) { return q_from_z(z_from_i64(n)); }
static R N(int n) { return r_from_c(c_const(qi(n))); }
static R L(const char *s) { return r_from_c(c_letter(letter_index(s,strlen(s)))); }
static void check(int ok, const char *s) { checks++; if (!ok) { printf("FAIL: %s\n",s); fail++; } }
static void reduced(Mat m) {
    for (int i=0;i<m.r*m.c;i++) {
        lowest++; check(c_equal(poly_gcd(m.a[i].num,m.a[i].den),c_const(qi(1))),"matrix entry in lowest terms");
    }
}
static void inverse_check(Mat a) {
    Mat inv=mat_inverse(a), id=mat_identity(a.r);
    putbacks++; check(mat_equal(mat_mul(a,inv),id),"A times inverse");
    putbacks++; check(mat_equal(mat_mul(inv,a),id),"inverse times A");
    reduced(inv);
}
static Q at(C p,Q *v) {
    Q sum=qi(0);
    for(int t=0;t<p.nt;t++) {
        Q m=p.t[t].k;
        for(int l=0;l<letter_count();l++) if(q_sign(p.t[t].e[l])) {
            int64_t e=0; z_fits_i64(p.t[t].e[l].num,&e); m=q_mul(m,q_pow(v[l],e));
        }
        sum=q_add(sum,m);
    }
    return sum;
}
static uint32_t seed=65537;
static Q sample(void) {
    seed=1664525u*seed+1013904223u; int n=(int)(seed%13)-6;
    seed=1664525u*seed+1013904223u;
    return q_make(z_from_i64(n),z_from_i64(seed%4+1));
}
int main(void) {
    clock_t start=clock();
    if(setjmp(nm_on_error)) { printf("error: %s\n",nm_error_msg); return 1; }
    R a=L("a"),b=L("b"),c=L("c"),d=L("d");
    int ia=letter_index("a",1),ib=letter_index("b",1),ic=letter_index("c",1),id=letter_index("d",1);
    Mat m=mat_new(2,2); m.a[0]=a; m.a[1]=b; m.a[2]=c; m.a[3]=d;
    inverse_check(m);
    check(!strcmp(mat_verdict(m,mat_new(0,0),mat_rank(m)),
        "[exact, for general a, b, c, d; rank drops where ad - bc = 0]"),"full rank locus excludes spurious pivot a=0");
    Mat sym=mat_new(2,2); sym.a[0]=N(1);sym.a[1]=a;sym.a[2]=a;sym.a[3]=N(1);
    inverse_check(sym);
    Mat rhs=mat_new(2,1);rhs.a[0]=N(1);
    Mat ns,x=mat_solve(sym,rhs,&ns);
    putbacks++;check(mat_equal(mat_mul(sym,x),rhs),"linear solution put back"); reduced(x);
    Mat row=mat_new(1,2);row.a[0]=a;row.a[1]=b;
    ns=mat_nullspace(row);putbacks++;check(mat_equal(mat_mul(row,ns),mat_new(1,1)),"null basis put back");reduced(ns);
    check(!strcmp(mat_verdict(row,ns,1),
        "[exact, for general a, b; rank drops where a = 0 and b = 0; displayed formula undefined where a = 0]"),"rank locus and basis pole are distinct");
    row.a[1]=N(1);
    check(strstr(mat_verdict(row,mat_nullspace(row),1),"rank never drops on the input domain")!=NULL,"constant minor prevents rank drop");
    Mat rational=mat_new(2,2);rational.a[0]=r_div(N(1),a);rational.a[3]=N(1);
    inverse_check(rational);
    check(!strcmp(mat_verdict(rational,mat_new(0,0),2),
        "[exact, for general a; input undefined where a = 0; rank never drops on the input domain]"),"domain poles are not rank drops");
    Mat always=mat_new(2,2);always.a[0]=a;always.a[1]=a;always.a[2]=N(1);always.a[3]=N(1);
    jmp_buf saved;memcpy(saved,nm_on_error,sizeof saved);
    if(!setjmp(nm_on_error)) { (void)mat_inverse(always);check(0,"singular inverse refused"); }
    else check(strstr(nm_error_msg,"singular")!=NULL,"singular refusal reason");
    memcpy(nm_on_error,saved,sizeof saved);
    Mat rect=mat_new(2,3);rect.a[0]=a;rect.a[1]=b;rect.a[4]=c;rect.a[5]=d;
    check(!strcmp(mat_verdict(rect,mat_new(0,0),2),
        "[exact, for general a, b, c, d; rank drops where ac = 0 and ad = 0 and bd = 0]"),"all rectangular maximal minors");
    Mat vand=mat_new(3,3); R letters[3]={a,b,c};
    for(int i=0;i<3;i++)for(int j=0;j<3;j++)vand.a[3*i+j]=r_pow_int(letters[i],j);
    inverse_check(vand);
    Mat inv=mat_inverse(m);
    R stored[4]; for(int i=0;i<4;i++)stored[i]=r_persist(inv.a[i]);
    /* Original 2x2 inverse computed directly in Q, not using any matrix or R operation. */
    int sets=0,skips=0;
    while(sets<64) {
        arena_reset();
        Q A=sample(),B=sample(),Cc=sample(),D=sample();
        if(sets==0) { A=qi(0);B=qi(1);Cc=qi(2);D=qi(3); }
        Q det=q_sub(q_mul(A,D),q_mul(B,Cc));
        if(!q_sign(det)){skips++;continue;}
        Q v[NM_MAXL];for(int l=0;l<NM_MAXL;l++)v[l]=qi(0);
        v[ia]=A;v[ib]=B;v[ic]=Cc;v[id]=D;
        Q want[4]={q_div(D,det),q_div(q_neg(B),det),q_div(q_neg(Cc),det),q_div(A,det)};
        Q got[4];for(int i=0;i<4;i++) {
            got[i]=q_div(at(stored[i].num,v),at(stored[i].den,v));
            numeric++;check(q_cmp(got[i],want[i])==0,"inverse against independent fractions");
        }
        Q orig[4]={A,B,Cc,D};
        for(int i=0;i<2;i++)for(int j=0;j<2;j++) {
            Q sum=qi(0);for(int k=0;k<2;k++)sum=q_add(sum,q_mul(orig[i*2+k],got[k*2+j]));
            numeric++;check(q_cmp(sum,qi(i==j))==0,"numeric original A times symbolic inverse is I");
        }
        sets++;
    }
    printf("ratmat: %d put-backs, %d lowest-terms, %d Q comparisons (%d sets, %d skipped); %.2f ms\n",
        putbacks,lowest,numeric,sets,skips,1000.0*(clock()-start)/CLOCKS_PER_SEC);
    printf("unit ratmat: %d checks, %d failed\n",checks,fail);
    return fail!=0;
}
