/* Finite areas of rational curves. Hermite's reduction is solved by exact
 * coefficients; only the square-free remainder needs conic factorization. */
#include "nm.h"
#include <stdio.h>
#include <string.h>

static Q qi(int n) { return q_from_z(z_from_i64(n)); }
static C N(int n) { return c_const(qi(n)); }
static C mono(int v,int n) { return c_pow_int(c_letter(v),n); }
static int degree(C p,int v) { return up_from(p,v).deg; }
static C coeff(C p,int v,int n) { return c_coeff_of(p,v,qi(n)); }
static C from_r(R r) {
    if (!c_is_monomial(r.den)) nm_fail("nonpolynomial coefficient in conic reduction: later");
    return c_div(r.num,r.den);
}
C integ_diff_poly(C p,int v) {
    UP a=up_from(p,v); C d=c_zero();
    for(int i=1;i<=a.deg;i++) d=c_add(d,c_mul(c_scale(a.c[i],qi(i)),mono(v,i-1)));
    return d;
}
C integ_subst_poly(C p,int v,C x) {
    UP a=up_from(p,v); C y=c_zero();
    for(int i=a.deg;i>=0;i--) y=c_add(c_mul(y,x),a.c[i]);
    return y;
}
static C fluent_poly(C p,int v) {
    UP a=up_from(p,v); C y=c_zero();
    for(int i=0;i<=a.deg;i++) y=c_add(y,c_mul(c_scale(a.c[i],q_div(qi(1),qi(i+1))),mono(v,i+1)));
    return y;
}
static void divide(C a,C b,int v,C *q,C *r) {
    int db=degree(b,v); C lead=coeff(b,v,db); *q=c_zero();
    if(db<0) nm_fail("division by zero");
    while(degree(a,v)>=db && !c_is_zero(a)) {
        int da=degree(a,v);
        C t=c_mul(c_div(coeff(a,v,da),lead),mono(v,da-db));
        *q=c_add(*q,t); a=c_sub(a,c_mul(t,b));
    }
    *r=a;
}
/* Coefficient comparison over the existing exact rational matrix engine. */
static C *compare(C *columns,int n,C rhs,int v) {
    Mat m=mat_new(n,n), b=mat_new(n,1), null;
    for(int row=0;row<n;row++) {
        b.a[row]=r_from_c(coeff(rhs,v,row));
        for(int col=0;col<n;col++) m.a[row*n+col]=r_from_c(coeff(columns[col],v,row));
    }
    Mat sol=mat_solve(m,b,&null);
    if(null.c) nm_fail("internal check failed: conic coefficients are not unique (vitiose)");
    C *c=arena_alloc((size_t)n*sizeof(C));
    for(int i=0;i<n;i++) c[i]=from_r(sol.a[i]);
    return c;
}
static int factors(C d,int v,C *f,int *mult) {
    C g=poly_gcd(d,integ_diff_poly(d,v)), s=poly_exact_div(d,g);
    int n;
    if(degree(s,v)<=1) { n=1;f[0]=s; }
    else n=elim_conic_factors(s,v,f,64);
    for(int i=0;i<n;i++) {
        /* Monic in the integration variable, irrespective of letter order. */
        f[i]=c_div(f[i],coeff(f[i],v,degree(f[i],v)));
        mult[i]=0;
        C q,r;
        for(;;) {
            divide(d,f[i],v,&q,&r);
            if(!c_is_zero(r)) break;
            mult[i]++; d=q;
            if(degree(d,v)<=0) break;
        }
    }
    if(degree(d,v)>0) nm_fail("internal check failed: factor product (vitiose)");
    return n;
}
Apart integ_apart(R f,int v) {
    int dd=degree(f.den,v);
    if(dd>64) nm_fail("rational integration degree limit is 64");
    C q,a; divide(f.num,f.den,v,&q,&a);
    Apart out={0,arena_alloc((size_t)(dd+2)*sizeof(R))};
    if(!c_is_zero(q)) out.part[out.n++]=r_from_c(q);
    if(!c_is_zero(a)) {
        C fs[64]; int mult[64],nf=factors(f.den,v,fs,mult);
        C *columns=arena_alloc((size_t)dd*sizeof(C));
        C *den=arena_alloc((size_t)dd*sizeof(C));
        int *pow=arena_alloc((size_t)dd*sizeof(int)),n=0;
        for(int i=0;i<nf;i++) for(int k=1;k<=mult[i];k++) {
            C dk=c_pow_int(fs[i],k), rest=poly_exact_div(f.den,dk);
            for(int j=0;j<degree(fs[i],v);j++) {
                columns[n]=c_mul(rest,mono(v,j)); den[n]=dk;pow[n++]=j;
            }
        }
        if(n!=dd) nm_fail("internal check failed: partial fraction dimensions (vitiose)");
        C *c=compare(columns,n,a,v);
        for(int i=0;i<n;) {
            C num=c_zero(),d=den[i];
            do { num=c_add(num,c_mul(c[i],mono(v,pow[i])));i++; }
            while(i<n && c_equal(d,den[i]));
            if(!c_is_zero(num)) out.part[out.n++]=r_make(num,d);
        }
    }
    if(!out.n) out.part[out.n++]=r_from_c(N(0));
    if(!r_equal(apart_sum(out),f)) nm_fail("internal check failed: apart added back (vitiose)");
    return out;
}
R apart_sum(Apart a) {
    R s=r_from_c(N(0)); for(int i=0;i<a.n;i++) s=r_add(s,a.part[i]); return s;
}
/* Raw fractions allow conjugate surd terms to cancel before the Q-rational
 * normalizer sees their common denominator. These are private, unreduced. */
static R raw_add(R a,R b) {
    return (R){c_add(c_mul(a.num,b.den),c_mul(b.num,a.den)),c_mul(a.den,b.den)};
}
R integ_derivative(Integral a) {
    R r=a.rational;
    R d={c_sub(c_mul(integ_diff_poly(r.num,a.var),r.den),c_mul(r.num,integ_diff_poly(r.den,a.var))),c_pow_int(r.den,2)};
    for(int i=0;i<a.n;i++) {
        AreaTerm t=a.term[i];
        C den=t.kind==AREA_LOG?t.poly:c_add(N(1),c_pow_int(t.poly,2));
        d=raw_add(d,(R){c_mul(t.coef,integ_diff_poly(t.poly,a.var)),den});
    }
    return r_make(d.num,d.den);
}
static void term(Integral *a,int kind,C coef,C poly) {
    if(c_is_zero(coef)) return;
    for(int i=0;i<a->n;i++) if(a->term[i].kind==kind && c_equal(a->term[i].poly,poly)) {
        a->term[i].coef=c_add(a->term[i].coef,coef); return;
    }
    a->term[a->n++]=(AreaTerm){coef,poly,kind};
}
Integral integ_rational(R f,int v) {
    int dd=degree(f.den,v);
    if(dd>64) nm_fail("rational integration degree limit is 64");
    Integral out={v,0,r_from_c(N(0)),arena_alloc((size_t)(2*dd+2)*sizeof(AreaTerm))};
    C q,a; divide(f.num,f.den,v,&q,&a);
    out.rational=r_from_c(fluent_poly(q,v));
    if(!c_is_zero(a)) {
        C g=poly_gcd(f.den,integ_diff_poly(f.den,v)),s=poly_exact_div(f.den,g);
        int dg=degree(g,v),ds=degree(s,v);
        C residual=a;
        if(dg>0) {
            /* A = S B' - (S G'/G) B + G C; integral(A/D)=B/G+integral(C/S). */
            C h=poly_exact_div(c_mul(s,integ_diff_poly(g,v)),g);
            C *cols=arena_alloc((size_t)dd*sizeof(C));
            for(int i=0;i<dg;i++) cols[i]=c_sub(c_mul(s,integ_diff_poly(mono(v,i),v)),c_mul(h,mono(v,i)));
            for(int i=0;i<ds;i++) cols[dg+i]=c_mul(g,mono(v,i));
            C *solution=compare(cols,dd,a,v),b=c_zero(); residual=c_zero();
            for(int i=0;i<dg;i++) b=c_add(b,c_mul(solution[i],mono(v,i)));
            for(int i=0;i<ds;i++) residual=c_add(residual,c_mul(solution[dg+i],mono(v,i)));
            out.rational=r_add(out.rational,r_make(b,g));
        }
        if(!c_is_zero(residual)) {
            Apart ap=integ_apart(r_make(residual,s),v);
            for(int i=0;i<ap.n;i++) {
                R r=ap.part[i]; int d=degree(r.den,v);
                if(d==0) { out.rational=r_add(out.rational,r_from_c(fluent_poly(from_r(r),v)));continue; }
                C lead=coeff(r.den,v,d),den=c_div(r.den,lead),num=c_div(r.num,lead);
                if(d==1) { term(&out,AREA_LOG,num,den);continue; }
                if(d!=2) nm_fail("remaining degree-%d factor; later: Rothstein-Trager",d);
                Q b,c;
                if(!c_const_value(coeff(den,v,1),&b) || !c_const_value(coeff(den,v,0),&c) || c_has_plain(coeff(num,v,0)) || c_has_plain(coeff(num,v,1)))
                    nm_fail("letters in a quadratic conic part: sign unknown; later");
                C lambda=c_scale(coeff(num,v,1),q_div(qi(1),qi(2)));
                C rem=c_sub(coeff(num,v,0),c_scale(lambda,b));
                term(&out,AREA_LOG,lambda,den);
                Q delta=q_sub(q_mul(b,b),q_mul(qi(4),c));
                if(q_sign(delta)<0) {
                    C root=c_radical_q(q_neg(delta),2);
                    C arg=c_div(c_add(c_scale(c_letter(v),qi(2)),c_const(b)),root);
                    term(&out,AREA_ATAN,c_div(c_scale(rem,qi(2)),root),arg);
                } else if(q_sign(delta)>0) {
                    C root=c_radical_q(delta,2), mid=c_add(c_letter(v),c_const(q_div(b,qi(2))));
                    C half=c_scale(root,q_div(qi(1),qi(2))),k=c_div(rem,root);
                    term(&out,AREA_LOG,k,c_sub(mid,half));
                    term(&out,AREA_LOG,c_neg(k),c_add(mid,half));
                } else nm_fail("internal check failed: repeated square-free factor (vitiose)");
            }
        }
    }
    if(!r_equal(integ_derivative(out),f)) nm_fail("internal check failed: integral differentiated back (vitiose)");
    return out;
}
C integ_value(Integral a,C x,C (*conic)(int,C)) {
    if(c_has_plain(x)) nm_fail("a conic area takes a number");
    C num=integ_subst_poly(a.rational.num,a.var,x),den=integ_subst_poly(a.rational.den,a.var,x);
    C y=c_div(num,den);
    for(int i=0;i<a.n;i++) {
        AreaTerm t=a.term[i];
        y=c_add(y,c_mul(t.coef,conic(t.kind,integ_subst_poly(t.poly,a.var,x))));
    }
    return y;
}
Integral integ_persist(Integral a) {
    a.rational=r_persist(a.rational); AreaTerm *t=perm_alloc((size_t)a.n*sizeof(AreaTerm));
    for(int i=0;i<a.n;i++) { t[i]=a.term[i];t[i].coef=c_persist(t[i].coef);t[i].poly=c_persist(t[i].poly); }
    a.term=t;return a;
}
Apart apart_persist(Apart a) {
    R *p=perm_alloc((size_t)a.n*sizeof(R));for(int i=0;i<a.n;i++) p[i]=r_persist(a.part[i]);a.part=p;return a;
}
/* Size the strings from the actual terms; no fixed display buffer. */
static char *join(char *a,const char *b) {
    size_t n=strlen(a),m=strlen(b);char *s=arena_alloc(n+m+4);
    snprintf(s,n+m+4,"%s%s%s",a,n?" + ":"",b);return s;
}
char *apart_to_str(Apart a) {
    char *s="";for(int i=0;i<a.n;i++) s=join(s,r_to_str(a.part[i]));return s;
}
char *integ_to_str(Integral a) {
    char *s=r_is_zero(a.rational)?"":r_to_str(a.rational);
    for(int i=0;i<a.n;i++) {
        AreaTerm t=a.term[i];if(c_is_zero(t.coef)) continue;
        char *k=c_to_str(t.coef),*p=c_to_str(t.poly);
        size_t n=strlen(k)+strlen(p)+24;char *b=arena_alloc(n);
        int one=c_equal(t.coef,N(1));
        snprintf(b,n,"%s%s%s%s%s%s",one?"":"(",one?"":k,one?"":")*",t.kind==AREA_LOG?"log|":"atan(",p,t.kind==AREA_LOG?"|":")");
        s=join(s,b);
    }
    return *s?s:"0";
}
