/* Finite areas of rational curves. Hermite's reduction is solved by exact
 * coefficients; only the square-free remainder needs conic factorization. */
#include "nm.h"
#include <ctype.h>
#include <stdio.h>
#include <string.h>

static Q qi(int n) { return q_from_z(z_from_i64(n)); }
static C N(int n) { return c_const(qi(n)); }
static C mono(int v,int n) { return c_pow_int(c_letter(v),n); }
static C half_c(void) { return c_const(q_make(z_from_i64(1),z_from_i64(2))); }
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
        C t=c_mul(from_r(r_make(coeff(a,v,da),lead)),mono(v,da-db));
        *q=c_add(*q,t); a=c_sub(a,c_mul(t,b));
    }
    *r=a;
}
/* Polynomial division over the parameter field, including nonmonomial
 * coefficients such as a+b. The remainder's denominator is independent of v. */
static void divide_parameters(C a,C b,int v,R *q,R *r) {
    int db=degree(b,v); C lead=coeff(b,v,db);
    *q=r_from_c(N(0)); *r=r_from_c(a);
    while(!r_is_zero(*r) && degree(r->num,v)>=db) {
        int da=degree(r->num,v);
        R t=r_make(c_mul(coeff(r->num,v,da),mono(v,da-db)),c_mul(r->den,lead));
        *q=r_add(*q,t); *r=r_sub(*r,r_mul(t,r_from_c(b)));
        if(!r_is_zero(*r) && degree(r->num,v)>=da)
            nm_fail("internal check failed: parameter division did not lower the degree (vitiose)");
    }
}
/* Coefficient comparison over the existing exact rational matrix engine. */
static R *compare(C *columns,int n,R rhs,int v) {
    Mat m=mat_new(n,n), b=mat_new(n,1), null;
    for(int row=0;row<n;row++) {
        b.a[row]=r_make(coeff(rhs.num,v,row),rhs.den);
        for(int col=0;col<n;col++) m.a[row*n+col]=r_from_c(coeff(columns[col],v,row));
    }
    Mat sol=mat_solve(m,b,&null);
    if(null.c) nm_fail("internal check failed: conic coefficients are not unique (vitiose)");
    R *c=arena_alloc((size_t)n*sizeof(R));
    for(int i=0;i<n;i++) c[i]=sol.a[i];
    return c;
}
static int factors(C d,int v,C *f,int *mult) {
    C g=poly_gcd(d,integ_diff_poly(d,v)), s=poly_exact_div(d,g);
    int n;
    if(degree(s,v)<=1) { n=1;f[0]=s; }
    else n=elim_conic_factors(s,v,f,64);
    for(int i=0;i<n;i++) {
        /* Keep parameter leading coefficients: dividing by a+b would turn
         * this polynomial into a rational coefficient expression. */
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
    if(dd<=1) {
        Apart out={0,arena_alloc(2*sizeof(R)),arena_alloc(2*sizeof(C)),arena_alloc(2*sizeof(int))};
        if(dd==0) { out.pw[out.n]=0; out.part[out.n++]=f; }
        else {
            R q,a; divide_parameters(f.num,f.den,v,&q,&a);
            if(!r_is_zero(q)) { out.pw[out.n]=0; out.part[out.n++]=q; }
            if(!r_is_zero(a)) { out.base[out.n]=f.den; out.pw[out.n]=1; out.part[out.n++]=r_div(a,r_from_c(f.den)); }
            if(!out.n) { out.pw[out.n]=0; out.part[out.n++]=r_from_c(N(0)); }
        }
        if(!r_equal(apart_sum(out),f)) nm_fail("internal check failed: apart added back (vitiose)");
        return out;
    }
    R q,a; divide_parameters(f.num,f.den,v,&q,&a);
    Apart out={0,arena_alloc((size_t)(dd+2)*sizeof(R)),arena_alloc((size_t)(dd+2)*sizeof(C)),arena_alloc((size_t)(dd+2)*sizeof(int))};
    if(!r_is_zero(q)) { out.pw[out.n]=0; out.part[out.n++]=q; }
    if(!r_is_zero(a)) {
        C fs[64]; int mult[64],nf=factors(f.den,v,fs,mult);
        C *columns=arena_alloc((size_t)dd*sizeof(C));
        C *den=arena_alloc((size_t)dd*sizeof(C));
        int *pow=arena_alloc((size_t)dd*sizeof(int)),*fi=arena_alloc((size_t)dd*sizeof(int)),*fk=arena_alloc((size_t)dd*sizeof(int)),n=0;
        for(int i=0;i<nf;i++) for(int k=1;k<=mult[i];k++) {
            C dk=c_pow_int(fs[i],k), rest=poly_exact_div(f.den,dk);
            for(int j=0;j<degree(fs[i],v);j++) {
                columns[n]=c_mul(rest,mono(v,j)); den[n]=dk;fi[n]=i;fk[n]=k;pow[n++]=j;
            }
        }
        if(n!=dd) nm_fail("internal check failed: partial fraction dimensions (vitiose)");
        R *c=compare(columns,n,a,v);
        for(int i=0;i<n;) {
            R num=r_from_c(N(0));C d=den[i];int bi=fi[i],bk=fk[i];
            do { num=r_add(num,r_mul(c[i],r_from_c(mono(v,pow[i]))));i++; }
            while(i<n && c_equal(d,den[i]));
            if(!r_is_zero(num)) { out.base[out.n]=fs[bi]; out.pw[out.n]=bk; out.part[out.n++]=r_div(num,r_from_c(d)); }
        }
    }
    if(!out.n) { out.pw[out.n]=0; out.part[out.n++]=r_from_c(N(0)); }
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
        d=raw_add(d,(R){c_mul(t.coef.num,integ_diff_poly(t.poly,a.var)),c_mul(t.coef.den,den)});
    }
    return r_make(d.num,d.den);
}
static void term_r(Integral *a,int kind,R coef,C poly) {
    if(r_is_zero(coef)) return;
    for(int i=0;i<a->n;i++) if(a->term[i].kind==kind && c_equal(a->term[i].poly,poly)) {
        a->term[i].coef=r_add(a->term[i].coef,coef); return;
    }
    a->term[a->n++]=(AreaTerm){coef,poly,kind};
}
static void term(Integral *a,int kind,C coef,C poly) { term_r(a,kind,r_from_c(coef),poly); }

/* ---- a factor with no split over Q: Newton finds its roots and divides them out ----
 * The real roots are certified surds r_k (Sturm, elim.c); the rest is one quadratic, or for a quartic with no real
 * root two quadratics (Ferrari). Each linear factor gives a logarithm, each quadratic a logarithm and an arc. */
int (*integ_sign)(C a);

static void divide_c(C a,C b,int v,C *q,C *r) {   /* b monic in v, coefficients exact numbers and surds */
    int db=degree(b,v); *q=c_zero();
    while(!c_is_zero(a) && degree(a,v)>=db) {
        int da=degree(a,v);
        C t=c_mul(coeff(a,v,da),mono(v,da-db));
        *q=c_add(*q,t); a=c_sub(a,c_mul(t,b));
    }
    *r=a;
}
static C monic(C p,int v) { return c_mul(p,c_inv(coeff(p,v,degree(p,v)))); }

/* the numerator over G (monic, degree 1 or 2) in N/(G H): N = (b x + c) H mod G */
static void numerator(C n,C g,C h,int v,C *b,C *c) {
    C q,nm,hm; divide_c(n,g,v,&q,&nm); divide_c(h,g,v,&q,&hm);
    if(degree(g,v)==1) { *b=c_zero(); *c=c_div(coeff(nm,v,0),coeff(hm,v,0)); return; }
    C p=coeff(g,v,1),qq=coeff(g,v,0),h1=coeff(hm,v,1),h0=coeff(hm,v,0),n1=coeff(nm,v,1),n0=coeff(nm,v,0);
    /* b(h0 - h1 p) + c h1 = n1,  -b h1 qq + c h0 = n0 */
    C a11=c_sub(h0,c_mul(h1,p)),a12=h1,a21=c_neg(c_mul(h1,qq)),a22=h0;
    C det=c_sub(c_mul(a11,a22),c_mul(a12,a21));
    *b=c_div(c_sub(c_mul(n1,a22),c_mul(a12,n0)),det);
    *c=c_div(c_sub(c_mul(a11,n0),c_mul(n1,a21)),det);
}

static void real_split(Integral *out,C num,C den,int v) {
    if(degree(den,v)>4) nm_fail("a factor of degree %d with no split over the rationals: later",degree(den,v));
    if(!integ_sign) nm_fail("internal: no sign for real numbers");
    C lead=coeff(den,v,degree(den,v)); num=c_div(num,lead); den=monic(den,v);
    C g[64]; int ng=0, dd=degree(den,v);
    C quad[2]; int nq=0;
    if(dd==3) {                                  /* one certified real root, then the quadratic left */
        C rr[4]; int nr=elim_real_roots(den,v,rr,4);
        if(nr<1) nm_fail("internal: a cubic without a real root (vitiose)");
        g[ng++]=c_sub(c_letter(v),rr[0]);
        C q,r; divide_c(den,g[0],v,&q,&r);
        if(!c_is_zero(r)) nm_fail("internal check failed: a certified root does not divide (vitiose)");
        quad[nq++]=q;
    } else if(dd==4) {                           /* Ferrari: x^4 + a x^3 + b x^2 + c x + d as two real quadratics */
        Q a,b,c,d;
        if(!c_const_value(coeff(den,v,3),&a)||!c_const_value(coeff(den,v,2),&b)||!c_const_value(coeff(den,v,1),&c)||!c_const_value(coeff(den,v,0),&d))
            nm_fail("internal: Ferrari needs rational coefficients");
        /* (x^2 + a x/2 + y/2)^2 - ((a^2/4 - b + y) x^2 + (a y/2 - c) x + y^2/4 - d): a square when
         * (a y/2 - c)^2 = 4 (a^2/4 - b + y)(y^2/4 - d), a cubic in y */
        int Y=letter_index("#ferrari",8);
        C y=c_letter(Y),A=c_const(a),B=c_const(b),Cc=c_const(c),Dd=c_const(d),half=half_c(),quarter=c_const(q_make(z_from_i64(1),z_from_i64(4)));
        C k2=c_add(c_sub(c_mul(c_mul(A,A),quarter),B),y);
        C k1=c_sub(c_mul(c_mul(A,y),half),Cc), k0=c_sub(c_mul(c_mul(y,y),quarter),Dd);
        C cubic=c_sub(c_mul(k1,k1),c_mul(c_mul(N(4),k2),k0));
        C ys[8]; int ny=elim_real_roots(cubic,Y,ys,8), done=0;
        for(int i=0;i<ny && !done;i++) {
            C K2=integ_subst_poly(k2,Y,ys[i]), K1=integ_subst_poly(k1,Y,ys[i]), K0=integ_subst_poly(k0,Y,ys[i]), lin;
            C base=c_add(c_add(mono(v,2),c_mul(c_mul(A,half),c_letter(v))),c_mul(ys[i],half));
            if(c_is_zero(K2) && c_is_zero(K1) && integ_sign(K0)>0) lin=c_radical_c(K0,2);   /* x^4 - 2: (x^2)^2 - 2 */
            else {
                if(integ_sign(K2)<=0) continue;
                C al=c_radical_c(K2,2), be=c_div(K1,c_mul(N(2),al));
                lin=c_add(c_mul(al,c_letter(v)),be);
            }
            quad[nq++]=c_add(base,lin); quad[nq++]=c_sub(base,lin); done=1;
        }
        if(!done) nm_fail("internal: no real resolvent root for Ferrari");
        C prod=c_mul(quad[0],quad[1]);
        if(!c_equal(prod,den)) nm_fail("internal check failed: Ferrari's two quadratics (vitiose)");
    } else nm_fail("a factor of degree %d with no split over the rationals: later",dd);
    for(int i=0;i<nq;i++) {                      /* a quadratic with real roots splits by the formula */
        C p=coeff(quad[i],v,1),qq=coeff(quad[i],v,0), disc=c_sub(c_mul(p,p),c_mul(N(4),qq));
        if(integ_sign(disc)>0) {
            C sq=c_radical_c(disc,2);
            g[ng++]=c_add(c_letter(v),c_mul(c_add(p,sq),half_c()));
            g[ng++]=c_add(c_letter(v),c_mul(c_sub(p,sq),half_c()));
        } else g[ng++]=quad[i];
    }
    /* partial fractions over the real factors, put back exactly */
    C back=c_zero();
    C *bs=arena_alloc(64*sizeof(C)),*cs=arena_alloc(64*sizeof(C));
    for(int i=0;i<ng;i++) {
        C h,rm; divide_c(den,g[i],v,&h,&rm);
        if(!c_is_zero(rm)) nm_fail("internal check failed: a factor does not divide (vitiose)");
        numerator(num,g[i],h,v,&bs[i],&cs[i]);
        back=c_add(back,c_mul(c_add(c_mul(bs[i],c_letter(v)),cs[i]),h));
    }
    if(!c_equal(back,num)) nm_fail("internal check failed: the partial fractions over the real roots do not add back (vitiose)");
    for(int i=0;i<ng;i++) {
        if(degree(g[i],v)==1) { term(out,AREA_LOG,cs[i],g[i]); continue; }
        C p=coeff(g[i],v,1),qq=coeff(g[i],v,0);
        C w2=c_sub(c_mul(N(4),qq),c_mul(p,p));
        if(integ_sign(w2)<=0) nm_fail("internal check failed: a quadratic factor with real roots (vitiose)");
        C w=c_radical_c(w2,2), lam=c_mul(bs[i],half_c()), K=c_div(c_mul(N(2),c_sub(cs[i],c_mul(lam,p))),w);
        if(!c_equal(c_mul(w,w),w2) || !c_equal(c_add(c_div(c_mul(K,w),N(2)),c_mul(lam,p)),cs[i]))
            nm_fail("internal check failed: an arc over a real quadratic (vitiose)");
        if(!c_is_zero(lam)) term(out,AREA_LOG,lam,g[i]);
        if(!c_is_zero(K)) term(out,AREA_ATAN,K,c_div(c_add(c_mul(N(2),c_letter(v)),p),w));
    }
    out->root=out->root;                           /* (the derivative check over surds is the put-back above) */
}
Integral integ_rational(R f,int v) {
    int dd=degree(f.den,v);
    if(dd>64) nm_fail("rational integration degree limit is 64");
    Integral out={v,0,r_from_c(N(0)),arena_alloc((size_t)(2*dd+2)*sizeof(AreaTerm)),0,0,{0,0},{0,0}};
    if(dd==0) {
        out.rational=r_make(fluent_poly(f.num,v),f.den);
        if(!r_equal(integ_derivative(out),f)) nm_fail("internal check failed: integral differentiated back (vitiose)");
        return out;
    }
    if(dd==1) {
        R q,a; divide_parameters(f.num,f.den,v,&q,&a);
        out.rational=r_make(fluent_poly(q.num,v),q.den);
        term_r(&out,AREA_LOG,r_div(a,r_from_c(coeff(f.den,v,1))),f.den);
        if(!r_equal(integ_derivative(out),f)) nm_fail("internal check failed: integral differentiated back (vitiose)");
        return out;
    }
    R q,a; divide_parameters(f.num,f.den,v,&q,&a);
    out.rational=r_make(fluent_poly(q.num,v),q.den);
    int split=0;
    if(!r_is_zero(a)) {
        C g=poly_gcd(f.den,integ_diff_poly(f.den,v)),s=poly_exact_div(f.den,g);
        int dg=degree(g,v),ds=degree(s,v);
        R residual=r_div(a,r_from_c(g));
        if(dg>0) {
            /* A = S B' - (S G'/G) B + G C; integral(A/D)=B/G+integral(C/S). */
            C h=poly_exact_div(c_mul(s,integ_diff_poly(g,v)),g);
            C *cols=arena_alloc((size_t)dd*sizeof(C));
            for(int i=0;i<dg;i++) cols[i]=c_sub(c_mul(s,integ_diff_poly(mono(v,i),v)),c_mul(h,mono(v,i)));
            for(int i=0;i<ds;i++) cols[dg+i]=c_mul(g,mono(v,i));
            R *solution=compare(cols,dd,a,v),b=r_from_c(N(0)); residual=r_from_c(N(0));
            for(int i=0;i<dg;i++) b=r_add(b,r_mul(solution[i],r_from_c(mono(v,i))));
            for(int i=0;i<ds;i++) residual=r_add(residual,r_mul(solution[dg+i],r_from_c(mono(v,i))));
            out.rational=r_add(out.rational,r_div(b,r_from_c(g)));
        }
        if(!r_is_zero(residual)) {
            Apart ap=integ_apart(r_div(residual,r_from_c(s)),v);
            (void)0;
            for(int i=0;i<ap.n;i++) {
                R r=ap.part[i]; int d=degree(r.den,v);
                if(d==0) { out.rational=r_add(out.rational,r_from_c(fluent_poly(from_r(r),v)));continue; }
                C lead=coeff(r.den,v,d);
                if(d==1) { term_r(&out,AREA_LOG,r_make(r.num,lead),r.den);continue; }
                if(d>2) {
                    for(int l=0;l<letter_count();l++)
                        if(l!=v && !letter_is_surd(l) && !letter_is_named(l) && (c_uses(r.num,l)||c_uses(r.den,l)))
                            nm_fail("letters in a factor of degree %d: later",d);
                    real_split(&out,r.num,r.den,v); split=1; continue;
                }
                for(int l=0;l<letter_count();l++)
                    if(l!=v && !letter_is_surd(l) && !letter_is_named(l) && (c_uses(r.num,l)||c_uses(r.den,l)))
                        nm_fail("letters in a quadratic conic part: sign unknown; later");
                C den=c_div(r.den,lead),num=c_div(r.num,lead);
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
    if(!split && !r_equal(integ_derivative(out),f)) nm_fail("internal check failed: integral differentiated back (vitiose)");
    return out;
}
C integ_value(Integral a,C x,C (*conic)(int,C)) {
    if(c_has_plain(x)) nm_fail("a conic area takes a number");
    C num=integ_subst_poly(a.rational.num,a.var,x),den=integ_subst_poly(a.rational.den,a.var,x);
    C y=c_div(num,den);
    for(int i=0;i<a.n;i++) {
        AreaTerm t=a.term[i];
        y=c_add(y,c_mul(c_div(t.coef.num,t.coef.den),conic(t.kind,integ_subst_poly(t.poly,a.var,x))));
    }
    return y;
}
Integral integ_persist(Integral a) {
    a.rational=r_persist(a.rational); AreaTerm *t=perm_alloc((size_t)a.n*sizeof(AreaTerm));
    for(int i=0;i<a.n;i++) { t[i]=a.term[i];t[i].coef=r_persist(t[i].coef);t[i].poly=c_persist(t[i].poly); }
    a.term=t;if(a.euler) { a.ek=c_persist(a.ek); a.ec=c_persist(a.ec); }
    return a;
}
Apart apart_persist(Apart a) {
    R *p=perm_alloc((size_t)a.n*sizeof(R));C *b=perm_alloc((size_t)a.n*sizeof(C)+1);int *w=perm_alloc((size_t)a.n*sizeof(int)+1);
    for(int i=0;i<a.n;i++) { p[i]=r_persist(a.part[i]); w[i]=a.pw[i]; b[i]=a.pw[i]?c_persist(a.base[i]):c_zero(); }
    a.part=p;a.base=b;a.pw=w;return a;
}
/* Printing in the notation of the series: a sign between terms, the coefficient's denominator last
 * (log|x - 1|/2 - atan(x)/2), and every printed form reads back as the same value. */
static int lead_negative(C a) { return a.nt && q_sign(a.t[c_leading_index(a)].k) < 0; }
static int needs_paren(const char *s) { return strchr(s, ' ') != NULL || strchr(s, '/') != NULL; }
static char *cat(char *a, const char *b) {
    size_t n = strlen(a), m = strlen(b); char *s = arena_alloc(n + m + 1);
    memcpy(s, a, n); memcpy(s + n, b, m + 1); return s;
}
static R clear_fractions(R k) {                          /* (1/2)/(x+1) as 1/(2(x + 1)): whole coefficients on top */
    Z l = z_from_i64(1);
    for (int i = 0; i < k.num.nt; i++) { Z g = z_gcd(l, k.num.t[i].k.den), q, r; z_divmod(z_mul(l, k.num.t[i].k.den), g, &q, &r); l = q; }
    if (z_is_one(l)) return k;
    C m = c_const(q_from_z(l));
    k.num = c_mul(k.num, m); k.den = c_mul(k.den, m);
    return k;
}
static char *add_term(char *s, int neg, const char *body) {   /* body is printed without its sign */
    if (!*s) return neg ? cat("-", body) : cat("", body);
    return cat(cat(s, neg ? " - " : " + "), body);
}
/* k * body, body a name such as log|p| or atan(p); NULL body: k alone */
static char *scaled(R k, const char *body, int *neg) {
    k = clear_fractions(k);
    *neg = lead_negative(k.num);
    C num = *neg ? c_neg(k.num) : k.num;
    char *n = c_to_str(num), *d = c_to_str(k.den), *s = "";
    int one = !strcmp(n, "1"), done = !strcmp(d, "1");
    if (!body) s = needs_paren(n) && !done ? cat(cat("(", n), ")") : n;
    else if (one) s = cat("", body);
    else {
        s = needs_paren(n) ? cat(cat("(", n), ")") : n;
        if (isalpha((unsigned char)s[strlen(s) - 1]) || s[strlen(s) - 1] == '_') s = cat(s, "*");
        s = cat(s, body);
    }
    if (!done) s = cat(cat(s, "/"), needs_paren(d) || !isdigit((unsigned char)d[0]) || strspn(d, "0123456789") != strlen(d) ? cat(cat("(", d), ")") : d);
    return s;
}
static char *rat_term(char *s, R r) {
    int neg = lead_negative(r.num);
    if (neg) r.num = c_neg(r.num);
    return add_term(s, neg, r_to_str(r));
}
char *apart_to_str(Apart a) {
    char *s = "";
    for (int i = 0; i < a.n; i++) {
        if (!a.pw[i]) { s = rat_term(s, a.part[i]); continue; }
        C b = a.base[i]; int k = a.pw[i];
        R top = r_mul(a.part[i], r_from_c(c_pow_int(b, k)));      /* part = top / b^k, top free of the letter */
        if (lead_negative(b)) { b = c_neg(b); if (k & 1) top.num = c_neg(top.num); }
        top = clear_fractions(top);
        int neg = lead_negative(top.num);
        if (neg) top.num = c_neg(top.num);
        char *n = c_to_str(top.num), *dk = c_to_str(top.den), *bs = c_to_str(b);
        char *base = b.nt > 1 ? cat(cat("(", bs), ")") : bs;
        if (k > 1) { char e[16]; snprintf(e, sizeof e, "^%d", k); base = cat(b.nt > 1 ? base : cat(cat("(", bs), ")"), e); if (b.nt == 1 && !strchr(bs, '^') && strspn(bs, "abcdefghijklmnopqrstuvwxyz") == strlen(bs) && strlen(bs) == 1) base = cat(bs, e); }
        char *den = strcmp(dk, "1") ? cat(needs_paren(dk) ? cat(cat("(", dk), ")") : dk, base) : base;
        char *body = cat(cat(needs_paren(n) ? cat(cat("(", n), ")") : n, "/"), strcmp(dk, "1") || (b.nt > 1 && k == 1) ? (strcmp(dk, "1") ? cat(cat("(", den), ")") : den) : den);
        s = add_term(s, neg, body);
    }
    return *s ? s : "0";
}
char *integ_to_str(Integral a) {
    char *s = r_is_zero(a.rational) ? "" : rat_term("", a.rational);
    for (int i = 0; i < a.n; i++) {
        AreaTerm t = a.term[i]; if (r_is_zero(t.coef)) continue;
        char *p = c_to_str(t.poly), *body = arena_alloc(strlen(p) + 8);
        size_t pl = strlen(p);                       /* log|((s - 1)/x)|: one pair of brackets is enough */
        if (pl > 2 && p[0] == '(' && p[pl - 1] == ')') {
            int dep = 0, outer = 1;
            for (size_t i = 0; i < pl - 1 && outer; i++) { dep += p[i] == '(' ? 1 : p[i] == ')' ? -1 : 0; if (dep == 0) outer = 0; }
            if (outer) { p = arena_alloc(pl); memcpy(p, c_to_str(t.poly) + 1, pl - 2); p[pl - 2] = 0; }
        }
        sprintf(body, t.kind == AREA_LOG ? "log|%s|" : t.kind == AREA_ASIN ? "asin(%s)" : "atan(%s)", p);
        int neg; char *b = scaled(t.coef, body, &neg);
        s = add_term(s, neg, b);
    }
    return *s ? s : "0";
}
