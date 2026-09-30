/* Newton's parallelogram: the start of a root y of F(x, y) = 0 as a power of x, then the rest.
 *
 * Mark each term x^i y^j at the point (j, i). Lay a ruler under the marks and turn it from the lowest mark in the
 * left column until it touches others; the terms on the ruler, set equal to nothing, give the first term c x^g
 * of the root ("quorum quemlibet pro primo termino Quotientis accipere liceat"). Put y = x^g (c + z) and
 * continue: while c is a simple root, each next term follows from the lowest terms; if it is a multiple root, the
 * parallelogram is applied again. Finally the root is substituted back into F. */
#include "nm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static Q q0(void) { return q_from_z(z_zero()); }
static Q qi(int64_t v) { return q_from_z(z_from_i64(v)); }

typedef struct { int j; Q i; C k; } Pt;          /* the mark of the terms k x^i y^j */

static Q gcdq_lcm_den(Q a, Q l) {               /* lcm(l, den(a)) as an integer rational */
    Z g = z_gcd(l.num, a.den), q, r;
    z_divmod(z_mul(l.num, a.den), g, &q, &r);
    return q_from_z(q);
}

/* the marks of F in letters (T, Z) */
static int marks(C F, int ti, int zi, Pt **out) {
    Pt *p = arena_alloc((size_t)(F.nt + 1) * sizeof(Pt));
    int np = 0;
    for (int a = 0; a < F.nt; a++) {
        CT t = F.t[a];
        int64_t j;
        if (!q_is_int(t.e[zi]) || !z_fits_i64(t.e[zi].num, &j) || j < 0)
            nm_fail("the equation must be a polynomial in %s (whole powers)", letter_name(zi));
        Q i = t.e[ti];
        CT rest = t; rest.e[ti] = q0(); rest.e[zi] = q0();
        C one; one.nt = 1; one.t = arena_alloc(sizeof(CT)); one.t[0] = rest;
        int found = -1;
        for (int b = 0; b < np; b++) if (p[b].j == j && q_cmp(p[b].i, i) == 0) found = b;
        if (found >= 0) p[found].k = c_add(p[found].k, one);
        else { p[np].j = (int)j; p[np].i = i; p[np].k = one; np++; }
    }
    int w = 0;
    for (int b = 0; b < np; b++) if (!c_is_zero(p[b].k)) p[w++] = p[b];
    *out = p;
    return w;
}

typedef struct { int j1, j2; Q gamma; } Edge;

/* the ruler: the lower edges of the marks, from the left column to the right */
static int ruler(Pt *p, int np, Edge *edges) {
    int maxj = 0;
    for (int a = 0; a < np; a++) if (p[a].j > maxj) maxj = p[a].j;
    Q *low = arena_alloc((size_t)(maxj + 1) * sizeof(Q));
    int *has = arena_alloc((size_t)(maxj + 1) * sizeof(int));
    memset(has, 0, (size_t)(maxj + 1) * sizeof(int));
    for (int a = 0; a < np; a++) {
        int j = p[a].j;
        if (!has[j] || q_cmp(p[a].i, low[j]) < 0) { low[j] = p[a].i; has[j] = 1; }
    }
    int *hull = arena_alloc((size_t)(maxj + 1) * sizeof(int));
    int h = 0;
    for (int j = 0; j <= maxj; j++) {
        if (!has[j]) continue;
        while (h >= 2) {                         /* keep the chain convex from below */
            int a = hull[h - 2], b = hull[h - 1];
            Q cross = q_sub(q_mul(q_sub(low[b], low[a]), qi(j - a)), q_mul(q_sub(low[j], low[a]), qi(b - a)));
            if (q_sign(cross) >= 0) h--; else break;
        }
        hull[h++] = j;
    }
    int ne = 0;
    for (int a = 0; a + 1 < h; a++) {
        Edge e;
        e.j1 = hull[a]; e.j2 = hull[a + 1];
        e.gamma = q_neg(q_div(q_sub(low[e.j2], low[e.j1]), qi(e.j2 - e.j1)));
        edges[ne++] = e;
    }
    return ne;
}

typedef struct { Q gamma; C c; int mult; int exact; char *note; } Start;

static int64_t small_int(Z a) { int64_t v; return z_fits_i64(a, &v) && v < 1000000000000LL && v > -1000000000000LL ? v : 0; }

/* the roots of the ruler's equation  sum k_j c^j = 0,  c != 0 */
static int ruler_roots(Pt *p, int np, Edge e, Start *st, int max) {
    C *k = arena_alloc((size_t)(e.j2 + 1) * sizeof(C));
    for (int j = 0; j <= e.j2; j++) k[j] = c_zero();
    /* the marks on the ruler: those at the ruler's level i + g j */
    Q level, lv = q0(); int got = 0;
    for (int a = 0; a < np; a++) if (p[a].j == e.j1) {
        Q v = q_add(p[a].i, q_mul(e.gamma, qi(p[a].j)));
        if (!got || q_cmp(v, lv) < 0) { lv = v; got = 1; }
    }
    level = lv;
    for (int a = 0; a < np; a++)
        if (q_cmp(q_add(p[a].i, q_mul(e.gamma, qi(p[a].j))), level) == 0) k[p[a].j] = c_add(k[p[a].j], p[a].k);
    /* Newton's reduction (y = v sqrt(ax)): c = v M, where M makes every term the same in the letters */
    for (int j = e.j1; j <= e.j2; j++)
        if (!c_is_zero(k[j]) && !c_is_monomial(k[j]))
            nm_fail("the terms on the ruler have coefficients that are not single terms (%s); this comes later", c_to_str(k[j]));
    CT L1 = k[e.j1].t[0], L2 = k[e.j2].t[0];
    Q eM[NM_MAXL];
    for (int l = 0; l < NM_MAXL; l++) eM[l] = q_div(q_sub(L1.e[l], L2.e[l]), qi(e.j2 - e.j1));
    Q *r = arena_alloc((size_t)(e.j2 - e.j1 + 1) * sizeof(Q));
    for (int j = e.j1; j <= e.j2; j++) {
        r[j - e.j1] = q0();
        if (c_is_zero(k[j])) continue;
        CT t = k[j].t[0];
        for (int l = 0; l < NM_MAXL; l++)
            if (q_cmp(q_add(t.e[l], q_mul(eM[l], qi(j))), q_add(L1.e[l], q_mul(eM[l], qi(e.j1)))) != 0)
                nm_fail("the terms on the ruler do not reduce to numbers alone; this comes later");
        r[j - e.j1] = t.k;
    }
    int deg = e.j2 - e.j1;
    /* whole-number coefficients */
    Q l = qi(1);
    for (int j = 0; j <= deg; j++) if (q_sign(r[j])) l = gcdq_lcm_den(r[j], l);
    Z *a = arena_alloc((size_t)(deg + 1) * sizeof(Z));
    for (int j = 0; j <= deg; j++) a[j] = q_mul(r[j], l).num;
    int ns = 0;
    /* rational roots v = p/q: p divides a[0], q divides a[deg] */
    int64_t a0 = small_int(z_abs(a[0])), an = small_int(z_abs(a[deg]));
    int64_t pd[4096], qd[4096]; int npd = 0, nqd = 0;
    for (int64_t d = 1; d * d <= a0 && npd < 4000; d++) if (a0 % d == 0) { pd[npd++] = d; if (d != a0 / d) pd[npd++] = a0 / d; }
    for (int64_t d = 1; d * d <= an && nqd < 4000; d++) if (an % d == 0) { qd[nqd++] = d; if (d != an / d) qd[nqd++] = an / d; }
    if (!a0 || !an) { npd = 1; pd[0] = 1; nqd = 1; qd[0] = 1; }
    int left = deg;
    Z *cur = a;
    for (int x = 0; x < npd && left > 0; x++)
        for (int y = 0; y < nqd && left > 0; y++)
            for (int sg = -1; sg <= 1 && left > 0; sg += 2) {
                Q v = q_make(z_from_i64(sg * pd[x]), z_from_i64(qd[y]));
                int dup = 0;
                for (int s = 0; s < ns; s++) if (q_cmp(st[s].c.t[0].k, v) == 0 && q_cmp(st[s].gamma, e.gamma) == 0) dup = 1;
                if (dup) continue;
                int m = 0;
                for (;;) {                           /* divide out (c - v) while it is a root */
                    Q acc = q0();
                    for (int j = left; j >= 0; j--) acc = q_add(q_mul(acc, v), q_from_z(cur[j]));
                    if (q_sign(acc)) break;
                    Z *nx = arena_alloc((size_t)(left + 1) * sizeof(Z));
                    Q carry = q0();
                    Q *qq = arena_alloc((size_t)(left + 1) * sizeof(Q));
                    for (int j = left; j >= 1; j--) { carry = q_add(q_mul(carry, v), q_from_z(cur[j])); qq[j - 1] = carry; }
                    Q ll = qi(1);
                    for (int j = 0; j < left; j++) if (q_sign(qq[j])) ll = gcdq_lcm_den(qq[j], ll);
                    for (int j = 0; j < left; j++) nx[j] = q_mul(qq[j], ll).num;
                    cur = nx; left--; m++;
                }
                if (m && ns < max) {
                    C c; c.nt = 1; c.t = arena_alloc(sizeof(CT));
                    c.t[0].k = v;
                    for (int l2 = 0; l2 < NM_MAXL; l2++) c.t[0].e[l2] = eM[l2];
                    st[ns].gamma = e.gamma; st[ns].c = c; st[ns].mult = m; st[ns].exact = 1; st[ns].note = NULL;
                    ns++;
                }
            }
    if (left > 0 && ns < max) {                   /* roots that are not rational numbers */
        Poly rest; rest.deg = left; rest.var = "v";
        rest.c = arena_alloc((size_t)(left + 1) * sizeof(Q));
        for (int j = 0; j <= left; j++) rest.c[j] = q_from_z(cur[j]);
        C M; M.nt = 1; M.t = arena_alloc(sizeof(CT)); M.t[0].k = qi(1);
        for (int l2 = 0; l2 < NM_MAXL; l2++) M.t[0].e[l2] = eM[l2];
        st[ns].gamma = e.gamma; st[ns].c = M; st[ns].mult = left; st[ns].exact = 0;
        st[ns].note = p_to_str(rest);
        ns++;
    }
    return ns;
}

/* F(T^q, T^p (c + z)) / T^m as a quantity in the new letters t, z */
static C substitute(C F, int Ti, int Zi, int ti, int zi, int64_t q, int64_t p, C c) {
    C out = c_zero();
    C cz = c_add(c, c_letter(zi));
    for (int a = 0; a < F.nt; a++) {
        CT t = F.t[a];
        int64_t j; z_fits_i64(t.e[Zi].num, &j);
        Q te = q_add(q_mul(t.e[Ti], qi(q)), qi(p * j));
        if (!q_is_int(te)) nm_fail("internal: a fractional power after substitution");
        CT rest = t; rest.e[Ti] = q0(); rest.e[Zi] = q0(); rest.e[ti] = te;
        C one; one.nt = 1; one.t = arena_alloc(sizeof(CT)); one.t[0] = rest;
        out = c_add(out, c_mul(one, c_pow_int(cz, j)));
    }
    Q m = q0(); int first = 1;
    for (int a = 0; a < out.nt; a++) if (first || q_cmp(out.t[a].e[ti], m) < 0) { m = out.t[a].e[ti]; first = 0; }
    C tm; tm.nt = 1; tm.t = arena_alloc(sizeof(CT)); tm.t[0].k = qi(1);
    for (int l = 0; l < NM_MAXL; l++) tm.t[0].e[l] = q0();
    tm.t[0].e[ti] = m;
    return c_div(out, tm);
}

/* F(t, Z(t)) as a series in t to n terms (whole powers of t) */
static Ser eval_series(C F, int ti, int zi, Ser Zs, int n) {
    Ser acc = s_const(c_zero(), n);
    for (int a = 0; a < F.nt; a++) {
        CT t = F.t[a];
        int64_t ex, j;
        z_fits_i64(t.e[ti].num, &ex); z_fits_i64(t.e[zi].num, &j);
        if (ex >= n) continue;
        CT rest = t; rest.e[ti] = q0(); rest.e[zi] = q0();
        C k; k.nt = 1; k.t = arena_alloc(sizeof(CT)); k.t[0] = rest;
        Ser term = s_pow_int(s_trunc(Zs, n), j);
        Ser sh = s_const(c_zero(), n);
        for (int i = 0; i + ex < n && i < term.n; i++) sh.c[i + ex] = c_mul(term.c[i], k);
        if (term.n + ex < n) sh.n = (int)(term.n + ex);
        acc = s_add(acc, sh);
    }
    return acc;
}

static char *start_str(Start s, const char *x, const char *y, Q offset, Q scale) {
    Q e = q_add(offset, q_mul(s.gamma, scale));
    char *body = s.exact ? c_term_str(s.c, x, e, 1) : NULL;
    char *out = arena_alloc(2048 + (body ? strlen(body) : 0) + (s.note ? strlen(s.note) : 0));
    if (s.exact) sprintf(out, "%s = %s%s", y, body, s.mult > 1 ? (s.mult == 2 ? "  (a double root: give the next term too)" : "  (a multiple root: give the next terms too)") : "");
    else sprintf(out, "%s = v%s with %s = 0  (not rational numbers; later)", y, c_term_str(s.c, x, e, 1), s.note);
    return out;
}

typedef struct { Q e; C c; } Term;

/* the root y of F = 0 in powers of x (letters xi, yi); start: the terms the user gave, or NULL */
char *parallelogram(C F, int xi, int yi, C *start, int have_start, int64_t order) {
    const char *x = letter_name(xi), *y = letter_name(yi);
    Term *found = arena_alloc(64 * sizeof(Term));
    int nfound = 0;
    int T = xi, Zl = yi;
    Q G = q0();                                   /* y = sum of found terms + x^G z, and T = x^(1/Qs) */
    Q Qs = qi(1);
    C S = have_start ? *start : c_zero();         /* the user's terms, in the letter T */
    C cur = F;
    for (int level = 1; level <= 12; level++) {
        Pt *p; int np = marks(cur, T, Zl, &p);
        if (np == 0) nm_fail("the equation is identically zero");
        Edge edges[64];
        int ne = ruler(p, np, edges);
        Start st[64]; int ns = 0;
        for (int a = 0; a < ne && ns < 60; a++) {
            if (level > 1 && q_sign(edges[a].gamma) <= 0) continue;   /* the rest must be smaller than what came before */
            ns += ruler_roots(p, np, edges[a], st + ns, 60 - ns);
        }
        /* which start: the user's lowest remaining term */
        int pick = -1;
        if (!c_is_zero(S)) {
            Q emin = q0(); int f = 1;
            for (int a = 0; a < S.nt; a++) if (f || q_cmp(S.t[a].e[T], emin) < 0) { emin = S.t[a].e[T]; f = 0; }
            C lead = c_coeff_of(S, T, emin);
            for (int a = 0; a < ns; a++) if (st[a].exact && q_cmp(st[a].gamma, emin) == 0 && c_equal(st[a].c, lead)) pick = a;
            if (pick < 0) {
                char *u = arena_alloc(64 + strlen(c_to_str(lead)));
                sprintf(u, "%s", c_term_str(lead, x, q_add(G, q_div(emin, Qs)), 1));
                nm_fail("%s is not a start of the root here; ask without 'starting' to see the starts", u);
            }
        }
        if (pick < 0) {                            /* list the starts, as Newton listed the four */
            size_t cap = 256;
            char **lines = arena_alloc((size_t)(ns + 1) * sizeof(char *));
            for (int a = 0; a < ns; a++) { lines[a] = start_str(st[a], x, y, G, q_div(qi(1), Qs)); cap += strlen(lines[a]) + 64; }
            char *pre = arena_alloc(cap + 64 * (size_t)nfound + 1024);
            char *o = pre;
            if (nfound) {
                o += sprintf(o, "after %s = ", y);
                for (int a = 0; a < nfound; a++) o += sprintf(o, "%s", c_term_str(found[a].c, x, found[a].e, a == 0));
                o += sprintf(o, " + ..., ");
            }
            o += sprintf(o, "%d start%s; choose one with 'starting %s = ...':", ns, ns == 1 ? "" : "s", y);
            for (int a = 0; a < ns; a++) o += sprintf(o, "\n  %s", lines[a]);
            return pre;
        }
        Start s = st[pick];
        /* record the term, and put y = T^g (c + z) with T = t^q */
        Q eg = q_add(G, q_div(s.gamma, Qs));
        found[nfound].e = eg; found[nfound].c = s.c; nfound++;
        Q qq = gcdq_lcm_den(s.gamma, qi(1));
        for (int a = 0; a < cur.nt; a++) qq = gcdq_lcm_den(cur.t[a].e[T], qq);
        int64_t q, pp;
        z_fits_i64(qq.num, &q);
        z_fits_i64(q_mul(s.gamma, qq).num, &pp);
        char nm[32];
        snprintf(nm, sizeof nm, "#t%d", level); int ti = letter_index(nm, strlen(nm));
        snprintf(nm, sizeof nm, "#z%d", level); int zi = letter_index(nm, strlen(nm));
        C next = substitute(cur, T, Zl, ti, zi, q, pp, s.c);
        /* the user's remaining terms, as z in the new letter */
        if (!c_is_zero(S)) {
            C xg; xg.nt = 1; xg.t = arena_alloc(sizeof(CT)); xg.t[0].k = qi(1);
            for (int l = 0; l < NM_MAXL; l++) xg.t[0].e[l] = q0();
            xg.t[0].e[T] = s.gamma;
            C rest = c_sub(c_div(S, xg), s.c);
            for (int a = 0; a < rest.nt; a++) { rest.t[a].e[ti] = q_mul(rest.t[a].e[T], qi(q)); rest.t[a].e[T] = q0(); }
            S = rest;
        }
        G = eg;
        Qs = q_mul(Qs, qi(q));
        if (s.mult > 1) { cur = next; T = ti; Zl = zi; continue; }
        /* a simple root: the rest term by term, z(0) = 0, each coefficient from the lowest term left */
        C one; one.nt = 1; one.t = arena_alloc(sizeof(CT)); one.t[0].k = qi(1);
        for (int l = 0; l < NM_MAXL; l++) one.t[0].e[l] = q0();
        C dz = c_zero();                         /* dF/dz at t = 0, z = 0 */
        for (int a = 0; a < next.nt; a++)
            if (q_sign(next.t[a].e[ti]) == 0 && q_cmp_one(next.t[a].e[zi]) == 0) {
                CT r = next.t[a]; r.e[zi] = q0();
                C cc; cc.nt = 1; cc.t = arena_alloc(sizeof(CT)); cc.t[0] = r;
                dz = c_add(dz, cc);
            }
        if (c_is_zero(dz)) nm_fail("internal: a simple start with no slope");
        /* how many terms of z: its terms are x^(G + k/Qs) */
        Q need = q_mul(q_sub(qi(order), G), Qs);
        Z nq, nr;
        z_divmod(need.num, need.den, &nq, &nr);
        int64_t nz; z_fits_i64(nq, &nz);
        if (q_sign(need) < 0) nz = -1;
        nz += 2;
        if (nz < 1) nz = 1;
        if (nz > 400) nm_fail("order too large for the parallelogram in this version");
        Ser Zs = s_const(c_zero(), (int)nz);
        for (int d = 1; d < nz; d++) {
            Ser Fz = eval_series(next, ti, zi, s_trunc(Zs, d + 1), d + 1);
            C Fd = Fz.n > d ? Fz.c[d] : c_zero();
            Zs.c[d] = c_neg(c_div(Fd, dz));
        }
        /* write y = found terms + x^G z(x^(1/Qs)), and substitute it back into F */
        size_t cap = 256;
        char **parts = arena_alloc((size_t)(nfound + nz + 2) * sizeof(char *));
        int np2 = 0;
        Term *all = arena_alloc((size_t)(nfound + nz + 2) * sizeof(Term));
        int na = 0;
        for (int a = 0; a < nfound; a++) all[na++] = found[a];
        for (int k = 1; k < nz; k++) if (!c_is_zero(Zs.c[k])) { all[na].e = q_add(G, q_div(qi(k), Qs)); all[na].c = Zs.c[k]; na++; }
        Q bound = q_add(G, q_div(qi(nz), Qs));
        for (int a = 0; a < na; a++) {
            if (q_cmp(all[a].e, qi(order)) > 0) continue;
            parts[np2] = c_term_str(all[a].c, x, all[a].e, np2 == 0);
            cap += strlen(parts[np2]); np2++;
        }
        Q shown = q_cmp(bound, qi(order + 1)) < 0 ? bound : qi(order + 1);
        for (int a = 0; a < na; a++) if (q_cmp(all[a].e, qi(order)) > 0 && q_cmp(all[a].e, shown) < 0) shown = all[a].e;
        /* the check by substitution, in u = x^(1/D) */
        {
            Q D = Qs;
            for (int a = 0; a < F.nt; a++) D = gcdq_lcm_den(F.t[a].e[xi], D);
            int64_t Dn; z_fits_i64(D.num, &Dn);
            int64_t vmin = 0; int f2 = 1;
            for (int a = 0; a < na; a++) {
                Q ue = q_mul(all[a].e, D); int64_t v; z_fits_i64(ue.num, &v);
                if (f2 || v < vmin) { vmin = v; f2 = 0; }
            }
            int64_t top; { Q ue = q_mul(bound, D); z_fits_i64(ue.num, &top); }
            int nW = (int)(top - vmin);
            if (nW < 1) nW = 1;
            Ser W = s_const(c_zero(), nW);
            for (int a = 0; a < na; a++) {
                Q ue = q_mul(all[a].e, D); int64_t v; z_fits_i64(ue.num, &v);
                if (v - vmin < nW) W.c[v - vmin] = c_add(W.c[v - vmin], all[a].c);
            }
            int64_t lo = 0, known = 0; int f3 = 1;
            for (int a = 0; a < F.nt; a++) {
                int64_t j; z_fits_i64(F.t[a].e[yi].num, &j);
                Q ie = q_mul(F.t[a].e[xi], D); int64_t iv; z_fits_i64(ie.num, &iv);
                int64_t e0 = iv + j * vmin;
                if (f3 || e0 < lo) lo = e0;
                if (f3 || e0 + nW < known) known = e0 + nW;
                f3 = 0;
            }
            int len = (int)(known - lo);
            if (len > 0) {
                C *sum = arena_alloc((size_t)len * sizeof(C));
                for (int i = 0; i < len; i++) sum[i] = c_zero();
                for (int a = 0; a < F.nt; a++) {
                    int64_t j; z_fits_i64(F.t[a].e[yi].num, &j);
                    Q ie = q_mul(F.t[a].e[xi], D); int64_t iv; z_fits_i64(ie.num, &iv);
                    CT r = F.t[a]; r.e[xi] = q0(); r.e[yi] = q0();
                    C k; k.nt = 1; k.t = arena_alloc(sizeof(CT)); k.t[0] = r;
                    Ser Wj = s_pow_int(W, j);
                    int64_t off = iv + j * vmin - lo;
                    for (int i = 0; i < Wj.n && off + i < len; i++) sum[off + i] = c_add(sum[off + i], c_mul(Wj.c[i], k));
                }
                for (int i = 0; i < len; i++)
                    if (!c_is_zero(sum[i])) nm_fail("internal check failed: the root does not satisfy the equation when substituted back (vitiose)");
            }
        }
        char *out = arena_alloc(cap + 128), *o = out;
        *o = 0;
        for (int a = 0; a < np2; a++) o += sprintf(o, "%s", parts[a]);
        CT one2; one2.k = qi(1);
        for (int l = 0; l < NM_MAXL; l++) one2.e[l] = q0();
        char *ob = ct_str(one2, x, shown, 1);
        sprintf(o, "%sO(%s)", np2 ? " + " : "", ob);
        return out;
    }
    nm_fail("the parallelogram did not separate the root after 12 steps");
    return NULL;
}
