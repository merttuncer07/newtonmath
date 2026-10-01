/* Independent forms, Hermite identities, conjugate surds and closed poles. */
#include "nm.h"
#include <stdio.h>
#include <string.h>
#include <setjmp.h>
extern jmp_buf nm_on_error;
extern char nm_error_msg[512];
static int checks,failed,back,added;
static Q qi(int n) { return q_from_z(z_from_i64(n)); }
static C N(int n) { return c_const(qi(n)); }
static void ck(int ok,const char *s) { checks++;if(!ok) {failed++;printf("FAIL %s\n",s);} }
static void refusal(R r,int x,const char *why) {
    jmp_buf save;memcpy(save,nm_on_error,sizeof save);
    if(!setjmp(nm_on_error)) { (void)integ_rational(r,x);ck(0,"expected refusal"); }
    else ck(strstr(nm_error_msg,why)!=NULL,"refusal reason");
    memcpy(nm_on_error,save,sizeof save);
}
int main(void) {
    if(setjmp(nm_on_error)) {printf("error: %s\n",nm_error_msg);return 1;}
    int v=letter_index("x",1);C x=c_letter(v),a=c_letter(letter_index("a",1));
    C xx=c_mul(x,x),xm=c_sub(x,N(1)),xp=c_add(x,N(1)),circle=c_add(xx,N(1));
    C quad=c_add(circle,x),cube=c_add(c_pow_int(x,3),N(1));
    C den[]={c_sub(xx,N(1)),circle,circle,quad,c_sub(xx,N(1)),c_pow_int(xm,2),
        c_mul(x,c_pow_int(xp,2)),c_add(c_add(xx,c_scale(x,qi(2))),N(5)),c_sub(xx,N(2)),
        c_sub(c_pow_int(x,4),N(1)),c_add(x,a),cube,
        c_pow_int(circle,2),c_pow_int(quad,3),c_mul(circle,c_add(xx,N(2))),
        c_pow_int(c_add(x,a),3),c_mul(c_add(xx,N(2)),c_add(c_add(xx,x),N(3)))};
    C num[]={N(1),x,N(1),N(1),cube,N(1),N(1),c_add(c_scale(x,qi(2)),N(3)),N(1),N(1),xx,N(1),N(1),N(1),N(1),xx,N(1)};
    Integral keep={0},keep_param={0};Apart saved={0},saved_param={0};
    for(int i=0;i<17;i++) {
        R r=r_make(num[i],den[i]); Integral z=integ_rational(r,v);
        printf("%d: %s\n",i+1,integ_to_str(z));
        back++;ck(r_equal(integ_derivative(z),r),"derivative back");
        Apart ap=integ_apart(r,v);added++;ck(r_equal(apart_sum(ap),r),"apart back");
        if(i==0) {ck(z.n==2,"two hyperbolic areas");keep=integ_persist(z);saved=apart_persist(ap);}
        if(i==2) ck(z.n==1 && z.term[0].kind==AREA_ATAN && r_equal(z.term[0].coef,r_from_c(N(1))) && c_equal(z.term[0].poly,x),"unit circle");
        if(i==5) ck(!z.n && r_equal(z.rational,r_make(N(-1),xm)),"repeated linear is rational");
        if(i==9) ck(ap.n==3,"quartic three partial fractions");
    }
    /* A repeated irreducible cubic may integrate rationally before factoring. */
    C cubic=c_add(c_add(c_pow_int(x,3),x),N(1));
    R rat=r_make(N(1),cubic);
    Integral primitive={v,0,rat,NULL}; R df=integ_derivative(primitive);
    Integral found=integ_rational(df,v);back++;
    ck(found.n==0 && r_equal(found.rational,rat),"Hermite precedes unsupported factorization");
    refusal(r_make(N(1),cubic),v,"Rothstein-Trager");
    refusal(r_make(N(1),c_add(xx,a)),v,"sign unknown");
    ck(elim_has_root_closed(x,v,qi(-1),qi(1)),"interior pole");
    ck(elim_has_root_closed(c_pow_int(xm,2),v,qi(0),qi(2)),"even multiplicity interior pole");
    ck(elim_has_root_closed(xm,v,qi(1),qi(2)),"left endpoint pole");
    ck(elim_has_root_closed(xm,v,qi(0),qi(1)),"right endpoint pole");
    ck(elim_has_root_closed(c_sub(xx,N(2)),v,qi(2),qi(0)),"surd pole and reversed interval");
    ck(!elim_has_root_closed(circle,v,qi(-2),qi(2)),"no real pole");
    ck(!elim_has_root_closed(xm,v,qi(2),qi(2)),"regular zero width interval");
    /* Parameter-field coefficients must remain R, including a+b in a
     * denominator. Check the primitives independently, not only by reversal. */
    int bi=letter_index("b",1); C u=c_add(a,c_letter(bi)),linear=c_add(c_mul(u,x),N(1));
    for(int i=0;i<3;i++) {
        R r=r_make(c_pow_int(x,i),u);Integral z=integ_rational(r,v);back++;
        R want=r_make(c_pow_int(x,i+1),c_scale(u,qi(i+1)));
        ck(z.n==0 && r_equal(z.rational,want),"parameter-independent denominator primitive");
        ck(r_equal(integ_derivative(z),r),"parameter-independent derivative back");
        Apart ap=integ_apart(r,v);added++;ck(ap.n==1 && r_equal(apart_sum(ap),r),"parameter-independent apart");
    }
    for(int i=0;i<3;i++) {
        R r=r_make(c_pow_int(x,i),linear);Integral z=integ_rational(r,v);back++;
        R logcoef=r_make(N(i==1?-1:1),c_pow_int(u,i+1));
        R rational=r_from_c(N(0));
        if(i==1) rational=r_make(x,u);
        if(i==2) rational=r_sub(r_make(xx,c_scale(u,qi(2))),r_make(x,c_pow_int(u,2)));
        ck(z.n==1 && z.term[0].kind==AREA_LOG && c_equal(z.term[0].poly,linear) &&
            r_equal(z.term[0].coef,logcoef) && r_equal(z.rational,rational),"general linear parameter primitive");
        ck(r_equal(integ_derivative(z),r),"general linear parameter derivative back");
        Apart ap=integ_apart(r,v);added++;ck(r_equal(apart_sum(ap),r),"general linear parameter apart");
        if(i==2) {keep_param=integ_persist(z);saved_param=apart_persist(ap);}
    }
    C scales[]={a,u};
    for(int i=0;i<2;i++) {
        R r=r_make(c_scale(x,qi(-2)),c_mul(scales[i],c_pow_int(circle,2)));
        Integral z=integ_rational(r,v);back++;
        ck(z.n==0 && r_equal(z.rational,r_make(N(1),c_mul(scales[i],circle))),"Hermite retains parameter content");
        ck(r_equal(integ_derivative(z),r),"parameter content derivative back");
        Apart ap=integ_apart(r,v);added++;ck(r_equal(apart_sum(ap),r),"parameter content apart");
        refusal(r_make(N(1),c_mul(scales[i],circle)),v,"sign unknown");
    }
    R repeated=r_make(N(1),c_pow_int(linear,2));
    Integral linear_repeat=integ_rational(repeated,v);back++;
    ck(linear_repeat.n==0 && r_equal(linear_repeat.rational,r_make(N(-1),c_mul(u,linear))),"repeated general linear primitive");
    ck(r_equal(integ_derivative(linear_repeat),repeated),"repeated general linear derivative back");
    Apart repeated_apart=integ_apart(repeated,v);added++;
    ck(r_equal(apart_sum(repeated_apart),repeated),"repeated general linear apart");
    R improper=r_make(xx,c_pow_int(linear,2));
    Integral improper_area=integ_rational(improper,v);back++;
    R improper_part=r_sub(r_make(x,c_pow_int(u,2)),r_make(N(1),c_mul(c_pow_int(u,3),linear)));
    ck(improper_area.n==1 && r_equal(improper_area.rational,improper_part) &&
       r_equal(improper_area.term[0].coef,r_make(N(-2),c_pow_int(u,3))),"improper repeated general linear primitive");
    ck(r_equal(integ_derivative(improper_area),improper),"improper repeated general linear derivative back");
    Apart improper_apart=integ_apart(improper,v);added++;
    ck(r_equal(apart_sum(improper_apart),improper),"improper repeated general linear apart");
    /* Sample roots must not depend on the moderate-coefficient search bound. */
    C large=c_mul(c_sub(xx,N(1)),c_add(xx,c_const(q_from_z(z_from_i64(10000000000000LL)))));
    C fs[8],product=N(1);int nf=elim_conic_factors(large,v,fs,8);
    for(int i=0;i<nf;i++) product=c_mul(product,fs[i]);
    ck(nf==3,"large coefficients retain rational sample roots");
    ck(c_equal(product,large),"large coefficient factors multiplied back");
    arena_reset();
    x=c_letter(v);back++;ck(r_equal(integ_derivative(keep),r_make(N(1),c_sub(c_mul(x,x),N(1)))),"persistent integral");
    added++;ck(r_equal(apart_sum(saved),r_make(N(1),c_sub(c_mul(x,x),N(1)))),"persistent apart");
    u=c_add(c_letter(letter_index("a",1)),c_letter(bi));
    R param=r_make(c_mul(x,x),c_add(c_mul(u,x),N(1)));
    back++;ck(r_equal(integ_derivative(keep_param),param),"persistent parameter coefficient");
    added++;ck(r_equal(apart_sum(saved_param),param),"persistent parameter apart");
    printf("integ: %d differentiated back, %d added back\n",back,added);
    printf("unit integ: %d checks, %d failed\n",checks,failed);return failed!=0;
}
