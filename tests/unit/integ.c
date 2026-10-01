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
    Integral keep={0};Apart saved={0};
    for(int i=0;i<17;i++) {
        R r=r_make(num[i],den[i]); Integral z=integ_rational(r,v);
        printf("%d: %s\n",i+1,integ_to_str(z));
        back++;ck(r_equal(integ_derivative(z),r),"derivative back");
        Apart ap=integ_apart(r,v);added++;ck(r_equal(apart_sum(ap),r),"apart back");
        if(i==0) {ck(z.n==2,"two hyperbolic areas");keep=integ_persist(z);saved=apart_persist(ap);}
        if(i==2) ck(z.n==1 && z.term[0].kind==AREA_ATAN && c_equal(z.term[0].coef,N(1)) && c_equal(z.term[0].poly,x),"unit circle");
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
    arena_reset();
    x=c_letter(v);back++;ck(r_equal(integ_derivative(keep),r_make(N(1),c_sub(c_mul(x,x),N(1)))),"persistent integral");
    added++;ck(r_equal(apart_sum(saved),r_make(N(1),c_sub(c_mul(x,x),N(1)))),"persistent apart");
    printf("integ: %d differentiated back, %d added back\n",back,added);
    printf("unit integ: %d checks, %d failed\n",checks,failed);return failed!=0;
}
