/*! MINRESQLP solves the system of linear equations A*X=B or the least-squares problem min norm(B-A*X) if A is singular. The N-by-N matrix A must be symmetric or Hermitian, but need not be positive definite or nonsingular.  It may be double or single. The rhs vector B must have length N.  We assume it is pure real.
 
 This version has been adapted by T. Albash from the Matlab minresqlp.m file published on 28 Jun 2013 by Sou-Cheng (Terrya) Choi, CI, University of Chicago, and Michael Saunders, SOL, Stanford University and from the C++ version of tminres.hpp by Umberto Villa.
 
 @param A: NxN symmetric/Hermitian matrix
 
 @param x: is the initial guess (INPUT) and the solution of the linear system (OUTPUT).
 
 @param rtol: specifies a stopping tolerance.
 
 @param maxit: specifies the maximum number of iterations.
 
 @return flag:
 -1 (beta2=0)  B and X are eigenvectors of (A - SHIFT*I).
 0 (beta1=0)  B = 0.  The exact solution is X = 0.
 1 X solves the compatible (possibly singular) system (A - SHIFT*I)X = B  to the desired tolerance: RELRES = RNORM / (ANORM*XNORM + NORM(B)) <= RTOL, where   R = B - (A - SHIFT*I)X and RNORM = norm(R).
 2 X solves the incompatible (singular) system (A - SHIFT*I)X = B   to the desired tolerance:  RELARES = ARNORM / (ANORM * RNORM) <= RTOL, where  AR = (A - SHIFT*I)R and ARNORM = NORM(AR).
 3 Same as 1 with RTOL = EPS.
 4 Same as 2 with RTOL = EPS.
 5 X converged to an eigenvector of (A - SHIFT*I).
 6 XNORM exceeded MAXXNORM.
 7 ACOND exceeded ACONDLIM.
 8 MAXIT iterations were performed before one of the previous conditions was satisfied.
 9 The system appears to be exactly singular.  XNORM does not yet exceed MAXXNORM, but would if further iterations were performed.
 
 */

#include <limits>
#include <memory>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>



void SymOrtho(double a, double b, double &c, double &s, double &r)
{
    double absa  = abs(a);
    double absb  = abs(b);
    double t;
    int signa;
    int signb;
    
    if(a > 0){
        signa = 1;
    }
    else if(a<0){
        signa = -1;
    }
    else{
        signa = 0;
    }
    
    if(b > 0){
        signb = 1;
    }
    else if(b<0){
        signb = -1;
    }
    else{
        signb = 0;
    }
    
    //if a or b equals 0
    if(b==0){
        if(a==0){
            c=1;
        }
        else{
            c = signa;
        }
        s=0;
        r=absa;
    }
    else if(a==0){
        c = 0;
        s = signb;
        r = absb;
    }
    else{
        //both a and b are non-zero
        if(absb>absa){
            t = a/b;
            s = signb/(std::sqrt(1 + pow(t,2.0)));
            c = s*t;
            r = b/s; // computationally better than d = a / c since |c| <= |s|
        }
        else{
            t = b/a;
            c = signa/(std::sqrt(1 + pow(t,2.0)));
            s = c*t;
            r = a/c; // computationally better than d = b / s since |s| <= |c|
        }
    }
} //SymOrtho

template< typename Operator, typename Vector>
int
MINRESQLP(const Operator &A, Vector &x, const Vector &b, double &rtol, int maxit, int &iter)
{
    typedef Operator operator_Type;
    typedef Vector   vector_Type;
    
    bool done;
    int flag;
    
    //in what follows, we will assume no preconditioner M is provided
    
    
    // Set up {beta1, p, v} for the first Lanczos vector v1.
    x=0;
    const std::unique_ptr<Vector> r2(x.Clone());  //r2 = b
    const std::unique_ptr<Vector> r3(x.Clone());  //r3 = b
    (*r2) = b;
    (*r3) = b;
    
    
    double beta1(0.0);
    beta1 = InnerProduct(*r2,*r3); //beta1 = b'*b
    
    
    if(beta1 == 0.0){
        done = true;
    }
    else{
        done = false;
        beta1 = std::sqrt(beta1); // Normalize y to get v1 later
    }
    
    iter=0;
    int flag0(-2), QLPiter(0), lines(1),headlines (20);
    double beta(0.0), tau(0.0), taul(0.0), phi(beta1);
    double betan(beta1), gmin(0.0), cs(-1.0), sn(0.0);
    double cr1(-1.0), sr1(0.0), cr2(-1.0), sr2(0.0);
    double dltan(0.0), eplnn(0.0), gama(0.0), gamal(0.0);
    double gamal2(0.0), eta(0.0), etal(0.0), etal2(0.0);
    double vepln(0.0), veplnl(0.0), veplnl2(0.0), ul3(0.0);
    double ul2(0.0), ul(0.0), u(0.0), rnorm(betan);
    double xnorm(0.0), xl2norm(0.0), Axnorm(0.0);
    double Anorm(0.0), Acond(-1.0);
    double gama_QLP(0.0), gamal_QLP(0.0), vepln_QLP(0.0), gminl(0.0);
    double u_QLP(0.0), ul_QLP(0.0);
    double relres   = rnorm / (beta1 + 1e-50);       // Safeguard for beta1 = 0
    double relresl(0.0), relAresl(0.0);
    
    const std::unique_ptr<Vector> xl2(x.Clone()), w(x.Clone()), wl(x.Clone()), wl2(x.Clone()), r1(x.Clone()), y(x.Clone());
    ( *xl2 ) = (* w ) = ( *wl ) = ( *wl2) = ( *r1 ) = (*y) = 0.0;
    
    const std::unique_ptr<Vector> v(x.Clone());
    
    
    flag = flag0;
    
    double betal,alfa,pnorm,gamal3,gama_tmp, taul2,dbar,dlta,epln,gbar,dlta_QLP,dlta_tmp,ul4,xnorm_tmp,xnorml,gamal_tmp,t1,t2,epsx,Arnorml,rootl,Acondl,gminl2,abs_gama,Anorml,rnorml,Arnorm,relAres;
    
    int likeLS,Miter;
    //set default values
    double realmin = std::numeric_limits<double>::min( );
    double eps(std::numeric_limits<double>::epsilon());
    double maxxnorm = pow(10.0,7);
    double shift = 0;
    double Acondlim = pow(10.0,15);
    double TranCond = pow(10.0,7);
    
    if(!done)
    {
        /* Main Iteration */
        while(flag==flag0 && iter < maxit){
            
            // Lanczos
            {
                iter  = iter + 1;
                betal = beta;
                beta = betan;
                
                //v = r3*(1/beta);
                (*v)  = (*r3);
                v->Scale(1/beta);
                
                
                A.Apply(*v,*r3);
                //for us shift will always equal 0
                /*if(shift!=0){
                 add(*r3, -shift, *v, *r3);
                 }*/
                if(iter>1){
                    //r3 = r3 - (beta/betal)*r1;
                    add(*r3, -beta/betal, *r1, *r3);
                }
                
                alfa = InnerProduct(*v, *r3);    // alphak
                
                //r3   = r3 - (alfa/beta)*r2;     r1 = r2;     r2 = r3;
                add(*r3, -alfa/beta, *r2, *r3);
                (*r1) = (*r2);
                (*r2) = (*r3);
                
                //no preconditioner
                betan = InnerProduct(*r3,*r3);
                betan = std::sqrt(betan);
                if(iter==1){
                    if(betan==0){
                        if(alfa==0){
                            flag = 0;
                            break;
                        }
                        else{
                            flag = -1;
                            x = b;
                            x.Scale(1/alfa);
                        }
                    }
                }
                
                if(iter<=2){
                    pnorm = std::sqrt(alfa*alfa + betan*betan);
                }
                else{
                    pnorm = std::sqrt(betal*betal + alfa*alfa + betan*betan);
                    
                }
            }//Lanczos
            
            //Apply previous left reflection Q_{k-1}
            {
                dbar  = dltan;
                dlta  = cs*dbar + sn*alfa;
                epln     = eplnn;
                gbar  = sn*dbar - cs*alfa;
                eplnn    = sn*betan;
                dltan = -cs*betan;
                dlta_QLP = dlta;
                
            }
            //Compute the current left reflection Q_k
            {
                gamal3 = gamal2;
                gamal2 = gamal;
                gamal    = gama;
                SymOrtho(gbar, betan,cs,sn,gama);
                gama_tmp = gama;
                taul2  = taul;
                taul   = tau;
                tau      = cs*phi;
                Axnorm = std::sqrt(Axnorm*Axnorm + tau*tau);
                phi = sn*phi;
                
            }
            //Apply the previous right reflection P{k-2,k}
            {
                if(iter > 2){
                    veplnl2  = veplnl;
                    etal2 = etal;
                    etal = eta;
                    dlta_tmp = sr2*vepln - cr2*dlta;
                    veplnl   = cr2*vepln + sr2*dlta;
                    dlta     = dlta_tmp;
                    eta = sr2*gama;
                    gama = -cr2*gama;
                } //if
                
            }
            // Compute the current right reflection P{k-1,k}, P_12, P_23,...
            {
                if(iter > 1){
                    SymOrtho(gamal, dlta,cr1, sr1, gamal);
                    vepln =   sr1*gama;
                    gama  = - cr1*gama;
                }// if
            }
            //Update xnorm
            {
                xnorml = xnorm;
                ul4 = ul3;
                ul3   = ul2;
                if (iter > 2){
                    ul2 = (taul2 - etal2*ul4 - veplnl2*ul3) / gamal2;
                }
                if (iter > 1){
                    ul = ( taul  - etal *ul3 - veplnl *ul2) / gamal;
                }
                
                xnorm_tmp = std::sqrt(xl2norm*xl2norm + ul2*ul2 +  ul*ul);
                
                if(relresl >= relAresl){
                    likeLS = 1;
                }
                else{
                    likeLS = 0;
                }
                if (abs(gama) > realmin && xnorm_tmp < maxxnorm){
                    u = (tau - eta*ul2 - vepln*ul) / gama;
                    if (std::sqrt(xnorm_tmp*xnorm_tmp +u*u) > maxxnorm && likeLS==1){
                        u = 0;
                        flag = 6;
                    }
                }
                else{
                    u = 0;
                    flag = 9;
                }
                
                xl2norm = std::sqrt(xl2norm*xl2norm + ul2*ul2);
                xnorm   = std::sqrt(xl2norm*xl2norm + ul*ul + u*u);
            }
            
            //Update w. Update x except if it will become too big
            {
                if ((Acond < TranCond) && flag != flag0 && QLPiter==0){
                    //MINRES updates
                    (*wl2) = (*wl);
                    (*wl) = (*w);
                    //w   = (v - epln*wl2 - dlta_QLP*wl) * (1/gama_tmp);
                    add(-epln, *wl2, -dlta_QLP, *wl, *w);
                    add((1.0/gama_tmp), *v, *w, *w);
                    if (xnorm < maxxnorm){
                        //x = x + tau*w;
                        add(x, tau, *w, x);
                    }
                    else{
                        flag = 6;
                    }
                }//if
                else{
                    //MINRES-QLP updates
                    QLPiter = QLPiter + 1;
                    if (QLPiter == 1){
                        //xl2 = zeros(n,1);
                        if  (iter > 1){
                            // construct w_{k-3}, w_{k-2}, w_{k-1}
                            if (iter > 3){
                                //wl2 = gamal3*wl2 + veplnl2*wl + etal*w;
                                add(gamal3, *wl2, veplnl2, *wl, *wl);
                                add(*wl,etal,*w,*wl);
                            } // w_{k-3}
                            if (iter > 2){
                                //wl = gamal_QLP*wl + vepln_QLP*w;
                                add(gamal_QLP,*wl,vepln_QLP,*w,*wl);
                            }// w_{k-2}
                            //w = gama_QLP*w;
                            w->Scale(gama_QLP);
                            //xl2 = x - wl*ul_QLP - w*u_QLP;
                            add(x,-ul_QLP,*wl,*xl2);
                            add(*xl2,-u_QLP,*w,*xl2);
                            
                        } //if
                    }//if
                    
                    if (iter == 1){
                        (*wl2) = (*wl);
                        //wl = v*sr1;
                        (*wl) = (*v);
                        wl->Scale(sr1);
                        
                        //w  = -v*cr1;
                        (*w) = (*v);
                        w->Scale(-cr1);
                    }
                    else if(iter == 2){
                        (*wl2) = (*wl);
                        //wl  = w*cr1 + v*sr1;
                        add(cr1,*w,sr1,*v,*wl);
                        
                        //w   = w*sr1 - v*cr1;
                        add(sr1,*w,-cr1,*v,*w);
                    }
                    else{
                        (*wl2) = (*wl);
                        (*wl) = (*w);
                        
                        //w  = wl2*sr2 - v*cr2;
                        add(sr2,*wl2,-cr2,*v,*w);
                        
                        //wl2 = wl2*cr2 + v*sr2;
                        add(cr2,*wl2,sr2,*v,*wl2);
                        
                        //v  = wl *cr1 + w*sr1;
                        add(cr1,*wl,sr1,*w,*v);
                        
                        //w   = wl *sr1 - w*cr1;
                        add(sr1,*wl,-cr1,*w,*w);
                        
                        (*wl) = (*v);
                    }
                    //xl2 = xl2 + wl2*ul2;
                    add(*xl2,ul2,*wl2,*xl2);
                    
                    //x   = xl2 + wl *ul + w*u;
                    add(*xl2,ul,*wl,x);
                    add(x,u,*w,x);
                }
                
            }
            
            //Compute the next right reflection P{k-1,k+1}
            {
                gamal_tmp = gamal;
                SymOrtho(gamal,eplnn,cr2,sr2,gamal);
            }
            
            //Store quantities for transfering from MINRES to MINRES-QLP
            {
                gamal_QLP = gamal_tmp;
                vepln_QLP = vepln;
                gama_QLP = gama;
                ul_QLP = ul;
                u_QLP = u;
            }
            
            //Estimate various norms
            {
                abs_gama = abs(gama);
                Anorml = Anorm;
                Anorm = max(max(max(Anorm, gamal), abs_gama), pnorm);
                
                if (iter == 1){
                    gmin = gama;
                    gminl = gmin;
                }
                else if (iter > 1){
                    gminl2 = gminl;
                    gminl = gmin;
                    gmin = min(min(gminl2, gamal), abs_gama);
                }
                Acondl   = Acond;
                Acond   = Anorm/gmin;
                rnorml   = rnorm;
                relresl = relres;
                
                if (flag != 9){
                    rnorm = phi;
                }
                relres   = rnorm / (Anorm*xnorm + beta1);
                rootl    = std::sqrt(gbar*gbar + dltan*dltan);
                Arnorml  = rnorml*rootl;
                relAresl = rootl / Anorm;
            }
            
            // See if any of the stopping criteria are satisfied.
            epsx = Anorm*xnorm*eps;
            if ( (flag == flag0) || (flag == 9)){
                t1 = 1 + relres;
                t2 = 1 + relAresl;
                if (iter     >= maxit){
                    flag = 8;
                }// Too many itns
                if (Acond    >= Acondlim) {
                    flag = 7;
                    
                } // Huge Acond
                if (xnorm    >= maxxnorm) {
                    flag = 6;
                    
                } // xnorm exceeded its limit
                if (epsx     >= beta1) {
                    flag = 5;
                    
                } // x is an eigenvector
                if (t2       <= 1) {
                    flag = 4;
                    
                } // Accurate LS solution
                
                if (t1       <= 1) {
                    flag = 3;
                    
                } // Accurate Ax=b solution
                if (relAresl <= rtol) {
                    flag = 2;
                    
                } // Good enough LS solution
                if (relres   <= rtol ) {
                    flag = 1;
                    
                } // Good enough Ax=b solution
                
                if ((flag == 2) || (flag == 4) || (flag == 6 && likeLS) || (flag == 7)) {  // Possibly singular
                   iter  = iter - 1;
                   Acond = Acondl;   rnorm = rnorml;   relres = relresl;
                }
                
            }
            
            
            
        }//while
        
        Miter = iter - QLPiter;
        
        //Compute final quantities directly.
        //r1      = b - minresxxxA(A,x) + shift*x;   % r1 is workspace for residual vector
        A.Apply(x, *r1);
        add(*r1,-shift,x,*r1);
        subtract(b, *r1, *r1);
        
        rnorm   = std::sqrt(InnerProduct(*r1,*r1));
        
        //y = minresxxxA(A,r1) - shift*r1
        A.Apply(*r1, *y);
        add(*y,-shift,*r1,*y);
        
        Arnorm  = std::sqrt(InnerProduct(*y,*y));
        xnorm   = std::sqrt(InnerProduct(x,x));
        relres  = rnorm / (Anorm*xnorm + beta1);
        relAres = 0;
        if (rnorm > realmin){
          relAres = Arnorm / (Anorm*rnorm);
        }
    }
        
    return flag;
}

/*
 @param maxnorm: is an upper bound on NORM(X).  Default MAXXNORM=1e7.
 
 @param M: preconditioner.  M must be positive definite and symmetric or Hermitian. If M is a null pointer, preconditioning is not used.
 
 @param shift: If shift != 0, solves (A - shift*I)x = b, or the corresponding least-squares problem if shift is an
 eigenvalue of A (in the latter case, A - shift*I is singular).
 
 @param Acondlim: is an upper bound on ACOND, an estimate of COND(A). Default ACONDLIM=1e15.
 
 @param TranCond: is a real scalar >= 1.
 If TRANCOND>1,        a switch is made from MINRES iterations to MINRES-QLP iterationsd when ACOND >= TRANCOND.
 If TRANCOND=1,        all iterations will be MINRES-QLP iterations.
 If TRANCOND=ACONDLIM, all iterations will be conventional MINRES iterations (which are slightly cheaper).
 Default TRANCOND=1e7.
 
 @param show: specifies the printing option.  If SHOW=true,  an iteration log will be output.  If SHOW=false, the log is suppressed. Default SHOW=true.
 */
