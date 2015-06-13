#ifndef __microfact_bsdf_h__
#define __microfacet_bsdf_h__

#include "ggx_utils.h"

vector SphericalDirection(float sintheta;
                          float costheta;
                          float phi) {
    return set(sintheta*cos(phi), sintheta*sin(phi), costheta);
}

float f0_to_eta(float F0) {
    float f0sqrt = sqrt(clamp(F0, 0.0, 0.999));
    return ( 1.0 + f0sqrt) / ( 1.0 - f0sqrt);
}

vector f0_to_eta(vector F0) {
    vector f0sqrt = sqrt(clamp(F0, 0.0, 0.999));
    return ( 1.0 + f0sqrt) / ( 1.0 - f0sqrt);
}

vector f0_to_k(vector F0) {
    vector reflectance = clamp(F0, 0.0, 0.999);
    return 2.0 * sqrt(reflectance / (1.0 - reflectance));
}

vector eta_to_f0(vector eta;) {
    vector f0sqrt = (eta-1.0)/(eta+1.0);
    return f0sqrt*f0sqrt;
}


vector F_schlick(vector F0; float udoth) {
    float m = clamp(1.0-udoth, 0.0, 1.0);
    float m2 = m*m;
    return F0 + (1.0-F0)*m2*m2*m; 
}

vector F_cookTorrance(vector eta; float udoth) {
    vector c = abs(udoth);
    vector g = sqrt(eta*eta + c*c -1.0);
    vector tmp1 = (g-c)/(g+c);
    vector tmp2 = ((g+c)*c-1.0)/((g-c)*c+1.0);
    return 0.5*tmp1*tmp1*(1.0+tmp2*tmp2);
}

vector F_dielectric(vector eta; float udoth) {

    if (udoth < 0.0) eta = 1.0 / eta;
    vector c = abs(udoth);
    vector g = eta * eta - 1.0 + c * c;
    g = sqrt(g);
    vector A = (g - c) / (g + c);
    vector B = (c * (g + c) - 1.0) / (c * (g - c) + 1.0);
    vector ret = 0.5 * A * A * (1.0 + B * B);

    if (g.x <= 0.0) ret.x = 1.0;
    if (g.y <= 0.0) ret.y = 1.0;
    if (g.z <= 0.0) ret.z = 1.0;
    return ret;
}

vector F_conductor(vector eta; vector k; float udoth) {
    vector tmp = (eta*eta + k*k) * udoth*udoth;
    vector Rparl2 = (tmp - (2.f * eta * udoth) + 1.0) /
                      (tmp + (2.f * eta * udoth) + 1.0);
    vector tmp_f = eta*eta + k*k;
    vector Rperp2 =
        (tmp_f - (2.f * eta * udoth) + udoth*udoth) /
        (tmp_f + (2.f * eta * udoth) + udoth*udoth);
    return (Rparl2 + Rperp2) / 2.f;
}



struct F_parms {
    int fmeth = 0;
    vector F0 = 0.0;
    vector eta = 0.0;
    vector k = 0.0;
}

vector F_selector(const int mode; const F_parms fval; const float udoth) {

    if (mode == 0) {
        return 1.0;
    }
    if (mode == 1) {
        vector input = (fval.fmeth==0) ? fval.F0 : eta_to_f0(fval.eta);
        return F_schlick(input, udoth);
    }
    if (mode == 2) {
        vector input = (fval.fmeth==0) ? f0_to_eta(fval.F0) : fval.eta;
        return F_cookTorrance(input, udoth);
    }
    if (mode == 3) {
        vector input = (fval.fmeth==0) ? f0_to_eta(fval.F0) : fval.eta;
        return F_dielectric(input, udoth);
    }
    if (mode == 4) {
        vector input = (fval.fmeth==0) ? f0_to_eta(fval.F0) : fval.eta;
        vector input_k = (fval.fmeth==0) ? f0_to_k(fval.F0) : fval.k;
        return F_conductor(input, input_k, udoth);
    }
    return 1.0;
}

/*
struct G_parms {
    vector wo = 0.0;
    vector wi = 0.0;
    vector wh = 0.0;
    vector ng = 0.0;
    float alpha = 0.0;
}
*/

struct G_parms {
    vector wo = 0.0;
    vector wi = 0.0;
    vector wh = 0.0;
    vector ng = 0.0;
	float NdotWo = 0.0;
	float NdotWi = 0.0;
	float WodotWh = 0.0;
	float cosTheta = 0.0;
	float alpha = 0.0;
}

float G_smith1(float ndot; float alpha) {
    if (ndot <= 0.0) {
        return 0.0;
    }


    float a = alpha*alpha;
    float b = ndot*ndot;
    return 1.0/(ndot + sqrt(a + b - a*b));
}

float G_smith(const G_parms gval) {
    //float alpha = max(0.01, gval.alpha);
    return G_smith1(abs(gval.NdotWo), gval.alpha) * G_smith1(abs(gval.NdotWi), gval.alpha);
}

float G_schlick1(const float ndot; const float cosThetaSqr) {
    float k = sqrt(2.0*cosThetaSqr*0.31830988618379067154 );
    if (ndot <= 0.0)
        return 0.0;
    return ndot/(ndot - k*ndot + k);
}

float G_schlick(const G_parms gval) {
    float cosThetaSqr = gval.cosTheta;
    cosThetaSqr *= cosThetaSqr;
    float NdotWo = abs(gval.NdotWo);
    float NdotWi = abs(gval.NdotWi);
    return G_schlick1(NdotWo, cosThetaSqr) * G_schlick1(NdotWi, cosThetaSqr) /
        (4.0*NdotWo*NdotWi);
}


float G_cookTorrance(const G_parms gval) {
    float NdotWh = abs(dot(gval.ng, gval.wh));
    float NdotWo = abs(dot(gval.ng, gval.wo));
    float NdotWi = abs(dot(gval.ng, gval.wi));
    float WodotWh = abs(dot(gval.wo, gval.wh));
    return min(1.0, 2.0*(NdotWh/WodotWh)*min(NdotWi, NdotWo)) / 
            (4.0*max(0.01,NdotWo*NdotWi));
}

float G_neumann(const G_parms gval) {
    return 0.25/max(dot(gval.ng, gval.wo), 
                    dot(gval.ng, gval.wi));
}

float G_ward(const G_parms gval) {
    return 0.25/sqrt(dot(gval.ng, gval.wo) * 
                     dot(gval.ng, gval.wi));
}

float G_ashikhminShirley(const G_parms gval) {
    return 0.25/(abs(dot(gval.wi, gval.wh)) * 
                    max(abs(dot(gval.ng, gval.wi)),
                        abs(dot(gval.ng, gval.wo))));
}

float G_ashikhminPremoze(const G_parms gval) {
    float NdotWo = abs(dot(gval.ng, gval.wo));
    float NdotWi = abs(dot(gval.ng, gval.wi));
    return 0.25/(NdotWi+NdotWo - NdotWi*NdotWo);
}

float G_kurt(const G_parms gval) {
    float alpha = max(0.01, gval.alpha);
    return 1.0/(4.0*dot(gval.wi, gval.wh) * 
                    pow(abs(dot(gval.ng, gval.wi))* 
                        abs(dot(gval.ng, gval.wo)), alpha));
}

float G_kelemen(const G_parms gval) {
    return 1.0/(2.0*(1.0+abs(dot(gval.wo, gval.wi))));
}





float G_duer(const G_parms gval) {
    float cosThetaSqr = dot(gval.wh, gval.ng);
    cosThetaSqr *= cosThetaSqr;
    return dot(gval.wh, gval.wh) / (3.14159265358979323846 *cosThetaSqr*cosThetaSqr);
}

float G_beckmann(const G_parms gval) {
    float NdotWo = abs(dot(gval.ng, gval.wo));
    float NdotWi = abs(dot(gval.ng, gval.wi));
    float alpha = max(0.01, gval.alpha);
    float c = NdotWo/(alpha * sqrt(1.0 - (NdotWo*NdotWo)));
    if ( c >= 1.6 ) {
        return 1.0 / (4.0*NdotWo*NdotWi);
    }
    return (3.535*c + 2.181*c*c)/(1.0 + 2.276*c + 2.577*c*c) /
            (4.0*max(0.01,NdotWo*NdotWi));
}

float G_ggx(const G_parms gval) {
    float NdotWo = abs(dot(gval.ng, gval.wo));
    float NdotWi = abs(dot(gval.ng, gval.wi));
    float alpha = max(0.01, gval.alpha);
    float alphaSqr = alpha*alpha;

    return (2.0           ) / (NdotWo+sqrt(alphaSqr+(1.0-alphaSqr) *
                               NdotWo*NdotWo)) /
            (4.0           *max(0.01,NdotWi));
}

float G_schlickBeckmann(const G_parms gval)
{
    float NdotWo = abs(dot(gval.ng, gval.wo));
    float NdotWi = abs(dot(gval.ng, gval.wi));
    float alpha = max(0.01, gval.alpha);
    float k = alpha*0.79788456080286535588 ;

    return            1.0/(NdotWo*(1.0-k)+k) /
            (4.0           *max(0.01,NdotWi));
}

float G_schlickGGX(const G_parms gval)
{
    float NdotWo = abs(dot(gval.ng, gval.wo));
    float NdotWi = abs(dot(gval.ng, gval.wi));
    float alpha = max(0.01, gval.alpha);
    float k = alpha/2.0;

    return            1.0/(NdotWo*(1.0-k)+k) /
            (4.0           *max(0.01,NdotWi));
}

float G_selector(const int mode; const G_parms gval) {

    if (mode == 0)
        return 0.25f; 
    if (mode == 1)
        return G_smith(gval);
    if (mode == 2)
        return G_schlick(gval);
    if (mode == 3)
        return G_cookTorrance(gval);
    if (mode == 4)
        return G_neumann(gval);
    if (mode == 5)
        return G_ward(gval);
    if (mode == 6)
        return G_ashikhminShirley(gval);
    if (mode == 7)
        return G_ashikhminPremoze(gval);
    if (mode == 8)
        return G_kurt(gval);
    if (mode == 9)
        return G_kelemen(gval);
    if (mode == 10)
        return G_duer(gval);
    if (mode == 11)
        return G_beckmann(gval);
    if (mode == 12)
        return G_ggx(gval);
    if (mode == 13)
        return G_schlickBeckmann(gval);
    if (mode == 14)
        return G_schlickGGX(gval);
    return 1.0f;
}



struct D_blinn {

    float m_alpha = 0.0;
    float m_shiny = 0.0;
    float c_cosTheta = 0.0;

    float alphaToExpon(float alpha) {
        return 2.0/(alpha*alpha) - 2.0;
    }

    void init(float alpha) {
        this.m_alpha = alpha;
        this.m_shiny = this->alphaToExpon(alpha);
    }

    float D() {
        return                          pow(this.c_cosTheta, this.m_shiny);
    }

    float D(float cosTheta) {
        return                          pow(cosTheta, this.m_shiny);
    }


    float Rho() {



        return 1.0/(this.m_shiny+2.0);
    }

    float pdf(float cosTheta) {
        return (this.m_shiny+1) * this->D(cosTheta);
    }

    float pdf() {
        return (this.m_shiny+1) * this->D();
    }

    vector sample(float sx; float sy;) {

         float  cosTheta = pow(1-sx,1.0/(this.m_shiny+1));
         float  sinTheta = sqrt(max(0,1.0-(cosTheta*cosTheta)));
         float  phi = sy * 6.28318530717958647692 ;
         this.c_cosTheta = cosTheta;
         return SphericalDirection(sinTheta,cosTheta,phi);
    }
}



struct D_ggx {

    float m_alpha = 0.0;
    float m_alphaSqr = 0.0;
    float c_cosTheta = 0.0;

    void init(float alpha) {
        this.m_alpha = alpha;
        this.m_alphaSqr = alpha*alpha;
    }

    float D(float cosTheta) {
        float alphaSqrM1 = this.m_alphaSqr - 1.0;
        float cosThetaSqr = cosTheta*cosTheta;
        float tmp = 1.0+(alphaSqrM1)*cosThetaSqr;
        return 1.0 / (tmp*tmp);
    }

    float D() {
        return this->D(this.c_cosTheta);
    }

    float Rho() {


        return  0.5/this.m_alphaSqr;
    }

    float pdf(float cosTheta) {
        return this->D(cosTheta) / this->Rho();
    }

    float pdf() {
        return this->pdf(this.c_cosTheta);
    }

    vector sample(float sx; float sy) {
        float tanThetaSqr = this.m_alphaSqr * sx/(1.0-sx);
        float cosTheta = 1.0/sqrt(1.0+tanThetaSqr);
        float sinTheta = cosTheta*sqrt(tanThetaSqr);
        float phi = sy * 6.28318530717958647692 ;
        this.c_cosTheta = cosTheta;
        return SphericalDirection(sinTheta, cosTheta, phi);
    }
}



struct D_gtr {

    float m_alpha = 0.0;
    float m_alphaSqr = 0.0;
    float m_gamma = 0.0;
    float c_cosTheta = 0.0;

    void init(float alpha; float gamma) {






        this.m_alpha = max(0.0025, alpha);
        this.m_gamma = gamma;
        this.m_alphaSqr = alpha*alpha;
    }

    float D(float cosTheta) {
        float cosThetaSqr = cosTheta*cosTheta;
        return 1.0/pow(cosThetaSqr*(this.m_alphaSqr-1.0)+1.0,this.m_gamma);
    }

    float D() {
        return this->D(this.c_cosTheta);
    }

    float Rho() {


        if (this.m_gamma == 1.0) {
            return log(this.m_alpha)/(-1.0+this.m_alphaSqr);
        } 
        return -0.5 * (-1.0+pow(this.m_alpha,2.0*(1.0-this.m_gamma))) /
               ((-1.0+this.m_alphaSqr) * (-1.0+this.m_gamma));
    }

    float pdf(float cosTheta) {
        return this->D(cosTheta) / this->Rho();
    }

    float pdf() {
        return this->pdf(this.c_cosTheta);
    }

    vector sample(float sx; float sy) {
        float cosThetaSqr = 1.0/(1.0-this.m_alphaSqr);
        if ( this.m_gamma == 1.0 ) {
            cosThetaSqr *= 1.0-pow(this.m_alphaSqr, sx); 
        } else {
            float c_gamma = 1.0 - this.m_gamma;
            cosThetaSqr *= 1.0-pow(pow(this.m_alphaSqr,
                                   c_gamma)*(1.0-sx)+sx, 1.0/c_gamma);
        }
        float cosTheta = sqrt(cosThetaSqr);
        float sinTheta = sqrt(max(0.0,1.0-cosThetaSqr));
        float phi = sy * 6.28318530717958647692 ;
        this.c_cosTheta = cosTheta;
        return SphericalDirection(sinTheta, cosTheta, phi);
    }
}



struct D_beckmann {

    float m_alpha = 0.0;
    float m_alphaSqr = 0.0;


    float c_cosTheta = 0.0;

    void init(float alpha) {
        this.m_alpha = alpha;
        this.m_alphaSqr = alpha*alpha;
    }

    float D(float cosTheta) {
        float cosThetaSqr = cosTheta*cosTheta;
        float tanThetaSqr = (1.0 - cosThetaSqr)/cosThetaSqr;
        float cosTheta4th = cosThetaSqr*cosThetaSqr;
        return exp(-tanThetaSqr/this.m_alphaSqr) / (this.m_alphaSqr*cosTheta4th);
     }

    float D() {
        return this->D(this.c_cosTheta);
    }

    float Rho() {


        return  0.5;
    }

    float pdf(float cosTheta) {
        return this->D(cosTheta) / this->Rho();
    }

    float pdf() {
        return this->pdf(this.c_cosTheta);
    }

    vector sample(float sx; float sy) {
        float tanThetaSqr = -(this.m_alphaSqr) * log(1.0-sx);
        float cosTheta = 1.0/sqrt(1.0+tanThetaSqr);
        float tanTheta = sqrt(tanThetaSqr);
        float sinTheta = cosTheta * tanTheta;
        float phi = sy * 6.28318530717958647692 ;
        this.c_cosTheta = cosTheta;
        return SphericalDirection(sinTheta, cosTheta, phi);
    }
}





























struct D_parms {
    float alpha = 0.0;
    float gamma = 0.0;
}






struct D_generic {
    int m_type = 0;
    D_parms m_parms = {0.0, 0.0};

    D_blinn m_blinn = {0.0, 0.0, 0.0};
    D_gtr m_gtr = {0.0, 0.0, 0.0, 0.0};
    D_beckmann m_beckmann = {0.0, 0.0, 0.0};

    void init(int dtype; D_parms dparms) {
        this.m_type = dtype;
        this.m_parms = dparms;
        this.m_blinn = {0.0, 0.0, 0.0};
        this.m_gtr = {0.0, 0.0, 0.0, 0.0};
        this.m_beckmann = {0.0, 0.0, 0.0};

        if ( dtype == 0 ) {
            this.m_blinn->init(dparms.alpha);
            return;
        } else 
        if ( dtype == 1 ) {
            this.m_gtr->init(dparms.alpha, dparms.gamma);
            return;
        }
        if ( dtype == 2 ) {
            this.m_beckmann->init(dparms.alpha);
            return;
        }
        return;
    }

    float D(float cosTheta) {
                    if (this.m_type == 0 ) return this.m_blinn->D(cosTheta);
                    if (this.m_type == 1 ) return this.m_gtr->D(cosTheta);
                    if (this.m_type == 2 ) return this.m_beckmann->D(cosTheta);
        return 0.0;
    }

    float D() {
                    if (this.m_type == 0 ) return this.m_blinn->D();
                    if (this.m_type == 1 ) return this.m_gtr->D();
                    if (this.m_type == 2 ) return this.m_beckmann->D();
        return 0.0;
    }

    float Rho() {
                    if (this.m_type == 0 ) return this.m_blinn->Rho();
                    if (this.m_type == 1 ) return this.m_gtr->Rho();
                    if (this.m_type == 2 ) return this.m_beckmann->Rho();
        return 0.0;
    }

    float pdf(float cosTheta) {
                    if (this.m_type == 0 ) return this.m_blinn->pdf(cosTheta);
                    if (this.m_type == 1 ) return this.m_gtr->pdf(cosTheta);
                    if (this.m_type == 2 ) return this.m_beckmann->pdf(cosTheta);
        return 0.0;
    }

    float pdf() {
                    if (this.m_type == 0 ) return this.m_blinn->pdf();
                    if (this.m_type == 1 ) return this.m_gtr->pdf();
                    if (this.m_type == 2 ) return this.m_beckmann->pdf();
        return 0.0;
    }

    vector sample(float sx; float sy) {
        if ( this.m_type == 0 ) return this.m_blinn->sample(sx, sy);
        if ( this.m_type == 1 ) return this.m_gtr->sample(sx, sy);
        if ( this.m_type == 2 ) return this.m_beckmann->sample(sx, sy);
        return {0,0,1};
    }
}

#endif