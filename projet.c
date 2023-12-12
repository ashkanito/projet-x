#include <stdio.h>
#include <math.h>
#include <string.h>

//let's define universal constant :
#define G 6.67408*1e-11 //units : [m^3*kg^-1*s^-2]
#define M_earth 5.972*1e24 //units : [kg]
#define R_earth 6371000 // units : [m]
#define g0 9.81 // units : [m*s^-2]
#define I_sp 304 // units : [s] c'est pour un falcon 1 j'ai défini comme une constante mais on peut changer plus tard
#define P0 101325 //standard pressure (Pa), sea level
#define T0 288.15 //unités: K, niveau de la mer
#define R_sp 287.0 // gas constant units [J/(kg*K)]

//1)THRUST FORCE:
double F_thrust(double t, double t_burn, double thrust) {
    if (t < t_burn) {
        return thrust;
    }
    else {
        return 0;
    }
} 
//valeurs à utiliser pour la fusée Antares (Orbital Sciences Corporation):
//F_thrust = 3 265 000 N, t_burn = 235 s, mass= 260 700 kg

//2)ATMOSPHERIC DRAG
//2.1.Densité deuxième version (comme le code et pas le site)
double density(double z) {
    double h = 6356.766*z/(6356.766+z);
    double P;
    double Tm; //Tm=T*M0/M avec T: Temperature, M: masse molaire à cette hauteur et M0: masse molaire au niveau océan
    if (0<h && h<=11) {
        Tm = 288.15-6.5*h;
        P = 101325.0*pow((288.15/(288.15-6.5*h)), 34.1632/-6.5);
    } if (11<h && h<=20) {
        Tm = 216.65;
        P = 22632.06*exp(-34.1632*(h-11)/216.65);
    } if (20<h && h<=32) {
        Tm = 196.65+h;
        P = 5474.889*pow(216.65/(216.65+(h-20)), 34.1632);
    } if (32<h && h<=47) {
        Tm = 139.05+2.8*h;
        P = 868.0187*pow(228.65/(228.65+2.8*(h-32)),34.1632/2.8);
    } if (47<h && h<=51) {
        Tm = 270.65;
        P = 110.9063*exp(-34.1632*(h-47)/270.65);
    } if (51<h && h<=71) {
        Tm = 413.45/2.8*h;
        P = 66.93887*pow(270.65/(270.65-2.8*(h-51)),34.1632/-2.8);
    } if (71<h && h<=84.852) {
        Tm = 365.65-2.0*h;
        P = 3.956420*pow(214.65/(214.65-2*(h-71)),34.1632/-2);
    } if (86<z<=91) {
        return exp(-3.322622e-06*pow(z,3) + 9.111460E-04*pow(z,2) -0.2609971*z + 5.944694);
    } if (91<z<=100) {
        return exp(2.873405E-05*pow(z,3) -0.008492037*pow(z,2) + 0.6541179*z -23.62010);
    } if (100<z<=110) {
        return exp(-1.240774E-05*pow(z,4) + 0.005162063*pow(z,3) -0.8048342*pow(z,2) + 55.55996*z -1443.338);
    } if (110<z<=120) {
        return exp(-8.854164E-05*pow(z,3) + 0.03373254*pow(z,2) -4.390837*z + 176.5294);
    } if (120<z<=150) {
        return exp(3.661771E-07*pow(z,4) -2.154344E-04*pow(z,3) + 0.04809214*pow(z,2) -4.884744*z + 172.3597);
    } if (150<z<=200) {
        return exp(1.906032E-08*pow(z,4) -1.527799E-05*pow(z,3) + 0.004724294*pow(z,2) -0.6992340*z + 20.50921);
    } if (200<z<=300) {
        return exp(1.199282E-09*pow(z,4) -1.451051E-06*pow(z,3) + 6.910474E-04*pow(z,2) -0.1736220*z -5.321644);
    } if (300<z<=500) {
        return exp(1.140564E-10*pow(z,4) -2.130756E-07*pow(z,3) + 1.570762E-04*pow(z,2) -0.07029296*z -12.89844);
    }
    printf("Temperature: %f, pression: %f\n", Tm,P);
    return P/(287.053*Tm);  //R=287.053 J/kg.K = 8314/M avec M masse molaire de l'air
    //pareil a ce quil fait sur son code parceque 0.029/8.314 = 1/287.053 (a peu pres)
}

//2.2. Force Fdrag
double F_drag(double z, double v, double diametre, double Cd) {
    double Aire = 3.1415 *pow(diametre/2, 2); //cross-sectional area of the rocket
    return (1/2) * density(z) * pow(v,2) *Cd * Aire;
}
//valeurs à utiliser:
//Cd0 = 0.75, diametre = 3.07 (diametre du module Cygnus (fusee qui livre à ISS) en[m])
//QU?: vitesse?

//3) WAVE DRAG
//local velocity of sound : 
double v_sound(double ratio, double T){
    //printf("v_sound : %f ", sqrt(ratio * T * R_sp) );
    return sqrt(ratio * T * R_sp);
}
//mach number :
double M_fct(double rocket_speed, double sound_speed){
    return rocket_speed/sound_speed;
}
//prandtl-Glauert factor : 
double beta(double mach_number){
    // with the mach_number > 1 case for not having an error
    if(mach_number<1){
    return sqrt(1-pow(mach_number,2));
    }
    else if (mach_number == 1){
        mach_number = 0.9999;
        return sqrt(1-pow(mach_number,2));
    }
    else{
        return sqrt(pow(mach_number,2)-1);
    }
}
//prandtl-Glauert rule :
double c_d(double c_d0, double beta){
    return c_d0/beta;
}


double Lapse(double z) {
    double h = 6371*z/(6371+z); //aller voir site "comprehensive mathematical model of the Earth’s atmosphere"
    double L =0;
    if (h<11) {
        L = -6.5;	
    }
    else if (11<=h && h<20) {
        L = 0;	
    }
    else if (20<=h && h<32) {
        L = 1.0;	
    }
    else if (32<=h && h<47) {
        L = 2.8;	
    }
    else if (47<=h && h<51) {
        L = 0;	
    }
    else if (51<=h && h<71) {
       L = -2.8;	
    }
    else if (71<=h && h<=84.852) {
       L = -2.0;	
    }
    return L;
}

double temperature(double z) {
    double L = Lapse(z);
    return T0 - L*z;
}
//4) GRAVTATIONAL FORCE
//gravity in function of altitude :
double gravity(double altitude){
    return ((G*M_earth)/(pow(altitude+R_earth,2)));
}

//the mass changes because of the  : 
double m_variable(double time, double m0, double t_burn, double thrust){
    if(time<t_burn){  // on sait pas encore ce que c'est t_burn
        return (m0-((F_thrust(time,t_burn,thrust)/(g0*I_sp))*time)); // définir F_thrust
    }
    else{
        return (m0-((thrust/(g0*I_sp))*t_burn));
    }
}

//gravitational force : 
double F_g(double time, double altitude, double m0, double t_burn, double thrust){
    return (m_variable(time,m0,t_burn, thrust)*gravity(altitude));
}

// We are using Antares 200 series characteristics for our rocket.
int main(){
    // Euler-Cromer Method
    double dt = 0.1; // time step
    double z0 = 1; // initial altitude
    double v0 = 0; // initial velocity
    double m0 = 262600; // initial mass
    double FT = 3265000; //thrust of rocket
    double P = P0;
    double diameter = 3.9; // in meter
    double cd0 = 0.75;

    double v = v0;
    double z = z0;
    double m = m0;
    double M = m;
    double thrust = FT/m;
    double drag = 0;
    double ratio = 1.4; //heat capacity ratio
    double vsound;
    double Mach;
    double b;
    double cd = cd0;


    FILE *file; // File pointer
    char text[1000000] = "This is a sample text."; // Text to write to the file
    // Reset the text array to an empty string
    text[0] = '\0';


    printf("la gravité à altitude 200 mètre est %f\n", gravity(z));
    double t_final = 235; //tburn c'est 215
    double t_burn = 215;

    for (double t = 0; t < t_final; t+=0.1)
    {
        double g = gravity(z);
        m = m_variable(t, m0, t_burn, F_thrust(t, t_final, FT));
        double d = density(z);
        thrust = F_thrust(t, t_burn,FT)/m;              //thrust acceleration 
        vsound = v_sound(ratio, temperature(z));
        Mach = M_fct(v,vsound);
        b = beta(Mach);
        cd = c_d(cd0, b);
        drag = F_drag(z,v,diameter,cd);
        if(drag<0){         //flip drag force vector if rocket falls
            drag = drag*(-1);
        }

        v = v + (thrust - drag - g)*dt;
        z = z + v * dt;
        if(z<0){            //if rocket crashes
            break;
        }
        printf("altitude : %f meters, ", z);
        printf("mass : %f kg", m);
        printf(" speed : %f m/s\n", v);

        char temp[1000]; // Temporary string buffer
        sprintf(temp, "%f ; %f ; %f ; %f\n", z, m, v, t); // first one being altitude, then mass, then speed, then time
        strcat(text, temp);
        
        
    }


    // Open a file in write mode ('w')
    // w mode mean it writes a new file even if it already exists
    file = fopen("output.txt", "w");

    // Check if the file is successfully opened
    if (file == NULL) {
        printf("File could not be opened.");
        return 1; // Return an error code
    }

    // Write to the file
    fprintf(file, "%s\n", text);

    // Close the file
    fclose(file);

    printf("File created and data written successfully.\n");
    
 
}