#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>




//let's define universal constant :
#define G 6.67408*1e-11 //units : [m^3*kg^-1*s^-2]
#define M_earth 5.972*1e24 //units : [kg]
#define R_earth 6371000 // units : [m]
#define g0 9.81 // units : [m*s^-2]
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


//2)ATMOSPHERIC DRAG
//The calculations are based on the formulas provided by http://www.braeunig.us/space/atmmodel.htm
double density(double z) {
    double h = 6356.766*z/(6356.766+z);
    double P;
    double Tm; //Tm=T*M0/M with T: Temperature, M: specific molar mass at given altitude, M0: molar mass at the level of the ocean
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
    return P/(287.053*Tm);   //R=287.053 J/kg.K = 8314/M with M the air molar mass
}

//2.2. Force Fdrag
double F_drag(double z, double v, double diametre, double Cd) {
    double Aire = 3.1415 *pow(diametre/2, 2); //cross-sectional area of the rocket
    return (1/2) * density(z) * pow(v,2) *Cd * Aire;
}


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
    double h = 6371*z/(6371+z); //http://www.braeunig.us/space/atmmodel.htm
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
double m_variable(double time, double FT, double I_sp){
    return (FT/(g0*I_sp))*time; 
}



//gravitational force : 
//double F_g(double time, double altitude, double m0, double t_burn, double thrust, double I_sp){
  //  return (m_variable(time,m0,t_burn, thrust, I_sp)*gravity(altitude));
//}

// We are using Antares 130 series characteristics for our rocket.
// second stage is for the Castor30B characteristics


//5)coordinates of ISS, read on a text file
struct Coordinates_min_dist {
    double distance;
    double lat;
    double lon;
    double alt;
    int hour;
    int min;
    int sec;
};

int readfile(char *filename, int length, struct Coordinates_min_dist *conc) { //see training midterm
  /*
    Reads in data file called "filename"
    up to "length" number of rows
    into an array of structures of struct TimeSeries called "conc".

    Returns an integer corresponding to number of rows read in.
  */

  FILE *fid = fopen(filename, "r");
  if(fid == NULL) return -1;

  int n=0;    
  char buffer[100];  
  while (fgets(buffer, 100, fid) != NULL) {
    if(n >= length) break;
    sscanf(buffer, "%lf %lf %lf %lf %d %d %d", &conc[n].distance, &conc[n].lat, &conc[n].lon, &conc[n].alt, &conc[n].hour, &conc[n].min, &conc[n].sec);
    n++;
  }

  return n;
}

struct time {
    int hour;
    int min;
    int sec;
};

struct time when_to_launch(int hour, int min, int sec, int travel_time){ //travel time as an integer
    struct time launch_time;
    int sec_iss = hour*3600 + min*60 + sec;
    int launch_in_sec = sec_iss - travel_time; //launch in seconds as an integer

    launch_time.hour = launch_in_sec/3600;
    int p = launch_in_sec%3600;
    launch_time.min = p/60;
    launch_time.sec = p%60;

    return launch_time;
}

int main(){
    // We asked ChatGPT3.5 to run our python filed (coord_skyfield) from our c file
    // Run Python script using system()
    int status = system("python coord_skyfield.py");
    
    if (status == -1) {
        printf("Failed to execute the command.\n");
        return 1;
    } else {
        printf("Python script executed successfully.\n");
    }

    //extract the coordinates of ISS at minimal distance
    struct Coordinates_min_dist conc[50];
    double n = readfile("output_coords.txt", 1, conc);
    double dist_min = conc[0].distance; 
    double lat_min = conc[0].lat; //the notation "min" is used to describe the variables when the distance is minimal, it is not the minmimum value of each variable
    double lon_min = conc[0].lon;
    double alt_min = conc[0].alt;
    int hour_min = conc[0]. hour;
    int min_min = conc[0].min;
    int sec_min = conc[0].sec;

    double alt_min_km=alt_min*1000;

    // Euler-Cromer Method
    double dt = 0.5; // time step
    double z0 = 1; // initial altitude
    double v0 = 0; // initial velocity
    double m_fuel1MAX = 242000; // initial mass of fuel and oxidizer of stage one
    double m_fuel2MAX = 12887; // initial mass of fuel and oxidizer of stage two
    double m1 = 18700; // dry mass of stage one
    double m2 = 1083; //dry mass of second stage 
    double m_charge = 8000; // mass of the payload
    double m_finale = 13330; // final dry mass of the rocket
    double m_tot = m_fuel1MAX + m_fuel2MAX + m1 + m2 + m_charge + m_finale; //total mass of the rocket (imaginons qu'on rajoute 8 tonnes de marchandises)
    double FT1 = 3000357.2; //calculated thrust force of rocket stage 1
    double FT2 = 302615.2; //calculated thrust force of rocket stage 2
    double FT1_MAX = 3265000; //maximum thrust of rocket stage 1 
    double FT2_MAX = 396300; //maximum thrust of rocket stage 2 
    double P = P0;
    double diameter = 3.9; // in meter
    double cd0 = 0.75; //aerodynamic coefficients (approximately for a rocket)
    double I_sp1 = 297; // units : [s] specific impulse of rocket (sea level)
    double I_sp2 = 304; // units : [s] specific impulse of castor30B
    double t_final = 600; // choice of the final time
    double t_min = t_final; // to find the minimum time
    
    // The values of thrust and the fuel mass that we are going to use for the plot:
        double FT1_plot;
        double FT2_plot;
        double m_fuel1_plot;
        double m_fuel2_plot;


    //create a file to write the values to plot in
    FILE *file; // File pointer
    char text[1000000] = "This is a sample text."; // Text to write to the file
    // Reset the text array to an empty string
    text[0] = '\0';


    double thrust;
    double vsound;
    double Mach;
    double b;
    double dz=500; //we will evaluate the rocket's behaviour between alt-dz and alt+dz, 500 is set arbitrarily       
    //for loops that simulates the rocket trajectory in function of the thrusts and the mass of the fuels
    for(double m_fuel2 = 0.1; m_fuel2 <= m_fuel2MAX; m_fuel2+=2500){

        for(double m_fuel1 = 0.1; m_fuel1 <= m_fuel1MAX; m_fuel1+=25000){  

            for(double FT2 = 0.1; FT2 <= FT2_MAX; FT2+=10000){

                for(double FT1 = 0.1; FT1 <= FT1_MAX; FT1+=100000){

                    double drag = 0;
                    double ratio = 1.4; //heat capacity ratio
                    double cd = cd0;
                    double v = v0;
                    double z = z0;
                    double m_tot = m_fuel1 + m_fuel2 + m1 + m2 + m_charge + m_finale; //total mass of the rocket (imaginons qu'on rajoute 8 tonnes de marchandises)
                    double m = m_tot;


                    //the t burns are in function of the thrusts and the mass of the fuels
                    int t_burn1;
                    if (FT1 != 0) {
                        t_burn1 = (m_fuel1*g0*I_sp1)/FT1; // tburn of stage 1
                    } else {
                        // handle the error, e.g., by setting t_burn1 to t_final
                        t_burn1 = t_final; // or any other value that indicates an error
                    }
                    
                    int t_burn2;
                    if (FT2 != 0) {
                        t_burn2 = (m_fuel2*g0*I_sp2)/FT2; // tburn of stage 1
                    } else {
                        // handle the error, e.g., by setting t_burn2 to t_final
                        t_burn2 = t_final; // or any other value that indicates an error
                    }
                    double t_burn_total = t_burn1 + t_burn2; //Pre-calculate t_burn1 + t_burn2 outside the loop to avoid calculating it in every iteration
                    double epsilon = 1e-10; // to adjust for floating point errors
                    
                    //iteration on time
                    for (double t = 0; t < t_final; t+=dt)
                    {
                        double g = gravity(z);
                        
                        //before the first fuel is finished
                        if (t < t_burn1 + epsilon){
                            m -= m_variable(dt, FT1, I_sp1);
                            thrust = F_thrust(t, t_burn1,FT1)/m;              //thrust acceleration 
                            if(fabs(t-t_burn1)<dt/2){                   //the case t == t_burn1 
                                m -= m1;
                            }
                        }
                        //before the second fuel is finished
                        else if ((t < t_burn_total + epsilon) && (t > t_burn1 + epsilon)){
                            m -= m_variable(dt, FT2, I_sp2);
                            thrust = F_thrust(t, t_burn_total,FT2)/m;              //thrust acceleration 
                            if(fabs(t-t_burn_total)<dt/2){
                                m -= m2;
                            }
                            
                        }
                        else{
                            thrust  = 0;
                            
                        }
                        
                        double d = density(z);
                        vsound = v_sound(ratio, temperature(z));
                        Mach = M_fct(v,vsound);
                        b = beta(Mach);
                        cd = c_d(cd0, b);
                        drag = F_drag(z,v,diameter,cd);
                        if(drag<0){         //flip drag force vector if rocket falls
                            drag = drag*(-1);
                        }

                        //change of speed and altitude
                        v = v + (thrust - drag - g)*dt;
                        z = z + v * dt;

                        if(z<0){            //if rocket crashes
                            //printf("crash du rocket...\n");
                            break;   
                        }

                        // we want to find the minimum time to reach 417km at a slow speed
                        if((z>alt_min_km-dz)&&(z<alt_min_km+dz)&&(fabs(v)<5)&&(t<t_min)){ 
                                t_min = t;
                                //printf("FT1  : %f newtons, ", FT1);
                                //printf("FT2  : %f newtons, ", FT2);
                                //printf("mass fuel 1 : %f kg, ", m_fuel1);
                                //printf("mass fuel 2 : %f kg, ", m_fuel2);
                                FT1_plot = FT1;
                                FT2_plot = FT2;
                                m_fuel1_plot = m_fuel1;
                                m_fuel2_plot = m_fuel2;
                                //printf("temps  : %f secondes, ", t);
                                //printf("altitude : %f meters, ", z);
                                //printf("mass : %f kg", m);
                                //printf(" speed : %f m/s\n", v);
                                //printf("t_min : %f\n", t_min);
                                break;
                        
                        }

                        
                    }    
                }
            }
        
        }
    }

//We should have done this differently for a better optimisation
    //for loop for plotting the values
    double v = v0;
    double z = z0;
    double m_tot_plot = m_fuel1_plot + m_fuel2_plot + m1 + m2 + m_charge + m_finale; //total mass of the rocket (imaginons qu'on rajoute 8 tonnes de marchandises)
    double m = m_tot_plot;
    double drag = 0;
    double ratio = 1.4; //heat capacity ratio
    double cd = cd0;

    int t_burn1;
    if (FT1_plot != 0) {
        t_burn1 = (m_fuel1_plot*g0*I_sp1)/FT1_plot; // tburn of stage 1
    } else {
        // handle the error, e.g., by setting t_burn1 to t_final
        t_burn1 = t_final; // or any other value that indicates an error
    }
    
    int t_burn2;
    if (FT2_plot != 0) {
        t_burn2 = (m_fuel2_plot*g0*I_sp2)/FT2_plot; // tburn of stage 1
    } else {
        // handle the error, e.g., by setting t_burn2 to t_final
        t_burn2 = t_final; // or any other value that indicates an error
    }
    double t_burn_total = t_burn1 + t_burn2; //Pre-calculate t_burn1 + t_burn2 outside the loop to avoid calculating it in every iteration
    double epsilon = 1e-10; // to adjust for floating point errors

    for (double t = 0; t < t_min + dt; t+=dt)
    {
        double g = gravity(z);
        if (t < t_burn1 + epsilon){
            m -= m_variable(dt, FT1_plot, I_sp1);
            thrust = F_thrust(t, t_burn1,FT1_plot)/m;              //thrust acceleration 
            if(fabs(t-t_burn1)<dt/2){                   //the case t == t_burn1 
                m -= m1;
            }
        }
        else if ((t < t_burn_total + epsilon) && (t > t_burn1 + epsilon)){
            m -= m_variable(dt, FT2_plot, I_sp2);
            thrust = F_thrust(t, t_burn_total,FT2_plot)/m;              //thrust acceleration 
            if(fabs(t-t_burn_total)<dt/2){
                m -= m2;
            }
            
        }
        else{
            thrust  = 0;
            
        }
        double d = density(z);
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
        //printf("temps  : %f secondes, ", t);
        //printf("altitude : %f meters, ", z);
        //printf("mass : %f kg", m);
        //printf(" speed : %f m/s\n", v);

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
    
    struct time launch_time_rocket = when_to_launch(hour_min, min_min, sec_min, t_min);
    printf("The launch time will need to be at %02d:%02d:%02d.\n", launch_time_rocket.hour, launch_time_rocket.min, launch_time_rocket.sec);

    //Print the final output on a separate text file
    
    // Open a file in write mode ('w')
    file = fopen("final_output.txt", "w");

    // Check if the file is successfully opened
    if (file == NULL) {
        printf("File could not be opened.");
        return 1; // Return an error code
    }

    // Write to the file
    fprintf(file, "Travel time (in seconds): %f\n", t_min);
    fprintf(file, "Launch time : %02d:%02d:%02d\n", launch_time_rocket.hour, launch_time_rocket.min, launch_time_rocket.sec);

    // Close the file
    fclose(file);

    // We asked ChatGPT3.5 again
    // Run Python script using system()
    int stats = system("python graphs.py");
    
    if (stats == -1) {
        printf("Failed to execute the command.\n");
        return 1;
    } else {
        printf("Python script executed successfully.\n");
    }
}

// The code takes a little time to run (around 1 minute)
