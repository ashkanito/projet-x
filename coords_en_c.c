#include <stdio.h>

//coordinates of ISS, read on a text file
struct Coordinates_min_dist {
    double distance;
    double lat;
    double lon;
    double alt;
    int hour;
    int min;
    int sec;
};

int readfile(char *filename, int length, struct Coordinates_min_dist *conc) { //NOTEPERSO: SOURCE: midterm d'entrainement
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

int main() {struct Coordinates_min_dist conc[50];
  double n = readfile("output_coords.txt", 1, conc);
  double dist_min = conc[0].distance;
  double lat_min = conc[0].lat;
  double lon_min = conc[0].lon;
  double alt_min = conc[0].alt;
  int hour_min = conc[0]. hour;
  int min_min = conc[0].min;
  int sec_min = conc[0].sec;
  return 0;
}
