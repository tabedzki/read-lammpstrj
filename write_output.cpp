#include "globals.h"
#include <sys/stat.h>


void write_vector_fields(bool separate_files){
    
  if (sites_per_LC == 0 ) return;

    if (myrank == 0) {
        std::string vec_dir = "vectors";
        FILE* otp;
        char file_path[100];
        int ind,  id;
        if (not separate_files) {
          if (step == 0) {
            otp = fopen("vectors.dat", "w");
            fprintf(otp, "step id x y ");
            if (Dim == 3) fprintf(otp, "z ");
            fprintf(otp, "dx dy ");
            if (Dim == 3) fprintf(otp, "dz ");
            fprintf(otp, "\n");
          } else {
            otp = fopen("vectors.dat", "a");
          }
        } else {
            mkdir(vec_dir.c_str(),  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            std::string file = "vector_step_" + std::to_string(step) + "_.dat";
            std::string path = vec_dir + "/" + file;
            // strcpy(file_path, vec_dir.c_str)); 
            otp = fopen(path.c_str(), "w");
            fprintf(otp, "step id x y ");
            if (Dim == 3) fprintf(otp, "z ");
            fprintf(otp, "dx dy ");
            if (Dim == 3) fprintf(otp, "dz ");
            fprintf(otp, "\n");

        }


        for (id = 0; id < nstot; id++){
            if (!LC_has_orientation(id))continue; 

            fprintf(otp, "%d ", step);
            fprintf(otp, "%d ", id);
            for (int  m=0 ; m<Dim ; m++ ) { fprintf(otp,"%lf ", x[id][m]); 
            }
              fprintf(otp,"%lf", mono_u[id][0]); 
            for (int  m=1 ; m<Dim ; m++ ) { fprintf(otp," %lf", mono_u[id][m]); 
            } 
            fprintf(otp,"\n");
        }
        fclose(otp);

    }  // if ( myrank == 0 )

    MPI_Barrier(MPI_COMM_WORLD);
}

void write_vector_fields(int local ){
    if (local==0) {
        write_vector_fields(false);
    } else {
        write_vector_fields(true);
    }
}

void write_vector_fields(void){
    write_vector_fields(multi_vec_files);
}